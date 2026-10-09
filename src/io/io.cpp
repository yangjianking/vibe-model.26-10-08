// =============================================================================
//  数据 I/O 实现（NetCDF 输出 + 重启 + 观测空间）
// =============================================================================

#include "cubed_sph/io/io.hpp"

#include <fstream>
#include <cmath>
#include <cstring>
#include <algorithm>

#if defined(CS_HAVE_NETCDF)
#include <netcdf.h>
#endif

namespace cubed_sph::io {

// =============================================================================
// NetCDFWriter
// =============================================================================
NetCDFWriter::NetCDFWriter(const Config& cfg) {
#if defined(CS_HAVE_NETCDF)
    netcdf_available_ = true;
#else
    netcdf_available_ = false;
#endif
    (void)cfg;
}

NetCDFWriter::~NetCDFWriter() { close(); }

void NetCDFWriter::open(const std::string& path, bool parallel) {
    path_ = path;
    opened_ = true;
    (void)parallel;
#if defined(CS_HAVE_NETCDF)
    // TODO(io): nc_create + 定义维度（panel, n, nlev, time）+ UGRID 元数据。
    // 完整 NetCDF-4 并行（HDF5 collective）见后续阶段，此处探测可用性。
#endif
}

void NetCDFWriter::write_fields(const std::vector<OutputField>& fields,
                                StateReal time) {
    if (!opened_) return;
#if defined(CS_HAVE_NETCDF)
    // TODO(io): 变量写入 nc_put_vara_*（并行 collective）。
    (void)fields; (void)time;
#else
    // 降级：二进制占位输出（调试用，非生产格式）
    std::ofstream ofs(path_ + ".bin", std::ios::app | std::ios::binary);
    ofs.write(reinterpret_cast<const char*>(&time), sizeof(StateReal));
    for (const auto& f : fields) {
        const size_t n = f.data ? f.data->size() : 0;
        ofs.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
        if (f.data && n > 0) {
            ofs.write(reinterpret_cast<const char*>(f.data->data()),
                      n * sizeof(StateReal));
        }
    }
#endif
}

void NetCDFWriter::write_grid_metadata(const dynamics::State& state) {
    // UGRID 立方球网格元数据：面板拓扑、度量、坐标。
    // 完整实现写出 CF-1.10 + UGRID 约定的网格变量（见 docs/theory/io.md）。
    (void)state;
}

void NetCDFWriter::close() { opened_ = false; }

// =============================================================================
// RestartIO
// =============================================================================
namespace {

// 序列化状态字段为扁平缓冲（顺序固定，保证位可重复）
void serialize_state(const dynamics::State& state,
                     std::vector<StateReal>& buf) {
    const IIndex total = state.ncell();
    buf.clear();
    buf.reserve(static_cast<size_t>(total) * 8);
    buf.insert(buf.end(), state.rho.begin(), state.rho.end());
    buf.insert(buf.end(), state.rho_u.begin(), state.rho_u.end());
    buf.insert(buf.end(), state.rho_v.begin(), state.rho_v.end());
    buf.insert(buf.end(), state.rho_w.begin(), state.rho_w.end());
    buf.insert(buf.end(), state.rho_theta.begin(), state.rho_theta.end());
    buf.insert(buf.end(), state.exner.begin(), state.exner.end());
    buf.insert(buf.end(), state.theta.begin(), state.theta.end());
    buf.insert(buf.end(), state.geopotential.begin(), state.geopotential.end());
}

void deserialize_state(dynamics::State& state,
                       const std::vector<StateReal>& buf) {
    const IIndex total = state.ncell();
    if (buf.size() != static_cast<size_t>(total) * 8) {
        return;  // 尺寸不符，交由调用方报错
    }
    size_t off = 0;
    auto copy = [&](std::vector<StateReal>& dst) {
        std::copy(buf.begin() + off, buf.begin() + off + total, dst.begin());
        off += total;
    };
    copy(state.rho);
    copy(state.rho_u);
    copy(state.rho_v);
    copy(state.rho_w);
    copy(state.rho_theta);
    copy(state.exner);
    copy(state.theta);
    copy(state.geopotential);
}

}  // namespace

// CRC32C（Castagnoli 多项式，硬件加速友好）
std::uint32_t RestartIO::crc32c(const std::uint8_t* data, std::size_t len) {
    // 软件实现（SSE4.2 硬件加速在 HPC 环境用 _mm_crc32_u64，此处用逐位）
    std::uint32_t crc = 0xFFFFFFFFu;
    for (std::size_t i = 0; i < len; ++i) {
        crc ^= data[i];
        for (int b = 0; b < 8; ++b) {
            const std::uint32_t mask =
                static_cast<std::uint32_t>(-(static_cast<std::int32_t>(crc & 1)));
            crc = (crc >> 1) ^ (0x82F63B78u & mask);  // Castagnoli 反射多项式
        }
    }
    return ~crc;
}

void RestartIO::serialize(const dynamics::State& state,
                          std::vector<StateReal>& buf) {
    serialize_state(state, buf);
}

void RestartIO::deserialize(dynamics::State& state,
                            const std::vector<StateReal>& buf) {
    deserialize_state(state, buf);
}

void RestartIO::write(const dynamics::State& state,
                      const RestartMetadata& meta,
                      const std::string& path) {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.good()) return;

    // 1. 写元数据头
    ofs.write(reinterpret_cast<const char*>(&meta.magic), sizeof(meta.magic));
    ofs.write(reinterpret_cast<const char*>(&meta.version), sizeof(meta.version));
    ofs.write(reinterpret_cast<const char*>(&meta.ncells), sizeof(meta.ncells));
    ofs.write(reinterpret_cast<const char*>(&meta.nhalo), sizeof(meta.nhalo));
    ofs.write(reinterpret_cast<const char*>(&meta.nlev), sizeof(meta.nlev));
    ofs.write(reinterpret_cast<const char*>(&meta.time), sizeof(meta.time));
    ofs.write(reinterpret_cast<const char*>(&meta.step), sizeof(meta.step));

    // 配置哈希（字符串，长度 + 内容）
    const std::uint32_t ch_len = static_cast<std::uint32_t>(meta.config_hash.size());
    ofs.write(reinterpret_cast<const char*>(&ch_len), sizeof(ch_len));
    ofs.write(meta.config_hash.data(), ch_len);

    // 2. 序列化状态字段
    std::vector<StateReal> buf;
    serialize(state, buf);
    ofs.write(reinterpret_cast<const char*>(buf.data()),
              buf.size() * sizeof(StateReal));

    // 3. CRC32C 校验（覆盖状态字段）
    const std::uint32_t crc = crc32c(
        reinterpret_cast<const std::uint8_t*>(buf.data()),
        buf.size() * sizeof(StateReal));
    ofs.write(reinterpret_cast<const char*>(&crc), sizeof(crc));
}

RestartMetadata RestartIO::read(dynamics::State& state,
                                const std::string& path) {
    RestartMetadata meta;
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.good()) {
        meta.magic = 0;  // 标记读取失败
        return meta;
    }

    // 1. 读元数据头
    ifs.read(reinterpret_cast<char*>(&meta.magic), sizeof(meta.magic));
    if (meta.magic != 0x43535253) {  // 魔数校验
        meta.magic = 0;
        return meta;
    }
    ifs.read(reinterpret_cast<char*>(&meta.version), sizeof(meta.version));
    ifs.read(reinterpret_cast<char*>(&meta.ncells), sizeof(meta.ncells));
    ifs.read(reinterpret_cast<char*>(&meta.nhalo), sizeof(meta.nhalo));
    ifs.read(reinterpret_cast<char*>(&meta.nlev), sizeof(meta.nlev));
    ifs.read(reinterpret_cast<char*>(&meta.time), sizeof(meta.time));
    ifs.read(reinterpret_cast<char*>(&meta.step), sizeof(meta.step));

    std::uint32_t ch_len = 0;
    ifs.read(reinterpret_cast<char*>(&ch_len), sizeof(ch_len));
    meta.config_hash.resize(ch_len);
    ifs.read(&meta.config_hash[0], ch_len);

    // 2. 读状态字段
    std::vector<StateReal> buf(static_cast<size_t>(state.ncell()) * 8);
    ifs.read(reinterpret_cast<char*>(buf.data()),
             buf.size() * sizeof(StateReal));
    deserialize(state, buf);

    // 3. CRC 校验
    std::uint32_t crc_read = 0;
    ifs.read(reinterpret_cast<char*>(&crc_read), sizeof(crc_read));
    const std::uint32_t crc_calc = crc32c(
        reinterpret_cast<const std::uint8_t*>(buf.data()),
        buf.size() * sizeof(StateReal));
    if (crc_read != crc_calc) {
        meta.magic = 0;  // CRC 失败，标记损坏
    }
    return meta;
}

// =============================================================================
// ObsSpace
// =============================================================================
void ObsSpace::add(const Observation& obs) {
    obs_.push_back(obs);
}

std::vector<const Observation*> ObsSpace::of_type(ObsType t) const {
    std::vector<const Observation*> out;
    for (const auto& o : obs_) {
        if (o.type == static_cast<int>(t)) out.push_back(&o);
    }
    return out;
}

std::vector<std::vector<size_t>> ObsSpace::time_slots(StateReal slot_width,
                                                      StateReal window_len) const {
    const size_t n_slots = static_cast<size_t>(window_len / slot_width) + 1;
    std::vector<std::vector<size_t>> slots(n_slots);
    for (size_t i = 0; i < obs_.size(); ++i) {
        const size_t slot = static_cast<size_t>(obs_[i].time / slot_width);
        if (slot < n_slots) slots[slot].push_back(i);
    }
    return slots;
}

void ObsSpace::gross_error_check(StateReal k) {
    // 背景检验（gross-error）：观测值超出合理范围即拒绝。
    // 完整实现需背景场 H(x_b)；此处用观测值自身的统计离群检验占位。
    if (obs_.empty()) return;
    double mean = 0.0;
    for (const auto& o : obs_) mean += o.value;
    mean /= static_cast<double>(obs_.size());
    double var = 0.0;
    for (const auto& o : obs_) {
        var += (o.value - mean) * (o.value - mean);
    }
    var /= static_cast<double>(obs_.size());
    const double std = std::sqrt(var);
    for (auto& o : obs_) {
        if (std::abs(o.value - mean) > k * std) {
            o.qc_flag = 1;  // 标记拒绝
        }
    }
}

}  // namespace cubed_sph::io
