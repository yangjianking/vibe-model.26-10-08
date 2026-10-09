// =============================================================================
//  数据 I/O 子系统接口（P4 并行 I/O）
// =============================================================================
//  模式数据进出（src/io/）。主格式 NetCDF-4/HDF5（CF-1.10 约定，UGRID
//  网格元数据），支持并行写出与重启。
//
//  职责划分：
//    - OutputWriter  ：预报场输出（NetCDF-4 并行，或退化为 ASCII/二进制）
//    - RestartIO     ：检查点/重启（自描述元数据 + CRC32C 校验 + N 进 M 出）
//    - ObsSpace      ：观测数据接口（与同化联动）
//
//  参考文献：
//    CF Conventions 1.10（气候预报元数据约定）
//    Klöwer et al. (2021), Nat. Comput. Sci. 1, 775–779（数据压缩）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/common/config.hpp"

#include <string>
#include <vector>
#include <cstdint>

namespace cubed_sph::io {

// ---------------------------------------------------------------------------
// 输出字段描述
// ---------------------------------------------------------------------------
struct OutputField {
    std::string name;         // 变量名（CF 标准名优先）
    std::string long_name;    // 描述
    std::string units;        // 单位
    std::vector<StateReal>* data = nullptr;  // 指向 State 中字段的指针
};

// ---------------------------------------------------------------------------
// 预报输出接口
// ---------------------------------------------------------------------------
class OutputWriter {
public:
    // 打开输出文件（路径、是否并行、压缩级别）
    virtual void open(const std::string& path, bool parallel = false) = 0;

    // 写出一个时间步的字段集
    virtual void write_fields(const std::vector<OutputField>& fields,
                              StateReal time) = 0;

    // 写出网格元数据（UGRID 立方球拓扑）
    virtual void write_grid_metadata(const dynamics::State& state) = 0;

    virtual void close() = 0;
    virtual ~OutputWriter() = default;
};

// ---------------------------------------------------------------------------
// NetCDF 输出（若 NetCDF 可用；否则降级为 ASCII/二进制占位）
// ---------------------------------------------------------------------------
class NetCDFWriter : public OutputWriter {
public:
    explicit NetCDFWriter(const Config& cfg);
    ~NetCDFWriter() override;

    void open(const std::string& path, bool parallel) override;
    void write_fields(const std::vector<OutputField>& fields,
                      StateReal time) override;
    void write_grid_metadata(const dynamics::State& state) override;
    void close() override;

    // 是否使用 NetCDF 后端（否则降级为二进制占位）
    bool netcdf_available() const { return netcdf_available_; }

private:
    std::string path_;
    bool netcdf_available_ = false;
    bool opened_ = false;
};

// ---------------------------------------------------------------------------
// 重启文件头部元数据（自描述，版本化）
// ---------------------------------------------------------------------------
// 含恢复精确积分所需的全部信息（提示词 P4 第 3 条）：
//   - 网格/坐标/halo/时间/配置哈希
//   - 随机种子、迭代计数
//   - 版本号（格式演进）
// ---------------------------------------------------------------------------
struct RestartMetadata {
    std::uint32_t magic = 0x43535253;   // "CSRS"（Cubed-Sphere ReStart）
    std::uint32_t version = 1;           // 重启格式版本
    int ncells = 0;
    int nhalo = 0;
    int nlev = 0;
    StateReal time = 0.0;
    IIndex step = 0;
    std::string config_hash;             // 配置内容哈希（可复现）
    std::string compiler_version;        // 编译器版本
};

// ---------------------------------------------------------------------------
// 重启 / 检查点接口
// ---------------------------------------------------------------------------
class RestartIO {
public:
    // 写重启文件（含元数据 + 全部状态字段 + CRC32C 校验）
    void write(const dynamics::State& state,
               const RestartMetadata& meta,
               const std::string& path);

    // 读重启文件（支持变更 rank 数；返回元数据）
    RestartMetadata read(dynamics::State& state, const std::string& path);

    // CRC32C 校验（断点续算、损坏检测）
    static std::uint32_t crc32c(const std::uint8_t* data, std::size_t len);

private:
    // 序列化/反序列化状态字段（二进制，进程无关分块）
    void serialize(const dynamics::State& state, std::vector<StateReal>& buf);
    void deserialize(dynamics::State& state, const std::vector<StateReal>& buf);
};

// ---------------------------------------------------------------------------
// 观测数据接口（与 P5 同化联动）
// ---------------------------------------------------------------------------
// 统一抽象为 ObsSpace：观测向量 + 元数据 + 质控标记。支持观测在时间窗内的
// 时隙（time slot）分箱。
// ---------------------------------------------------------------------------
struct Observation {
    StateReal value = 0.0;       // 观测值
    StateReal error = 1.0;       // 观测误差标准差 σ_o
    StateReal lon = 0.0;         // 经度 [rad]
    StateReal lat = 0.0;         // 纬度 [rad]
    StateReal height = 0.0;      // 高度 [m]（探空为观测高度）
    StateReal time = 0.0;        // 观测时间（相对同化窗起点）[s]
    int type = 0;                // 观测类型（TEMP/SYNOP/...，枚举见下）
    int qc_flag = 0;             // 质控标记（0=通过，非 0=拒绝）
};

// 观测类型枚举（提示词 P5 第 4 条首批观测类型）
enum class ObsType : int {
    TEMP = 0,      // 探空温度
    SYNOP = 1,     // 地面报
    Satellite = 2, // 卫星晴空辐射率
    GNSS_RO = 3,   // GPS 无线电掩星折射率
    Scatterometer = 4, // 散射计海面风
    AMV = 5,       // 大气运动矢量
};

// ---------------------------------------------------------------------------
// 观测空间：观测向量容器 + 时隙分箱
// ---------------------------------------------------------------------------
class ObsSpace {
public:
    // 添加一个观测
    void add(const Observation& obs);

    // 按观测类型筛选
    std::vector<const Observation*> of_type(ObsType t) const;

    // 按时间窗时隙分箱（slot_width 为时隙宽度 [s]）
    // 返回每个时隙内的观测索引列表
    std::vector<std::vector<size_t>> time_slots(StateReal slot_width,
                                                StateReal window_len) const;

    // 基础质控：背景检验（gross-error check）
    // 超出 [bg - k*σ, bg + k*σ] 的观测被标记为拒绝
    void gross_error_check(StateReal k = 5.0);

    size_t size() const { return obs_.size(); }
    const Observation& operator[](size_t i) const { return obs_[i]; }

private:
    std::vector<Observation> obs_;
};

}  // namespace cubed_sph::io
