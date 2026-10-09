// =============================================================================
//  配置系统（P9：YAML 注入 + 分层覆盖 + fail-fast 校验）
// =============================================================================
//  所有可调参数（网格、时间步长、扩散系数、精度、并行布局、I/O）一律经
//  YAML 配置注入，禁止在业务代码中硬编码。配置加载顺序（后者覆盖前者）：
//    default.yaml  <  site.yaml  <  case.yaml  <  命令行 --set key=value
//
//  非法配置在初始化阶段 fail-fast 并明确指出字段（P9 第 1 条）。
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

#include <string>
#include <map>
#include <vector>
#include <optional>
#include <stdexcept>

namespace cubed_sph {

// ---------------------------------------------------------------------------
// 配置异常（非法配置 fail-fast）
// ---------------------------------------------------------------------------
class ConfigError : public std::runtime_error {
public:
    explicit ConfigError(const std::string& msg) : std::runtime_error(msg) {}
};

// ---------------------------------------------------------------------------
// 配置值：统一的可变类型容器（简化 YAML 节点的内部表示）
// ---------------------------------------------------------------------------
class Config {
public:
    // 从一组 YAML 文件加载并合并（后加载的覆盖先加载的）
    static Config load(const std::vector<std::string>& paths);

    // 命令行 --set key=value 覆盖（点号分隔层级，如 grid.nx=48）
    void apply_cmdline(const std::vector<std::string>& overrides);

    // 类型化读取（缺失或类型不符则抛出 ConfigError）
    int    get_int(const std::string& key) const;
    double get_double(const std::string& key) const;
    bool   get_bool(const std::string& key) const;
    std::string get_string(const std::string& key) const;

    // 带默认值的读取
    int    get_int(const std::string& key, int def) const noexcept;
    double get_double(const std::string& key, double def) const noexcept;
    bool   get_bool(const std::string& key, bool def) const noexcept;
    std::string get_string(const std::string& key, const std::string& def) const;

    // 嵌套子配置（用于网格/物理/同化等分组）
    const Config& sub(const std::string& key) const;
    bool has(const std::string& key) const;

    // 序列化当前生效配置（写入输出目录 provenance/ 供完全复现）
    std::string dump() const;

    // 计算配置内容的哈希（可复现性：配置哈希纳入 provenance）
    std::string content_hash() const;

private:
    // 内部扁平键值存储（key 已展开为点号路径）
    std::map<std::string, std::string> values_;
    std::map<std::string, Config> children_;

    void set_raw(const std::string& key, const std::string& value);
};

}  // namespace cubed_sph
