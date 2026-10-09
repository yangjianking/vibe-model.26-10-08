// =============================================================================
//  配置系统实现（YAML 解析 + 分层覆盖 + fail-fast 校验）
// =============================================================================

#include "cubed_sph/common/config.hpp"

#include <yaml-cpp/yaml.h>

#include <fstream>
#include <sstream>
#include <functional>
#include <iomanip>

namespace cubed_sph {

namespace {

// 将 YAML 节点递归展开为扁平 key=value 与子 Config 树
void flatten(const std::string& prefix, const YAML::Node& node,
             Config& out) {
    if (node.IsScalar()) {
        out.set_raw(prefix, node.as<std::string>());
    } else if (node.IsMap()) {
        for (const auto& kv : node) {
            const std::string key = kv.first.as<std::string>();
            const std::string full = prefix.empty() ? key : prefix + "." + key;
            flatten(full, kv.second, out);
        }
    } else if (node.IsSequence()) {
        // 序列暂以逗号拼接存储（本阶段 I/O/观测清单等简单列表够用）
        std::ostringstream oss;
        bool first = true;
        for (const auto& item : node) {
            if (!first) oss << ",";
            oss << item.as<std::string>();
            first = false;
        }
        out.set_raw(prefix, oss.str());
    }
    // IsNull 节点跳过
}

// 简单的字符串哈希（FNV-1a，用于配置内容指纹）
std::string fnv1a(const std::string& s) {
    std::uint64_t h = 1469598103934665603ULL;
    for (char c : s) {
        h ^= static_cast<unsigned char>(c);
        h *= 1099511628211ULL;
    }
    std::ostringstream oss;
    oss << std::hex << std::setw(16) << std::setfill('0') << h;
    return oss.str();
}

}  // namespace

Config Config::load(const std::vector<std::string>& paths) {
    Config cfg;
    for (const auto& p : paths) {
        std::ifstream ifs(p);
        if (!ifs.good()) {
            throw ConfigError("无法打开配置文件: " + p);
        }
        try {
            YAML::Node root = YAML::Load(ifs);
            flatten("", root, cfg);
        } catch (const YAML::Exception& e) {
            throw ConfigError("YAML 解析失败 [" + p + "]: " + e.what());
        }
    }
    return cfg;
}

void Config::apply_cmdline(const std::vector<std::string>& overrides) {
    for (const auto& o : overrides) {
        const auto eq = o.find('=');
        if (eq == std::string::npos) {
            throw ConfigError("命令行覆盖格式应为 key=value: " + o);
        }
        const std::string key = o.substr(0, eq);
        const std::string val = o.substr(eq + 1);
        set_raw(key, val);
    }
}

void Config::set_raw(const std::string& key, const std::string& value) {
    values_[key] = value;
}

int Config::get_int(const std::string& key) const {
    auto it = values_.find(key);
    if (it == values_.end()) throw ConfigError("缺少配置项: " + key);
    try { return std::stoi(it->second); }
    catch (...) { throw ConfigError("配置项非整数: " + key + "=" + it->second); }
}

double Config::get_double(const std::string& key) const {
    auto it = values_.find(key);
    if (it == values_.end()) throw ConfigError("缺少配置项: " + key);
    try { return std::stod(it->second); }
    catch (...) { throw ConfigError("配置项非浮点: " + key + "=" + it->second); }
}

bool Config::get_bool(const std::string& key) const {
    auto it = values_.find(key);
    if (it == values_.end()) throw ConfigError("缺少配置项: " + key);
    const std::string& v = it->second;
    if (v == "true" || v == "True" || v == "1") return true;
    if (v == "false" || v == "False" || v == "0") return false;
    throw ConfigError("配置项非布尔: " + key + "=" + v);
}

std::string Config::get_string(const std::string& key) const {
    auto it = values_.find(key);
    if (it == values_.end()) throw ConfigError("缺少配置项: " + key);
    return it->second;
}

int Config::get_int(const std::string& key, int def) const noexcept {
    try { return get_int(key); } catch (...) { return def; }
}
double Config::get_double(const std::string& key, double def) const noexcept {
    try { return get_double(key); } catch (...) { return def; }
}
bool Config::get_bool(const std::string& key, bool def) const noexcept {
    try { return get_bool(key); } catch (...) { return def; }
}
std::string Config::get_string(const std::string& key,
                               const std::string& def) const {
    auto it = values_.find(key);
    return it == values_.end() ? def : it->second;
}

const Config& Config::sub(const std::string& key) const {
    // 子配置通过前缀匹配实现：返回一个共享底层存储的轻量视图
    // 为简化本阶段实现，直接返回 *this（读取时用 "key." 前缀即可）。
    static Config empty;
    return *this;
}

bool Config::has(const std::string& key) const {
    return values_.count(key) > 0;
}

std::string Config::dump() const {
    std::ostringstream oss;
    for (const auto& [k, v] : values_) {
        oss << k << " = " << v << "\n";
    }
    return oss.str();
}

std::string Config::content_hash() const {
    // 按 key 排序保证哈希确定性
    std::string canon;
    for (const auto& [k, v] : values_) {
        canon += k + "=" + v + ";";
    }
    return fnv1a(canon);
}

}  // namespace cubed_sph
