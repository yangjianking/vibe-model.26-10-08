// =============================================================================
//  配置系统单元测试
// =============================================================================
//  验证 YAML 解析、分层覆盖、fail-fast 校验与配置哈希确定性。
// =============================================================================

#include <catch2/catch_test_macros.hpp>

#include "cubed_sph/common/config.hpp"

using namespace cubed_sph;

TEST_CASE("配置缺失项抛出 ConfigError（fail-fast）", "[config]") {
    Config cfg = Config::load({"configs/default.yaml"});
    REQUIRE(cfg.get_int("grid.ncells") == 48);

    // 缺失项应抛出异常
    REQUIRE_THROWS_AS(cfg.get_int("nonexistent.key"), ConfigError);
}

TEST_CASE("带默认值读取回退", "[config]") {
    Config cfg = Config::load({"configs/default.yaml"});
    REQUIRE(cfg.get_int("grid.ncells", 99) == 48);
    REQUIRE(cfg.get_int("no.such.key", 99) == 99);
    REQUIRE(cfg.get_double("no.such.key", 1.5) == 1.5);
    REQUIRE(cfg.get_bool("no.such.key", true) == true);
}

TEST_CASE("命令行覆盖生效", "[config]") {
    Config cfg = Config::load({"configs/default.yaml"});
    cfg.apply_cmdline({"grid.ncells=96", "timeint.dt=150.0"});
    REQUIRE(cfg.get_int("grid.ncells") == 96);
    REQUIRE(cfg.get_double("timeint.dt") == Catch::Approx(150.0));
}

TEST_CASE("配置哈希确定性", "[config]") {
    Config a = Config::load({"configs/default.yaml"});
    Config b = Config::load({"configs/default.yaml"});
    REQUIRE(a.content_hash() == b.content_hash());
    REQUIRE(!a.content_hash().empty());
}
