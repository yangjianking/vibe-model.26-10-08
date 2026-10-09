# testcases 模块：标准算例（Williamson / Jablonowski / Held-Suarez）

add_library(cs_testcases OBJECT
    src/testcases/test_cases.cpp
)
target_link_libraries(cs_testcases PUBLIC cs_common cs_grid cs_dynamics)

add_library(cubed_sph::testcases ALIAS cs_testcases)
