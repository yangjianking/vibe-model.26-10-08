# tests：单元测试与组件测试（Catch2 + 自检程序）

# 单元测试（Catch2）
add_executable(unit_tests
    unit/test_grid_metric.cpp
    unit/test_config.cpp
    unit/test_soa_precision.cpp
)
target_link_libraries(unit_tests PRIVATE
    cs_common cs_grid cs_dynamics cs_timeint cs_physics cs_io cs_parallel
    Catch2::Catch2WithMain)

add_test(NAME unit_tests COMMAND unit_tests)

# 数据同化模块自检（观测算子 TL/AD 转置、协方差语义、4D-Var 梯度）
add_executable(da_selfcheck
    da_selfcheck.cpp
)
target_link_libraries(da_selfcheck PRIVATE
    cs_common cs_grid cs_dynamics cs_io cs_da)

add_test(NAME da_selfcheck COMMAND da_selfcheck)

# 检验评估自检（守恒诊断 + 标准算例解析解）
add_executable(p7_selfcheck
    p7_selfcheck.cpp
)
target_link_libraries(p7_selfcheck PRIVATE
    cs_common cs_grid cs_dynamics cs_diagnostics cs_testcases)

add_test(NAME p7_selfcheck COMMAND p7_selfcheck)

# 平流通量配对守恒性实测
add_executable(advection_conservation
    advection_conservation.cpp
)
target_link_libraries(advection_conservation PRIVATE
    cs_common cs_grid cs_dynamics cs_diagnostics)

add_test(NAME advection_conservation COMMAND advection_conservation)

# 组件测试（后续阶段扩展：时间积分、halo 交换）
# add_subdirectory(integration)
