# timeint 模块：时间积分（半隐式 + 显式 RK3）

add_library(cs_timeint OBJECT
    src/timeint/time_integrator.cpp
    src/timeint/bicgstab.cpp
    src/timeint/helmholtz_operator.cpp
)
target_link_libraries(cs_timeint PUBLIC cs_common cs_grid cs_dynamics)

add_library(cubed_sph::timeint ALIAS cs_timeint)
