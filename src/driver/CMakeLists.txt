# driver 模块：模式驱动与主程序

add_library(cs_driver OBJECT
    src/driver/driver.cpp
)
target_link_libraries(cs_driver PUBLIC
    cs_common cs_grid cs_dynamics cs_timeint cs_physics cs_io cs_parallel)

# 可执行文件
add_executable(cubed_sphere_nwp src/driver/main.cpp)
target_link_libraries(cubed_sphere_nwp PRIVATE
    cs_common cs_grid cs_dynamics cs_timeint cs_physics cs_io cs_parallel)
