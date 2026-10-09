# physics 模块：物理参数化接口与占位方案

add_library(cs_physics OBJECT
    src/physics/physics.cpp
)
target_link_libraries(cs_physics PUBLIC cs_common cs_grid cs_dynamics)

add_library(cubed_sph::physics ALIAS cs_physics)
