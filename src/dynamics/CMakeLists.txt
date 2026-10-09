# dynamics 模块：非静力动力核心

add_library(cs_dynamics OBJECT
    src/dynamics/state.cpp
    src/dynamics/equations.cpp
    src/dynamics/advection.cpp
    src/dynamics/reference_state.cpp
    src/dynamics/coriolis.cpp
)
target_link_libraries(cs_dynamics PUBLIC cs_common cs_grid)

add_library(cubed_sph::dynamics ALIAS cs_dynamics)
