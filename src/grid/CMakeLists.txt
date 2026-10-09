# grid 模块：立方球网格（Cubed-Sphere Grid）
# 核心水平离散方案：gnomonic equiangular 投影 + Arakawa C/D 交错。

add_library(cs_grid OBJECT
    src/grid/cubed_sphere_grid.cpp
    src/grid/panel_rotation.cpp
    src/grid/vertical_coordinate.cpp
    src/grid/topography.cpp
)
target_link_libraries(cs_grid PUBLIC cs_common)

add_library(cubed_sph::grid ALIAS cs_grid)
