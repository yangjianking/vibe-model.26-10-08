# parallel 模块：MPI 并行环境、域分解与 halo 交换

add_library(cs_parallel OBJECT
    src/parallel/parallel.cpp
)
target_link_libraries(cs_parallel PUBLIC cs_common cs_grid)

add_library(cubed_sph::parallel ALIAS cs_parallel)
