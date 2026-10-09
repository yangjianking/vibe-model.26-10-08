# io 模块：数据输入输出（NetCDF 并行 I/O + 降级路径）

add_library(cs_io OBJECT
    src/io/io.cpp
)
target_link_libraries(cs_io PUBLIC cs_common cs_grid cs_dynamics)

if(CS_ENABLE_NETCDF AND NetCDF_FOUND)
    target_compile_definitions(cs_io PRIVATE CS_HAVE_NETCDF)
    target_link_libraries(cs_io PUBLIC NetCDF::NetCDF)
endif()

add_library(cubed_sph::io ALIAS cs_io)
