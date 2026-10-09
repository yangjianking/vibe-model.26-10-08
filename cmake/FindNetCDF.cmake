# =============================================================================
#  NetCDF-4/HDF5 查找模块（P4 并行 I/O 主格式）
# =============================================================================
# 说明：探测系统 NetCDF-C 库（可选 MPI 并行版本），设置 NetCDF_FOUND、
#       NetCDF_INCLUDE_DIRS、NetCDF_LIBRARIES。若未找到则优雅降级：
#       I/O 子系统退化为 ASCII/二进制占位实现，模式其余部分不受影响。

find_path(NetCDF_INCLUDE_DIR
    NAMES netcdf.h
    HINTS ENV NETCDF_DIR
    PATH_SUFFIXES include)

find_library(NetCDF_LIBRARY
    NAMES netcdf
    HINTS ENV NETCDF_DIR
    PATH_SUFFIXES lib lib64)

set(NetCDF_INCLUDE_DIRS ${NetCDF_INCLUDE_DIR})
set(NetCDF_LIBRARIES ${NetCDF_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
    REQUIRED_VARS NetCDF_LIBRARY NetCDF_INCLUDE_DIR)

if(NetCDF_FOUND)
    add_library(NetCDF::NetCDF UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF PROPERTIES
        IMPORTED_LOCATION "${NetCDF_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIR}")
endif()
