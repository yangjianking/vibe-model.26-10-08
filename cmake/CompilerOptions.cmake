# =============================================================================
#  编译器选项统一管理
# =============================================================================
# 说明：所有编译警告与优化级别集中于此，避免各子目录重复设置。
#       零容忍警告遵循提示词 P9 硬性要求（-Wall -Wextra）。

# 优化级别：Release 走 -O3，Debug 保留符号
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-O3)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        add_compile_options(-O3)
    endif()
endif()

# 位置无关代码（便于链接共享库）
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 生成依赖文件
set(CMAKE_DEPFILE_FLAGS_CXX "-MMD -MT <OBJECT> -MF <DEPFILE>")
