# diagnostics 模块：守恒诊断

add_library(cs_diagnostics OBJECT
    src/diagnostics/diagnostics.cpp
)
target_link_libraries(cs_diagnostics PUBLIC cs_common cs_grid cs_dynamics)

add_library(cubed_sph::diagnostics ALIAS cs_diagnostics)
