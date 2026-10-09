# da 模块：数据同化（4D-Var / 切线性 / 伴随）

add_library(cs_da OBJECT
    src/da/observation_operator.cpp
    src/da/tlm_ad_utils.cpp
    src/da/tlm_ad.cpp
    src/da/four_dvar.cpp
)
target_link_libraries(cs_da PUBLIC cs_common cs_grid cs_dynamics cs_timeint cs_io)

add_library(cubed_sph::da ALIAS cs_da)
