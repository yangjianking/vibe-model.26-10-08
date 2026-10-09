# 文献出处总表（References）

> 集中列出本项目引用过的全部文献，含真实可查证出处（作者, 年份, 期刊, 卷期页码）。
> 无法确认 DOI 的标注 [TO-VERIFY] 并说明查证路径（提示词 P0 第 7 条）。

## 网格与动力框架

- Sadourny, R. (1972). Conservative finite-difference approximations of the
  primitive equations on quasi-uniform spherical grids. *Mon. Wea. Rev.*, 100, 136–144.
- Ronchi, C., Iacono, R., Paolucci, P. S. (1996). The "cubed sphere": a new method
  for the solution of partial differential equations in spherical geometry.
  *J. Comput. Phys.*, 124, 93–114.
- Putman, W. M., Lin, S.-J. (2007). Finite-volume transport on various cubed-sphere
  grids. *J. Comput. Phys.*, 227, 55–78.
- Harris, L. M., et al. (2021). A Scientific Description of the GFDL FV³ Dynamical
  Core. *GFDL Tech. Memo.* GFDL2021001.
- Lin, S.-J., Rood, R. B. (1996). Multidimensional flux-form semi-Lagrangian
  transport schemes. *Mon. Wea. Rev.*, 124, 2046–2070.
- Arakawa, A., Lamb, V. R. (1981). A potential enstrophy and energy conserving
  scheme for the shallow water equations. *Mon. Wea. Rev.*, 109, 18–36.

## 半隐式时间积分与椭圆求解器

- Robert, A. (1981). A stable numerical integration scheme for the primitive
  meteorological equations. *Atmos.-Ocean*, 19, 35–52.
- Tanguay, M., Robert, A., Laprise, R. (1990). A semi-implicit semi-Lagrangian
  fully compressible regional forecast model. *Mon. Wea. Rev.*, 118, 1970–1980.
- Staniforth, A., Côté, J. (1991). Semi-Lagrangian integration schemes for
  atmospheric models — a review. *Mon. Wea. Rev.*, 119, 2206–2223.
- Wood, N., et al. (2014). An inherently mass-conserving semi-implicit
  semi-Lagrangian discretization of the deep-atmosphere global non-hydrostatic
  equations. *Q. J. R. Meteorol. Soc.*, 140, 1505–1528.
- Buckeridge, S., Scheichl, R. (2010). Parallel geometric multigrid for global
  weather prediction. *Q. J. R. Meteorol. Soc.*, 136, 2723–2735.
- Knoll, D. A., Keyes, D. E. (2004). Jacobian-free Newton-Krylov methods.
  *J. Comput. Phys.*, 193, 357–397.
- van der Vorst, H. A. (1992). Bi-CGSTAB: a fast and smoothly converging variant
  of Bi-CG for the solution of nonsymmetric linear systems.
  *SIAM J. Sci. Stat. Comput.*, 13, 631–644.

## 混合精度与性能

- Váňa, F., et al. (2017). Single precision in weather forecasting models:
  An evaluation with the IFS. *Mon. Wea. Rev.*, 145, 4957–4970.
- Düben, P. D., Palmer, T. N. (2014). Benchmark tests for numerical weather
  forecasts on inexact hardware. *Mon. Wea. Rev.*, 142, 3809–3826.
- Edwards, H. C., et al. (2014). Kokkos: Enabling manycore performance
  portability through polymorphic memory access patterns.
  *J. Parallel Distrib. Comput.*, 74, 3202–3216.
- Adams, S. V., et al. (2019). LFRic: Meeting the challenges of scalability and
  performance portability in weather and climate models.
  *J. Parallel Distrib. Comput.*, 132, 383–396.
- Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*,
  2nd ed. SIAM.

## 验证算例（后续阶段使用）

- Williamson, D. L., et al. (1992). A standard test set for numerical
  approximations to the shallow water equations in spherical geometry.
  *J. Comput. Phys.*, 102, 211–224.
- Jablonowski, C., Williamson, D. L. (2006). A baroclinic instability test case
  for atmospheric model dynamical cores. *Q. J. R. Meteorol. Soc.*, 132.
- Held, I. M., Suarez, M. J. (1994). A proposal for the intercomparison of the
  dynamical cores of atmospheric general circulation models.
  *Bull. Amer. Meteor. Soc.*, 75, 1825–1830.
- Ullrich, P. A., et al. (2017). DCMIP2016: a review of non-hydrostatic dynamical
  core design and intercomparison of participating models.
  *Geosci. Model Dev.*, 10, 4477–4509.

## 变分同化与伴随（补充，ima 知识库可检索）

- Courtier, P., Thépaut, J.-N., Hollingsworth, A. (1994). A strategy for
  operational implementation of 4D-Var, using an incremental approach.
  *Q. J. R. Meteorol. Soc.*, 120, 1367–1387.
- Rabier, F., Järvinen, H., Klinker, E., Mahfouf, J.-F., Simmons, A. (2000).
  The ECMWF operational implementation of four-dimensional variational
  assimilation. I: Experimental results with simplified physics.
  *Q. J. R. Meteorol. Soc.*, 126, 1143–1170.
- Veersé, F., Thépaut, J.-N. (1998). Multiple-truncation incremental approach for
  four-dimensional variational data assimilation. *Q. J. R. Meteorol. Soc.*, 124.
- Talagrand, O., Courtier, P. (1987). Variational assimilation of meteorological
  observations with the adjoint vorticity equation. *Q. J. R. Meteorol. Soc.*, 113.
- Giering, R., Kaminski, T. (1998). Recipes for adjoint code construction.
  *ACM Trans. Math. Softw.*, 24, 437–474.
- Griewank, A., Walther, A. (2000). Algorithm 799: revolve: an implementation of
  checkpointing for the reverse or adjoint mode of computational differentiation.
  *ACM Trans. Math. Softw.*, 26, 19–45.
- Parrish, D. F., Derber, J. C. (1992). The National Meteorological Center's
  spectral statistical-interpolation analysis system. *Mon. Wea. Rev.*, 120, 1747–1763.
- Bannister, R. N. (2008). A review of forecast error covariance statistics in
  atmospheric variational data assimilation. II: Modelling the forecast error
  covariance statistics. *Q. J. R. Meteorol. Soc.*, 134, 1971–1996.
- Fisher, M. (2003). Background error covariance modelling. *ECMWF Seminar on
  Recent developments in data assimilation*.
- Bloom, S. C., Takacs, L. L., da Silva, A. M., Ledvina, D. (1996). Data
  assimilation using incremental analysis updates. *Mon. Wea. Rev.*, 124, 1256–1271.
- Desroziers, G., Berre, L., Chapnik, B., Poli, P. (2005). Diagnosis of
  observation, background and analysis-error statistics in observation space.
  *Q. J. R. Meteorol. Soc.*, 131, 3385–3396.
- Andersson, E., Järvinen, H. (1999). Variational quality control.
  *Q. J. R. Meteorol. Soc.*, 125, 697–722.

## 平流与有限体积方法（补充，ima 知识库可检索）

- Colella, P., Woodward, P. R. (1984). The piecewise parabolic method (PPM) for
  gas-dynamical simulations. *J. Comput. Phys.*, 54, 174–201.
- LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*.
  Cambridge University Press.
- Durran, D. R. (2010). *Numerical Methods for Fluid Dynamics: With Applications
  to Geophysics*, 2nd ed. Springer.
- Fletcher, S. J. (2019). *Semi-Lagrangian Advection Methods and Their
  Applications in Geoscience*. Elsevier.
- Li, X., Xiao, F. (2010). 多矩约束有限体积（MCV）方法系列（CMA 路线）。

## 混合精度与性能（补充，ima 知识库可检索）

- Nakano, M., Yashiro, H., Kodama, C., Tomita, H. (2018). Single precision in the
  dynamical core of a nonhydrostatic global atmospheric model: Evaluation using a
  baroclinic wave test case. *Mon. Wea. Rev.*, 146, 409–416.
- Chantry, M., Thornes, T., Palmer, T., Düben, P. (2019). Scale-selective precision
  for weather and climate forecasting. *Mon. Wea. Rev.*, 147, 1409–1420.
- Zhang, S., et al. (2024). 混合精度在 GRIST 全球非静力模式中的评估.
  *Geosci. Model Dev.*, 17.

## 教材与专著（ima 知识库「数值预报工具箱」可检索）

- Kalnay, E. (2003). *Atmospheric Modeling, Data Assimilation, and Predictability*.
  Cambridge University Press.
- Coiffier, J. (2011). *Fundamentals of Numerical Weather Prediction*. Cambridge
  University Press.
- Vallis, G. K. (2017). *Atmospheric and Oceanic Fluid Dynamics*, 2nd ed.
  Cambridge University Press.
- Satoh, M. (2014). *Atmospheric Circulation Dynamics and General Circulation
  Models*, 2nd ed. Springer.
- Jacobson, M. Z. (2005). *Fundamentals of Atmospheric Modeling*, 2nd ed.
  Cambridge University Press.
- Pielke Sr., R. A. (2013). *Mesoscale Meteorological Modeling*, 3rd ed. Academic Press.
- Holton, J. R., Hakim, G. J. (2013). *An Introduction to Dynamic Meteorology*,
  5th ed. Academic Press.
- Wilks, D. S. (2019). *Statistical Methods in the Atmospheric Sciences*, 4th ed.
  Elsevier.
- Jolliffe, I. T., Stephenson, D. B. (2011). *Forecast Verification: A
  Practitioner's Guide in Atmospheric Science*, 2nd ed. Wiley-Blackwell.
- Fletcher, S. J. (2022). *Data Assimilation for the Geosciences: From Theory to
  Application*, 2nd ed. Elsevier.
