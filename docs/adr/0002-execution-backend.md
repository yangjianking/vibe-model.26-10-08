# ADR-0002：执行后端抽象 —— Kokkos 首选 + 退化路径

- **状态**：已采纳
- **日期**：2026-10-08

## 背景

提示词 P3 硬性要求：GPU 与异构加速一律通过 Kokkos 抽象层接入，业务代码禁止
裸 CUDA/HIP；同时要求覆盖 CPU 众核（OpenMP + NUMA + SIMD）与 GPU 两条路径，
共用同一份内核源码。

## 决策

优先使用 Kokkos（`Kokkos::View` / `parallel_for` / `parallel_reduce`）。当构建
环境未安装 Kokkos 时，退化为**自研 C++ 模板执行策略层**（同语义），保证同一
份业务代码无需修改即可在纯 CPU 环境编译运行。

## 后果

- 业务代码通过 `parallel_for_range` / `parallel_reduce_range` / `View` 抽象，
  平台差异仅存在于 `common/execution_policy.hpp`；
- 退化路径用于本地开发与 CI，生产 HPC 集群用 Kokkos 后端（CUDA/HIP/SYCL）。

## 参考文献

- Edwards et al. (2014), *J. Parallel Distrib. Comput.*, 74, 3202–3216.
- Adams et al. (2019), *J. Parallel Distrib. Comput.*, 132, 383–396.
