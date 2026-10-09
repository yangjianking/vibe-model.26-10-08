#!/bin/bash
# =============================================================================
#  SLURM 作业模板：立方球 NWP 模式（生产预报）
# =============================================================================
#  使用：sbatch scripts/slurm_forecast.sh
#  根据集群调整：分区、节点数、CPU 数、模块加载。
# =============================================================================

#SBATCH --job-name=cubed_sph_nwp
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --partition=normal
#SBATCH --output=run_%j.out
#SBATCH --error=run_%j.err

# ---- 加载工具链模块（按集群实际模块名调整） ----
# module load gcc/11.2.0 openmpi/4.1.4 netcdf-c/4.9.0 cmake/3.24

set -euo pipefail

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# 执行（MPI 进程 × OpenMP 线程混合）
mpirun -np $((SLURM_NTASKS)) \
    ./build/src/driver/cubed_sphere_nwp \
    --config configs/default.yaml \
    --set run.nsteps=360

echo "预报完成：$?"
