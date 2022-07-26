module load 2021
module load GCC/10.3.0
export CHPL_COMM_SUBSTRATE=ibv CHPL_COMM=gasnet
export CHPL_LLVM=bundled CHPL_UNWIND=bundled CHPL_GMP=none CHPL_HWLOC=bundled CHPL_RE2=none
export CHPL_HOME="$HOME/chapel"
export CHPL_LAUNCHER=slurm-gasnetrun_ibv
