export CHPL_HOME="/vol/tcm28/westerhout_tom/chapel"
export GASNET_BACKTRACE=1
export CHPL_COMM=gasnet
export CHPL_COMM_SUBSTRATE=ofi
export CHPL_LAUNCHER=gasnetrun_ofi
export CHPL_GASNET_MORE_CFG_OPTIONS=--with-ofi-provider=tcp
export GASNET_OFI_SPAWNER=ssh

export CHPL_LLVM=bundled
export CHPL_UNWIND=bundled
export CHPL_HWLOC=bundled
export CHPL_GMP=none
export CHPL_RE2=none
