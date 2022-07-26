export CHPL_COMM_SUBSTRATE=mpi CHPL_COMM=gasnet
# export CHPL_COMM=none
export CHPL_LLVM=bundled CHPL_UNWIND=bundled
export CHPL_HOME="/vol/tcm28/westerhout_tom/chapel"
export MPI_CC="$CHPL_HOME/third-party/llvm/install/linux64-x86_64-gnu/bin/clang"
export MPI_CFLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include"
export MPI_LIBS="-L/usr/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi"
