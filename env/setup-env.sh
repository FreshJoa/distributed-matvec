# export PKG_CONFIG_PATH=/home/tom/src/lattice-symmetries/prefix/lib/pkgconfig
export LD_LIBRARY_PATH="$PWD/third_party/lib:$LD_LIBRARY_PATH"
# /home/tom/src/lattice-symmetries/prefix/lib:/home/tom/src/lattice-symmetries-haskell/build
export CHPL_COMM_SUBSTRATE=smp CHPL_COMM=gasnet
# export CHPL_COMM_SUBSTRATE=
# export CHPL_COMM=none
export CHPL_LLVM=system
export CHPL_UNWIND=none
export CHPL_TARGET_CPU=native
export CHPL_HWLOC=bundled
export CHPL_RE2=none
. "${HOME}/src/chapel-1.26.0/util/setchplenv.bash"
export CHPL_RT_OVERSUBSCRIBED=yes
# export CHPL_RT_NUM_THREADS_PER_LOCALE=1
export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
