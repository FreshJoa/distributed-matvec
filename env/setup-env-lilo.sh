. env/chpl-env-lilo.sh
export LD_LIBRARY_PATH="$PWD/third_party/lib:$LD_LIBRARY_PATH"
source "${CHPL_HOME}/util/setchplenv.bash"
# export CHPL_RT_OVERSUBSCRIBED=yes
# export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

example05() {
    case "$1" in
        -nl)
            NUM_LOCALES=$2
            shift
            shift
            ;;
        -nl*)
            NUM_LOCALES="${1#-nl}"
            shift
            ;;
        *)
            echo "Expected -nl argument specifying the number of locales"
            return 1
    esac
    mpirun -n $NUM_LOCALES -mca btl_base_warn_component_unused 0 -- \
        bin/Example05_real -nl$NUM_LOCALES "$@"
}
