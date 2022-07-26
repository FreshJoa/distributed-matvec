#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>

long matvec_get_thread_id() { return syscall(SYS_gettid); }
