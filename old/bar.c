#include "foo.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
    chpl_library_init(argc, argv);

    printf("%f\n", bar(7)); // Call into a library function

    chpl_library_finalize();

    return 0;
}
