#include <stdio.h>

void print_external_string(char const *str) {
  fprintf(stderr, "print_external_string received '%s'\n", str);
}
