#include "s21_matrix_test.h"

int main(void) {
  int no_failed = 0;
  no_failed += allocate_free_test();
  no_failed += eq_test();
  no_failed += sum_sub_test();
  no_failed += mul_test();
  no_failed += tr_test();
  no_failed += det_test();
  no_failed += cm_test();
  no_failed += minor_test();
  no_failed += inverse_test();

  printf("\n");

  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
