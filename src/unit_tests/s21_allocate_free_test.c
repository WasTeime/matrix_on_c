#include "s21_matrix_test.h"

START_TEST(test_1) {
  matrix_t m1 = {0};
  int rc = s21_create_matrix(3, 4, &m1);
  int r = m1.rows;
  int c = m1.columns;
  s21_remove_matrix(&m1);

  ck_assert_int_eq(r, 3);
  ck_assert_int_eq(c, 4);
  ck_assert_int_eq(rc, OK);
}
END_TEST

START_TEST(test_2) {
  matrix_t m1 = {0};
  int rc = s21_create_matrix(1000, 1000, &m1);

  int r = m1.rows;
  int c = m1.columns;
  s21_remove_matrix(&m1);

  ck_assert_int_eq(r, 1000);
  ck_assert_int_eq(c, 1000);
  ck_assert_int_eq(rc, OK);
}
END_TEST

START_TEST(test_1_n) {
  matrix_t m1 = {0};
  int rc = s21_create_matrix(0, 4, &m1);
  ck_assert_int_eq(rc, ERROR_INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_2_n) {
  int rc = s21_create_matrix(3, 8, NULL);
  ck_assert_int_eq(rc, ERROR_INCORRECT_MATRIX);
}
END_TEST

Suite *allocate_free_suite(void) {
  Suite *s;
  TCase *tc_pos, *tc_neg;

  s = suite_create("allocate/free");

  tc_pos = tcase_create("positives");

  tcase_add_test(tc_pos, test_1);
  tcase_add_test(tc_pos, test_2);

  tc_neg = tcase_create("negatives");

  tcase_add_test(tc_neg, test_1_n);
  tcase_add_test(tc_neg, test_2_n);

  suite_add_tcase(s, tc_pos);
  suite_add_tcase(s, tc_neg);

  return s;
}

int allocate_free_test(void) {
  int no_failed = 0;
  Suite *s;
  SRunner *runner;

  printf("\n");
  s = allocate_free_suite();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_VERBOSE);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
