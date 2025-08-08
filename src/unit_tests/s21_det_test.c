#include "s21_matrix_test.h"

START_TEST(test_1) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 2, &m1);

  m1.matrix[0][0] = 1.0;
  m1.matrix[0][1] = 2.0;
  m1.matrix[1][0] = 3.0;
  m1.matrix[1][1] = 4.0;

  double res = 0.0;

  ck_assert_int_eq(s21_determinant(&m1, &res), OK);
  ck_assert_double_eq_tol(res, -2.0, eps);

  s21_remove_matrix(&m1);
}
END_TEST

START_TEST(test_2) {
  matrix_t m1 = {0};
  s21_create_matrix(4, 4, &m1);

  m1.matrix[0][0] = 1.0;
  m1.matrix[0][1] = 2.0;
  m1.matrix[0][2] = 2.0;
  m1.matrix[0][3] = 1.0;
  m1.matrix[1][0] = 5.0;
  m1.matrix[1][1] = 4.0;
  m1.matrix[1][2] = 4.0;
  m1.matrix[1][3] = 5.0;
  m1.matrix[2][0] = 7.0;
  m1.matrix[2][1] = 8.0;
  m1.matrix[2][2] = 8.0;
  m1.matrix[2][3] = 7.0;
  m1.matrix[0][0] = 111.0;
  m1.matrix[0][1] = 222.0;
  m1.matrix[0][2] = 333.0;
  m1.matrix[0][3] = 444.0;

  double res = 0.0;

  ck_assert_int_eq(s21_determinant(&m1, &res), OK);
  ck_assert_double_eq_tol(res, 0.0, eps);

  s21_remove_matrix(&m1);
}
END_TEST

START_TEST(test_3) {
  matrix_t m1 = {0};
  s21_create_matrix(1, 1, &m1);

  m1.matrix[0][0] = -1.6784684687;

  double res = 0.0;

  ck_assert_int_eq(s21_determinant(&m1, &res), OK);
  ck_assert_double_eq_tol(res, -1.6784684687, eps);

  s21_remove_matrix(&m1);
}
END_TEST

START_TEST(test_1_n) {
  matrix_t m1 = {0};
  s21_create_matrix(1, 2, &m1);

  double res = 0.0;

  ck_assert_int_eq(s21_determinant(&m1, &res), ERROR_CALCULATE_MATRIX);

  s21_remove_matrix(&m1);
}
END_TEST

START_TEST(test_2_n) {
  matrix_t m2 = {0};
  m2.matrix = NULL;

  double res = 0.0;

  ck_assert_int_eq(s21_determinant(&m2, &res), ERROR_INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_3_n) {
  double res = 0.0;

  ck_assert_int_eq(s21_determinant(NULL, &res), ERROR_INCORRECT_MATRIX);
}
END_TEST

Suite *det_suite(void) {
  Suite *s;
  TCase *tc_pos, *tc_neg;

  s = suite_create("determinant");

  tc_pos = tcase_create("positives");

  tcase_add_test(tc_pos, test_1);
  tcase_add_test(tc_pos, test_2);
  tcase_add_test(tc_pos, test_3);

  tc_neg = tcase_create("negatives");

  tcase_add_test(tc_neg, test_1_n);
  tcase_add_test(tc_neg, test_2_n);
  tcase_add_test(tc_neg, test_3_n);

  suite_add_tcase(s, tc_pos);
  suite_add_tcase(s, tc_neg);

  return s;
}

int det_test(void) {
  int no_failed = 0;
  Suite *s;
  SRunner *runner;

  printf("\n");
  s = det_suite();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_VERBOSE);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
