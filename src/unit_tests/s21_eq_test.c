#include "s21_matrix_test.h"

START_TEST(test_1) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 2, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(2, 2, &m2);

  m1.matrix[0][0] = 1.1;
  m1.matrix[0][1] = 2.2;
  m1.matrix[1][0] = 3.3;
  m1.matrix[1][1] = 4.444444444449;

  m2.matrix[0][0] = 1.1;
  m2.matrix[0][1] = 2.2;
  m2.matrix[1][0] = 3.3;
  m2.matrix[1][1] = 4.4444444444444444;

  ck_assert_int_eq(s21_eq_matrix(&m1, &m2), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
}
END_TEST

START_TEST(test_2) {
  matrix_t m1 = {0};
  s21_create_matrix(1, 1, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(1, 1, &m2);

  m1.matrix[0][0] = 123.456789;

  m2.matrix[0][0] = 123.456789;

  ck_assert_int_eq(s21_eq_matrix(&m1, &m2), SUCCESS);

  s21_remove_matrix(&m1);

  s21_remove_matrix(&m2);
}
END_TEST

START_TEST(test_1_n) {
  matrix_t m1 = {0};
  s21_create_matrix(1, 1, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(1, 1, &m2);

  m1.matrix[0][0] = 123.456799;

  m2.matrix[0][0] = 123.456789;

  ck_assert_int_eq(s21_eq_matrix(&m1, &m2), FAILURE);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
}
END_TEST

START_TEST(test_2_n) {
  matrix_t m1 = {0};
  s21_create_matrix(1, 1, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(2, 1, &m2);

  m1.matrix[0][0] = 123.456789;

  m2.matrix[0][0] = 123.456789;
  m2.matrix[1][0] = 123.0;

  ck_assert_int_eq(s21_eq_matrix(&m1, &m2), FAILURE);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
}
END_TEST

START_TEST(test_3_n) {
  matrix_t m2 = {0};
  s21_create_matrix(2, 1, &m2);

  m2.matrix[0][0] = 123.456789;
  m2.matrix[1][0] = 123.0;

  ck_assert_int_eq(s21_eq_matrix(NULL, &m2), FAILURE);
  s21_remove_matrix(&m2);
}
END_TEST

START_TEST(test_4_n) {
  matrix_t m1 = {0};
  m1.matrix = NULL;
  m1.rows = 2;
  m1.columns = 1;
  matrix_t m2 = {0};
  s21_create_matrix(2, 1, &m2);

  m2.matrix[0][0] = 123.456789;
  m2.matrix[1][0] = 123.0;

  ck_assert_int_eq(s21_eq_matrix(&m1, &m2), FAILURE);

  s21_remove_matrix(&m2);
}
END_TEST

Suite *eq_suite(void) {
  Suite *s;
  TCase *tc_pos, *tc_neg;

  s = suite_create("eq");

  tc_pos = tcase_create("positives");

  tcase_add_test(tc_pos, test_1);
  tcase_add_test(tc_pos, test_2);

  tc_neg = tcase_create("negatives");

  tcase_add_test(tc_neg, test_1_n);
  tcase_add_test(tc_neg, test_2_n);
  tcase_add_test(tc_neg, test_3_n);
  tcase_add_test(tc_neg, test_4_n);

  suite_add_tcase(s, tc_pos);
  suite_add_tcase(s, tc_neg);

  return s;
}

int eq_test(void) {
  int no_failed = 0;
  Suite *s;
  SRunner *runner;

  printf("\n");
  s = eq_suite();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_VERBOSE);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
