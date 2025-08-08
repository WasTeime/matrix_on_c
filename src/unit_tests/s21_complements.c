#include "s21_matrix_test.h"

START_TEST(test_1) {
  matrix_t m1 = {0};
  s21_create_matrix(3, 3, &m1);

  matrix_t mr = {0};
  s21_create_matrix(3, 3, &mr);

  m1.matrix[0][0] = 1.0;
  m1.matrix[0][1] = 2.0;
  m1.matrix[0][2] = 3.0;
  m1.matrix[1][0] = 4.0;
  m1.matrix[1][1] = 7.0;
  m1.matrix[1][2] = 9.0;
  m1.matrix[2][0] = 87.0;
  m1.matrix[2][1] = 65.0;
  m1.matrix[2][2] = 43.0;

  mr.matrix[0][0] = -284.0;
  mr.matrix[0][1] = 611.0;
  mr.matrix[0][2] = -349.0;
  mr.matrix[1][0] = 109.0;
  mr.matrix[1][1] = -218.0;
  mr.matrix[1][2] = 109.0;
  mr.matrix[2][0] = -3.0;
  mr.matrix[2][1] = 3.0;
  mr.matrix[2][2] = -1.0;

  matrix_t res = {0};

  ck_assert_int_eq(s21_calc_complements(&m1, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&mr, &res), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&mr);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(test_2) {
  matrix_t m1 = {0};
  s21_create_matrix(3, 3, &m1);

  matrix_t mr = {0};
  s21_create_matrix(3, 3, &mr);

  m1.matrix[0][0] = -1.1;
  m1.matrix[0][1] = 2.0;
  m1.matrix[0][2] = 3.0;
  m1.matrix[1][0] = 4.0;
  m1.matrix[1][1] = -7.7;
  m1.matrix[1][2] = 9.0;
  m1.matrix[2][0] = 87.0;
  m1.matrix[2][1] = 65.0;
  m1.matrix[2][2] = -43.43;

  mr.matrix[0][0] = -250.589;
  mr.matrix[0][1] = 956.72;
  mr.matrix[0][2] = 929.9;
  mr.matrix[1][0] = 281.86;
  mr.matrix[1][1] = -213.227;
  mr.matrix[1][2] = 245.5;
  mr.matrix[2][0] = 41.1;
  mr.matrix[2][1] = 21.9;
  mr.matrix[2][2] = 0.47;

  matrix_t res = {0};

  ck_assert_int_eq(s21_calc_complements(&m1, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&mr, &res), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&mr);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(test_1_n) {
  matrix_t m1 = {0};
  s21_create_matrix(1, 2, &m1);

  matrix_t res = {0};

  ck_assert_int_eq(s21_calc_complements(&m1, &res), ERROR_CALCULATE_MATRIX);

  s21_remove_matrix(&m1);
}
END_TEST

START_TEST(test_2_n) {
  matrix_t m2 = {0};
  m2.matrix = NULL;

  matrix_t res = {0};

  ck_assert_int_eq(s21_calc_complements(&m2, &res), ERROR_INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_3_n) {
  matrix_t res = {0};

  ck_assert_int_eq(s21_calc_complements(NULL, &res), ERROR_INCORRECT_MATRIX);
}
END_TEST

Suite *cm_suite(void) {
  Suite *s;
  TCase *tc_pos, *tc_neg;

  s = suite_create("complements");

  tc_pos = tcase_create("positives");

  tcase_add_test(tc_pos, test_1);
  tcase_add_test(tc_pos, test_2);

  tc_neg = tcase_create("negatives");

  tcase_add_test(tc_neg, test_1_n);
  tcase_add_test(tc_neg, test_2_n);
  tcase_add_test(tc_neg, test_3_n);

  suite_add_tcase(s, tc_pos);
  suite_add_tcase(s, tc_neg);

  return s;
}

int cm_test(void) {
  int no_failed = 0;
  Suite *s;
  SRunner *runner;

  printf("\n");
  s = cm_suite();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_VERBOSE);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
