#include "s21_matrix_test.h"

START_TEST(test_1) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);

  matrix_t mr = {0};
  s21_create_matrix(3, 2, &mr);

  m1.matrix[0][0] = 1.0;
  m1.matrix[0][1] = 2.0;
  m1.matrix[0][2] = 3.0;
  m1.matrix[1][0] = 4.0;
  m1.matrix[1][1] = 5.0;
  m1.matrix[1][2] = 6.0;

  mr.matrix[0][0] = 1.0;
  mr.matrix[0][1] = 4.0;
  mr.matrix[1][0] = 2.0;
  mr.matrix[1][1] = 5.0;
  mr.matrix[2][0] = 3.0;
  mr.matrix[2][1] = 6.0;

  matrix_t res = {0};

  ck_assert_int_eq(s21_transpose(&mr, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&m1, &res), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&mr);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(test_2) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);

  matrix_t mr = {0};
  s21_create_matrix(3, 2, &mr);

  m1.matrix[0][0] = 1.0;
  m1.matrix[0][1] = 2.0;
  m1.matrix[0][2] = 3.0;
  m1.matrix[1][0] = 4.0;
  m1.matrix[1][1] = 5.0;
  m1.matrix[1][2] = 6.0;

  mr.matrix[0][0] = 1.0;
  mr.matrix[0][1] = 4.0;
  mr.matrix[1][0] = 2.0;
  mr.matrix[1][1] = 5.0;
  mr.matrix[2][0] = 3.0;
  mr.matrix[2][1] = 6.0;

  matrix_t res = {0};

  ck_assert_int_eq(s21_transpose(&m1, &res), OK);
  // ck_assert_int_eq(s21_eq_matrix(&mr, &res), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&mr);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(test_2_n) {
  matrix_t m2 = {0};
  m2.matrix = NULL;

  matrix_t res = {0};

  ck_assert_int_eq(s21_transpose(&m2, &res), ERROR_INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_3_n) {
  matrix_t res = {0};

  ck_assert_int_eq(s21_transpose(NULL, &res), ERROR_INCORRECT_MATRIX);
}
END_TEST

Suite *tr_suite(void) {
  Suite *s;
  TCase *tc_pos, *tc_neg;

  s = suite_create("transpose");

  tc_pos = tcase_create("positives");

  tcase_add_test(tc_pos, test_1);
  tcase_add_test(tc_pos, test_2);

  tc_neg = tcase_create("negatives");

  tcase_add_test(tc_neg, test_2_n);
  tcase_add_test(tc_neg, test_3_n);

  suite_add_tcase(s, tc_pos);
  suite_add_tcase(s, tc_neg);

  return s;
}

int tr_test(void) {
  int no_failed = 0;
  Suite *s;
  SRunner *runner;

  printf("\n");
  s = tr_suite();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_VERBOSE);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
