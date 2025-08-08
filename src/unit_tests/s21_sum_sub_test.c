#include "s21_matrix_test.h"

START_TEST(test_1) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(2, 3, &m2);
  matrix_t mr = {0};
  s21_create_matrix(2, 3, &mr);

  m1.matrix[0][0] = 1.1;
  m1.matrix[0][1] = 2.2;
  m1.matrix[1][0] = 3.3;
  m1.matrix[1][1] = 4.444444444449;
  m1.matrix[0][2] = -2.0;
  m1.matrix[1][2] = 2.0;

  m2.matrix[0][0] = 0.0;
  m2.matrix[0][1] = 2.2;
  m2.matrix[1][0] = 3.3;
  m2.matrix[1][1] = 1.0;
  m2.matrix[0][2] = 2.0;
  m2.matrix[1][2] = -2.0;

  mr.matrix[0][0] = 1.1;
  mr.matrix[0][1] = 4.4;
  mr.matrix[1][0] = 6.6;
  mr.matrix[1][1] = 5.444444444449;
  mr.matrix[0][2] = 0.0;
  mr.matrix[1][2] = 0.0;

  matrix_t res = {0};

  ck_assert_int_eq(s21_sum_matrix(&m1, &m2, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&mr, &res), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&mr);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(test_2) {
  matrix_t m1 = {0};
  s21_create_matrix(1, 1, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(1, 1, &m2);
  matrix_t mr = {0};
  s21_create_matrix(1, 1, &mr);

  m1.matrix[0][0] = -1.1;

  m2.matrix[0][0] = 0.0;

  mr.matrix[0][0] = -1.1;

  matrix_t res = {0};

  ck_assert_int_eq(s21_sum_matrix(&m1, &m2, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&mr, &res), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&mr);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(test_3) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 2, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(2, 2, &m2);

  ck_assert_int_eq(s21_sum_matrix(&m1, &m2, NULL), ERROR_INCORRECT_MATRIX);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
}

START_TEST(test_1_n) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(3, 2, &m2);

  m1.matrix[0][0] = 1.1;
  m1.matrix[0][1] = 2.2;
  m1.matrix[1][0] = 3.3;
  m1.matrix[1][1] = 4.444444444449;
  m1.matrix[0][2] = -2.0;
  m1.matrix[1][2] = 2.0;

  m2.matrix[0][0] = 0.0;
  m2.matrix[0][1] = 2.2;
  m2.matrix[1][0] = 3.3;
  m2.matrix[1][1] = 1.0;
  m2.matrix[2][0] = 2.0;
  m2.matrix[2][1] = -2.0;

  matrix_t res = {0};

  ck_assert_int_eq(s21_sum_matrix(&m1, &m2, &res), ERROR_CALCULATE_MATRIX);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
}
END_TEST

START_TEST(test_2_n) {
  matrix_t m1 = {0};

  s21_create_matrix(2, 3, &m1);
  matrix_t m2 = {0};
  m2.matrix = NULL;

  matrix_t res = {0};

  ck_assert_int_eq(s21_sum_matrix(&m1, &m2, &res), ERROR_INCORRECT_MATRIX);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
}
END_TEST

START_TEST(test_3_n) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);

  matrix_t res = {0};

  ck_assert_int_eq(s21_sum_matrix(&m1, NULL, &res), ERROR_INCORRECT_MATRIX);

  s21_remove_matrix(&m1);
}
END_TEST

START_TEST(test_1_) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(2, 3, &m2);
  matrix_t mr = {0};
  s21_create_matrix(2, 3, &mr);

  m1.matrix[0][0] = 1.1;
  m1.matrix[0][1] = 2.2;
  m1.matrix[1][0] = 3.3;
  m1.matrix[1][1] = 4.444444444449;
  m1.matrix[0][2] = -2.0;
  m1.matrix[1][2] = 2.0;

  m2.matrix[0][0] = 0.0;
  m2.matrix[0][1] = 2.2;
  m2.matrix[1][0] = 3.3;
  m2.matrix[1][1] = 1.0;
  m2.matrix[0][2] = 2.0;
  m2.matrix[1][2] = -2.0;

  mr.matrix[0][0] = 1.1;
  mr.matrix[0][1] = 0.0;
  mr.matrix[1][0] = 0.0;
  mr.matrix[1][1] = 3.444444444449;
  mr.matrix[0][2] = -4.0;
  mr.matrix[1][2] = 4.0;

  matrix_t res = {0};

  ck_assert_int_eq(s21_sub_matrix(&m1, &m2, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&mr, &res), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&mr);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(test_2_) {
  matrix_t m1 = {0};
  s21_create_matrix(1, 1, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(1, 1, &m2);
  matrix_t mr = {0};
  s21_create_matrix(1, 1, &mr);

  m1.matrix[0][0] = -1.1;

  m2.matrix[0][0] = 0.0;

  mr.matrix[0][0] = -1.1;

  matrix_t res = {0};

  ck_assert_int_eq(s21_sub_matrix(&m1, &m2, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&mr, &res), SUCCESS);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&mr);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(test_1_n_) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(3, 2, &m2);

  matrix_t res = {0};

  ck_assert_int_eq(s21_sub_matrix(&m1, &m2, &res), ERROR_CALCULATE_MATRIX);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
}
END_TEST

START_TEST(test_2_n_) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);
  matrix_t m2 = {0};
  m2.matrix = NULL;

  matrix_t res = {0};

  ck_assert_int_eq(s21_sub_matrix(&m1, &m2, &res), ERROR_INCORRECT_MATRIX);

  s21_remove_matrix(&m1);
}
END_TEST

START_TEST(test_3_n_) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 3, &m1);

  matrix_t res = {0};

  ck_assert_int_eq(s21_sub_matrix(&m1, NULL, &res), ERROR_INCORRECT_MATRIX);

  s21_remove_matrix(&m1);
}
END_TEST

START_TEST(test_4_n_) {
  matrix_t m1 = {0};
  s21_create_matrix(2, 2, &m1);
  matrix_t m2 = {0};
  s21_create_matrix(2, 2, &m2);

  ck_assert_int_eq(s21_sub_matrix(&m1, &m2, NULL), ERROR_INCORRECT_MATRIX);

  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
}

Suite *sum_suite(void) {
  Suite *s;
  TCase *tc_pos, *tc_neg;

  s = suite_create("sum");

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

Suite *sub_suite(void) {
  Suite *s;
  TCase *tc_pos, *tc_neg;

  s = suite_create("sub");

  tc_pos = tcase_create("positives");

  tcase_add_test(tc_pos, test_1_);
  tcase_add_test(tc_pos, test_2_);
  tcase_add_test(tc_pos, test_3);

  tc_neg = tcase_create("negatives");

  tcase_add_test(tc_neg, test_1_n_);
  tcase_add_test(tc_neg, test_2_n_);
  tcase_add_test(tc_neg, test_3_n_);
  tcase_add_test(tc_neg, test_4_n_);

  suite_add_tcase(s, tc_pos);
  suite_add_tcase(s, tc_neg);

  return s;
}

int sum_sub_test(void) {
  int no_failed = 0;
  Suite *s;
  SRunner *runner;

  printf("\n");
  s = sum_suite();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_VERBOSE);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  printf("\n");
  s = sub_suite();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_VERBOSE);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
