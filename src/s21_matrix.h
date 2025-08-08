#ifndef __S21_MATRIX_H__
#define __S21_MATRIX_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define OK 0
#define ERROR_INCORRECT_MATRIX 1
#define ERROR_CALCULATE_MATRIX 2

#define SUCCESS 1
#define FAILURE 0
#define eps 1e-7

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);

int s21_eq_matrix(matrix_t *A, matrix_t *B);

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

int create_minor(matrix_t *matrix, int row, int col, matrix_t *result);

#endif