#include "s21_matrix.h"

int is_size_more_zero(matrix_t *A) { return (A->rows > 0 && A->columns > 0); }

int is_matrix_not_null(matrix_t *matrix) {
  return (matrix && is_size_more_zero(matrix) && matrix->matrix != NULL);
}

int is_matrix_exists_and_size_eq(matrix_t *A, matrix_t *B) {
  if (!is_matrix_not_null(A) || !is_matrix_not_null(B))
    return ERROR_INCORRECT_MATRIX;
  if (A->rows != B->rows || A->columns != B->columns)
    return ERROR_CALCULATE_MATRIX;
  return OK;
}

int is_matrix_square_and_not_null(matrix_t *A) {
  if (!is_matrix_not_null(A)) return ERROR_INCORRECT_MATRIX;
  if (A->rows != A->columns) return ERROR_CALCULATE_MATRIX;
  return OK;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (rows < 1 || columns < 1) return ERROR_INCORRECT_MATRIX;

  if (!result) return ERROR_INCORRECT_MATRIX;

  result->rows = rows;
  result->columns = columns;
  result->matrix = calloc(rows, sizeof(double *));

  if (result->matrix == NULL) return ERROR_INCORRECT_MATRIX;

  for (int i = 0; i < rows; i++) {
    result->matrix[i] = calloc(columns, sizeof(double));
    if (result->matrix[i] == NULL) return ERROR_INCORRECT_MATRIX;
    for (int j = 0; j < columns; j++) result->matrix[i][j] = 0.0;
  }

  return OK;
}

void s21_remove_matrix(matrix_t *A) {
  if (!A || !A->matrix) return;

  for (int i = 0; i < A->rows; i++) {
    free(A->matrix[i]);
  }
  free(A->matrix);
  A->matrix = NULL;

  A->rows = 0;
  A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int resp = is_matrix_exists_and_size_eq(A, B);
  if (resp != OK) return FAILURE;

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) > eps) return FAILURE;
    }
  }

  return SUCCESS;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int resp = is_matrix_exists_and_size_eq(A, B);
  if (resp != OK) return resp;

  if (result) {
    int resp = s21_create_matrix(A->rows, A->columns, result);
    if (resp != OK) return resp;
  } else
    return ERROR_INCORRECT_MATRIX;

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }

  return OK;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int resp = is_matrix_exists_and_size_eq(A, B);
  if (resp != OK) return resp;

  if (result) {
    resp = s21_create_matrix(A->rows, A->columns, result);
    if (resp != OK) return resp;
  } else
    return ERROR_INCORRECT_MATRIX;

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }

  return OK;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (!is_matrix_not_null(A)) return ERROR_INCORRECT_MATRIX;

  if (!result) return ERROR_INCORRECT_MATRIX;

  int resp = s21_create_matrix(A->rows, A->columns, result);
  if (resp != OK) return resp;

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
    }
  }

  return OK;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (!is_matrix_not_null(A) || !is_matrix_not_null(B))
    return ERROR_INCORRECT_MATRIX;
  if (A->columns != B->rows) return ERROR_CALCULATE_MATRIX;

  if (result) {
    int resp = s21_create_matrix(A->rows, B->columns, result);
    if (resp != OK) return resp;
  } else
    return ERROR_INCORRECT_MATRIX;

  for (int i = 0; i < A->rows; i++) {
    for (int k = 0; k < B->columns; k++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][k] += A->matrix[i][j] * B->matrix[j][k];
      }
    }
  }

  return OK;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (!is_matrix_not_null(A)) return ERROR_INCORRECT_MATRIX;

  if (!result) return ERROR_INCORRECT_MATRIX;

  int resp = s21_create_matrix(A->columns, A->rows, result);
  if (resp != OK) return resp;

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[j][i] = A->matrix[i][j];
    }
  }

  return OK;
}

int create_minor(matrix_t *matrix, int row, int col, matrix_t *result) {
  int resp = is_matrix_square_and_not_null(matrix);
  if (resp != OK) return resp;
  if (row < 1 || row > matrix->rows || col < 1 || col > matrix->columns)
    return ERROR_INCORRECT_MATRIX;

  if (result && matrix->rows == 1) {
    resp = s21_create_matrix(1, 1, result);
    if (resp != OK) return resp;
    result->matrix[0][0] = matrix->matrix[0][0];
    return OK;
  }

  if (result) {
    int resp = s21_create_matrix(matrix->rows - 1, matrix->columns - 1, result);
    if (resp != OK) return resp;
  } else
    return ERROR_INCORRECT_MATRIX;

  row--;
  col--;
  for (int i = 0; i < matrix->rows; i++) {
    for (int j = 0; j < matrix->columns; j++) {
      if (i < row && j < col) result->matrix[i][j] = matrix->matrix[i][j];
      if (i < row && j > col) result->matrix[i][j - 1] = matrix->matrix[i][j];
      if (i > row && j < col) result->matrix[i - 1][j] = matrix->matrix[i][j];
      if (i > row && j > col)
        result->matrix[i - 1][j - 1] = matrix->matrix[i][j];
    }
  }

  return OK;
}

int s21_determinant(matrix_t *A, double *result) {
  int resp = is_matrix_square_and_not_null(A);
  if (resp != OK) return resp;

  if (A->rows == 2 && A->columns == 2) {
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    return OK;
  } else if (A->rows == 1 && A->columns == 1) {
    *result = A->matrix[0][0];
    return OK;
  }

  matrix_t minor = {0};
  double det = 0;
  for (int i = 0; i < A->columns; i++) {
    int resp = create_minor(A, 1, i + 1, &minor);
    if (resp != OK) return resp;

    double minor_det = 0;
    s21_determinant(&minor, &minor_det);
    det += (i % 2 == 0 ? 1 : -1) * A->matrix[0][i] * minor_det;
    s21_remove_matrix(&minor);
  }
  *result = det;

  return OK;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int resp = is_matrix_square_and_not_null(A);
  if (resp != OK) return resp;

  if (!result) return ERROR_INCORRECT_MATRIX;

  resp = s21_create_matrix(A->rows, A->columns, result);
  if (resp != OK) return resp;

  matrix_t minor = {0};

  double minor_det = 0;
  resp = OK;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      resp = create_minor(A, i + 1, j + 1, &minor);
      if (resp != OK) return resp;

      resp = s21_determinant(&minor, &minor_det);
      if (resp != OK) return resp;

      result->matrix[i][j] = pow(-1, i + j) * minor_det;
      s21_remove_matrix(&minor);
    }
  }

  return OK;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int resp = is_matrix_square_and_not_null(A);
  if (resp != OK) return resp;

  if (!result) return ERROR_INCORRECT_MATRIX;

  if (A->rows == 1) {
    resp = s21_create_matrix(1, 1, result);
    if (resp != OK) return resp;
    if (fabs(A->matrix[0][0]) < eps) return ERROR_CALCULATE_MATRIX;
    result->matrix[0][0] = 1 / A->matrix[0][0];
    return OK;
  }

  resp = s21_create_matrix(A->rows, A->columns, result);
  if (resp != OK) return resp;

  double det = 0;
  resp = s21_determinant(A, &det);
  if (resp != OK) return resp;
  if (fabs(det) < eps) return ERROR_CALCULATE_MATRIX;

  resp = s21_calc_complements(A, result);
  if (resp != OK) return resp;

  matrix_t temp = {0};
  resp = s21_create_matrix(A->rows, A->columns, &temp);
  if (resp != OK) {
    s21_remove_matrix(&temp);
    return resp;
  }
  resp = s21_transpose(result, &temp);
  if (resp != OK) {
    s21_remove_matrix(&temp);
    return resp;
  }
  resp = s21_mult_number(&temp, 1 / det, result);
  if (resp != OK) {
    s21_remove_matrix(&temp);
    return resp;
  }

  s21_remove_matrix(&temp);
  return OK;
}