#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int flag = OK;

  if (result == NULL || rows <= 0 || columns <= 0) flag = MATRIX_ERROR;

  if (!flag) {
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix == NULL) print_err();

    for (int i = 0; i < rows; i++) {
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
      if (result->matrix[i] == NULL) print_err();
    }

    result->rows = rows;
    result->columns = columns;
  }

  return flag;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i] != NULL) free(A->matrix[i]);
    }
    free(A->matrix);
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int flag = SUCCESS;

  if (is_matrix_invalid(A) || is_matrix_invalid(B) || is_size_diff(A, B))
    flag = FAILURE;
  if (flag) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= 1e-6) {
          flag = FAILURE;
          i = A->rows;
          break;
        }
      }
    }
  }

  return flag;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = OK;

  if (is_matrix_invalid(A) || is_matrix_invalid(B) || is_size_diff(A, B))
    flag = MATRIX_ERROR;

  if (!flag && !s21_create_matrix(A->rows, A->columns, result)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return flag;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = OK;

  if (is_matrix_invalid(A) || is_matrix_invalid(B) || is_size_diff(A, B))
    flag = MATRIX_ERROR;

  if (!flag && !s21_create_matrix(A->rows, A->columns, result)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return flag;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int flag = OK;

  if (is_matrix_invalid(A)) flag = MATRIX_ERROR;

  if (!flag && !s21_create_matrix(A->rows, A->columns, result)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return flag;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = OK;

  if (is_matrix_invalid(A) || is_matrix_invalid(B) || A->columns != B->rows ||
      A->rows != B->columns)
    flag = MATRIX_ERROR;
  if (!flag && !s21_create_matrix(A->rows, B->columns, result)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        result->matrix[i][j] = 0;
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return flag;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int flag = OK;

  if (is_matrix_invalid(A)) flag = MATRIX_ERROR;

  if (!flag && !s21_create_matrix(A->columns, A->rows, result)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return flag;
}

int s21_determinant(matrix_t *A, double *result) {
  int flag = OK;
  if (is_matrix_invalid(A))
    flag = MATRIX_ERROR;
  else if (A->rows != A->columns)
    flag = CALC_ERROR;
  else if (!flag)
    *result = minor_determinant(A);
  return flag;
}

int minor_matrix(matrix_t *A, int row, int column, matrix_t *result) {
  int flag = OK;

  if (s21_create_matrix(A->rows - 1, A->columns - 1, result))
    flag = MATRIX_ERROR;
  else {
    int m_i = 0, m_j = 0;

    for (int i = 0; i < A->rows; i++) {
      if (i == row) continue;

      for (int j = 0; j < A->columns; j++) {
        if (j == column) continue;

        result->matrix[m_i][m_j] = A->matrix[i][j];
        m_j++;
      }

      m_i++;
      m_j = 0;
    }
  }
  return flag;
}

double minor_determinant(matrix_t *A) {
  int flag = OK;
  double result = 0;

  if (A->rows == 1) {
    result = A->matrix[0][0];
    flag++;
  }

  if (!flag && A->rows == 2) {
    result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    flag++;
  }
  if (!flag) {
    for (int j = 0; j < A->columns; j++) {
      matrix_t minor;
      minor_matrix(A, 0, j, &minor);
      result += A->matrix[0][j] * pow(-1, j) * minor_determinant(&minor);
      s21_remove_matrix(&minor);
    }
  }
  return result;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  //(-1)^(i + j) × M[i][j].
  int flag = OK;
  if (is_matrix_invalid(A)) flag = MATRIX_ERROR;

  if (!flag && (A->rows != A->columns)) flag = CALC_ERROR;

  if (!flag && (A->columns == 1)) flag = CALC_ERROR;

  if (!flag) {
    double determinant = 0;
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->rows; j++) {
        matrix_t minor;
        minor_matrix(A, i, j, &minor);
        s21_determinant(&minor, &determinant);
        result->matrix[i][j] = pow(-1, i + j) * determinant;
        s21_remove_matrix(&minor);
      }
    }
  }

  return flag;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  // произведение равно единице
  int flag = OK;
  double determinant;

  if (is_matrix_invalid(A)) flag = MATRIX_ERROR;
  if (!flag && A->rows != A->columns) flag = CALC_ERROR;

  if (!flag && (A->rows == 1)) {
    if (!A->matrix[0][0])
      flag = CALC_ERROR;
    else {
      s21_create_matrix(1, 1, result);
      result->matrix[0][0] = 1 / A->matrix[0][0];
    }
  } else if (!flag) {
    s21_determinant(A, &determinant);

    if (!determinant)
      flag = CALC_ERROR;

    else {
      matrix_t compl ;       // алгебр дополнения
      matrix_t trans_compl;  // транспонированная

      s21_create_matrix(A->rows, A->columns, &compl );
      s21_calc_complements(A, &compl );
      s21_transpose(&compl, &trans_compl);
      s21_mult_number(&trans_compl, 1.0 / determinant, result);

      s21_remove_matrix(&compl );
      s21_remove_matrix(&trans_compl);
    }
  }
  return flag;
}

void print_err() {
  perror("Memory error");
  exit(MATRIX_ERROR);
}

int is_matrix_invalid(matrix_t *matrix) {
  int flag = OK;
  if (matrix == NULL)
    flag = MATRIX_ERROR;
  else if (matrix->matrix == NULL)
    flag = MATRIX_ERROR;
  else if (matrix->columns <= 0 || matrix->rows <= 0)
    flag = MATRIX_ERROR;
  return flag;
}

int is_size_diff(matrix_t *A, matrix_t *B) {
  int flag = OK;
  if (A->columns != B->columns || A->rows != B->rows) flag = MATRIX_ERROR;
  return flag;
}
