# s21_matrix

Implementation of matrix.h library using C programming language.

## Contents

1. [Information](#information) 
2. [Library functions](#library-functions) 
3. [Build](#build) 

### Library functions
| Name                   | Description                                              |
| ---------------------- | -------------------------------------------------------- |
| `s21_get_minor`        |  Get minor value at selected position in matrix          |
| `s21_eq_matrix`        |  Determine if matrices are equal.                        |
| `s21_transpose`        |  Switch around rows and columns of matrix into new one.  |
| `s21_sub_matrix`       |  Substruct matrices.                                     |
| `s21_sum_matrix`       |  Summarize matrices.                                     |
| `s21_mult_matrix`      |  Multiplicate matrices.                                  |
| `s21_mult_number`      |  Multiply matrix by number.                              |
| `s21_determinant`      |  Calculate determinant of matrix.                        |
| `s21_copy_matrix`      |  Copy matrix to new matrix.                              |
| `s21_matrix_to_str`    |  Creates table of matrix in str with given precision.    |
| `s21_create_matrix`    |  Create matrix of given size.                            |
| `s21_remove_matrix`    |  Free matrix.                                            |
| `s21_inverse_matrix`   |  Calculate inverse matrix.                               |
| `s21_calc_complements` |  Calculate matrix cofactor matrix.                       |


### Build
* Build library and launch tests

   ```
   make all
   ```
* Build and run tests

   ```
   make test
   ```
* Build library 

   ```
   make s21_matrix.a
   ```
* Clean library

   ```
   make clean
   ```
* Run test coverage check

   ```
   make gcov_report
   ```
