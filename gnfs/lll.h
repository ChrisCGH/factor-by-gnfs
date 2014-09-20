#ifndef LLL_H
#define LLL_H
int LLL_reduce_1_on_columns(Matrix<VeryLong>& b, Matrix<VeryLong>& H);
int LLL_reduce_1(Matrix<VeryLong>& b, Matrix<VeryLong>& H);
int LLL_reduce_2_on_columns(Matrix<VeryLong>& b);
int LLL_reduce_3_on_columns(Matrix<VeryLong>& b);
int LLL_reduce_2(Matrix<VeryLong>& b);
Matrix<VeryLong> kernel_over_Z_using_LLL(const Matrix<VeryLong>& A);
#endif
