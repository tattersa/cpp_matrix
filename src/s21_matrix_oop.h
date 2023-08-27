#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <stdexcept>
#include <utility>

#define EPS 1e-7

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other) noexcept;
  ~S21Matrix();

  bool EqMatrix(const S21Matrix &other) const;
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double value);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;

  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix operator*(const double number);
  bool operator==(const S21Matrix &other) const noexcept;
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double number);
  double &operator()(int row, int col);

  int getRows() const;
  int getCols() const;

  void setRows(int rows);
  void setCols(int cols);
  void fillMatrix(double *values);

 private:
  int rows_, cols_;
  double **matrix_;

  void memmoryAllocation();
  void checkMatrixSize() const;
  void copyMatrix(const S21Matrix &other);
  void checkTwoMatrixSumSub(const S21Matrix &other) const;
  void checkTwoMatrixMult(const S21Matrix &other) const;
  void checkSquareSize() const;
  double calcDeterminant() const;
  void calcMinor(int row, int column, const S21Matrix &other);
};
S21Matrix operator*(const double value, const S21Matrix &other);

#endif  // SRC_S21_MATRIX_OOP_H_
