#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  memmoryAllocation();
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  memmoryAllocation();
  copyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept {
  rows_ = std::exchange(other.rows_, 0);
  cols_ = std::exchange(other.cols_, 0);
  matrix_ = std::exchange(other.matrix_, nullptr);
}

S21Matrix::~S21Matrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; ++i) {
      if (matrix_[i] != nullptr) {
        delete[] matrix_[i];
      }
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
  rows_ = 0;
  cols_ = 0;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  checkTwoMatrixSumSub(other);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix copy = *this;
  copy += other;
  return copy;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  checkTwoMatrixSumSub(other);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix copy = *this;
  copy -= other;
  return copy;
}

void S21Matrix::checkTwoMatrixMult(const S21Matrix &other) const {
  if (cols_ != other.rows_) {
    throw std::logic_error("Size validation error for mult");
  }
}

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (std::fabs(matrix_[i][j] - other.matrix_[i][j]) > EPS) {
        return false;
      }
    }
  }
  return true;
}

double &S21Matrix::operator()(int row, int col) {
  if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
    throw std::out_of_range("Out of range");
  }
  return matrix_[row][col];
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  checkTwoMatrixMult(other);
  S21Matrix result(rows_, other.cols_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      for (int k = 0; k < cols_; ++k) {
        result(i, j) += (matrix_[i][k] * other.matrix_[k][j]);
      }
    }
  }
  *this = std::move(result);
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix copy = *this;
  copy.MulMatrix(other);
  return copy;
}

S21Matrix &S21Matrix::operator*=(const double number) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= number;
    }
  }
  return *this;
}

void S21Matrix::MulNumber(const double value) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= value;
    }
  }
}

S21Matrix S21Matrix::operator*(const double number) {
  S21Matrix copy = *this;
  copy.MulNumber(number);
  return copy;
}

S21Matrix operator*(const double value, S21Matrix &other) {
  S21Matrix copy = other;
  copy.MulNumber(value);
  return copy;
}

bool S21Matrix::operator==(const S21Matrix &other) const noexcept {
  return EqMatrix(other);
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this == &other) {
    return *this;
  }

  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; ++i) {
      if (matrix_[i] != nullptr) {
        delete[] matrix_[i];
      }
    }
    delete[] matrix_;
  }

  rows_ = other.rows_;
  cols_ = other.cols_;
  memmoryAllocation();
  copyMatrix(other);
  return *this;
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < result.rows_; ++i) {
    for (int j = 0; j < result.cols_; ++j) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

void S21Matrix::checkSquareSize() const {
  if (rows_ != cols_) {
    throw std::logic_error("Not square matrix, cant calculate");
  }
}
void S21Matrix::calcMinor(int row, int column, const S21Matrix &other) {
  for (int i = 0, minor_i = 0; i < other.rows_; i++) {
    if (i == row) continue;
    for (int j = 0, minor_j = 0; j < other.cols_; j++) {
      if (j == column) continue;
      matrix_[minor_i][minor_j] = other.matrix_[i][j];
      minor_j++;
    }
    minor_i++;
  }
}

double S21Matrix::calcDeterminant() const {
  double det = 0;

  if (rows_ == 1) {
    det = matrix_[0][0];
  } else {
    S21Matrix minor(rows_ - 1, cols_ - 1);

    for (int i = 0, sign = 1; i < rows_; i++, sign *= -1) {
      minor.calcMinor(i, 0, *this);
      det += sign * matrix_[i][0] * minor.calcDeterminant();
    }
  }

  return det;
}

double S21Matrix::Determinant() const {
  checkSquareSize();
  double result = 0;
  if (cols_ == 1) {
    result = matrix_[0][0];
  } else {
    result = calcDeterminant();
  }

  return result;
}

S21Matrix S21Matrix::CalcComplements() const {
  checkSquareSize();
  S21Matrix result(rows_, cols_);
  if (rows_ == 1) {
    result.matrix_[0][0] = matrix_[0][0];
  } else {
    S21Matrix minor(rows_ - 1, cols_ - 1);
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        minor.calcMinor(i, j, *this);
        result.matrix_[i][j] = minor.Determinant() * pow(-1, i + j);
      }
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  checkSquareSize();
  S21Matrix result(rows_, cols_);
  double det = Determinant();
  if (fabs(det) < EPS) {
    throw std::range_error("The determinant is 0");
  }
  if (rows_ == 1) {
    result.matrix_[0][0] = 1 / det;
  } else {
    result = CalcComplements();
    result = result.Transpose();
    result.MulNumber(1.0 / det);
  }
  return result;
}

int S21Matrix::getRows() const { return rows_; }

int S21Matrix::getCols() const { return cols_; }

void S21Matrix::setRows(int rows) {
  if (rows > 0 && rows_ != rows) {
    S21Matrix result(rows, cols_);
    result.copyMatrix(*this);
    *this = result;
  }
}

void S21Matrix::setCols(int cols) {
  if (cols > 0 && cols_ != cols) {
    S21Matrix result(rows_, cols);
    result.copyMatrix(*this);
    *this = result;
  }
}

void S21Matrix::memmoryAllocation() {
  checkMatrixSize();
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::checkMatrixSize() const {
  if (rows_ < 1 || cols_ < 1) {
    throw std::invalid_argument("Invalid matrix size");
  }
}

void S21Matrix::copyMatrix(const S21Matrix &other) {
  for (int i = 0; i < rows_ && i < other.rows_; ++i) {
    for (int j = 0; j < cols_ && j < other.cols_; ++j) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  checkTwoMatrixSumSub(other);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  checkTwoMatrixSumSub(other);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::checkTwoMatrixSumSub(const S21Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Sum/sub size validation error");
  }
}

void S21Matrix::fillMatrix(double *values) {
  for (int i = 0, iv = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++, iv++) {
      matrix_[i][j] = values[iv];
    }
  }
}
