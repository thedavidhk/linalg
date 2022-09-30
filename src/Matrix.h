#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>

#include "MatrixBase.h"
#include "Timer.h"

namespace LinAlg {

template <typename _Scalar, int _Rows, int _Cols>
struct traits<Matrix<_Scalar, _Rows, _Cols>> {
    typedef _Scalar Scalar;

    static constexpr int Rows_at_compile_time = _Rows;
    static constexpr int Cols_at_compile_time = _Cols;

    static constexpr bool is_matrix = true;
};

/// Fixed size matrix
template <typename _Scalar, int _Rows, int _Cols>
class Matrix : public Matrix_base<Matrix<_Scalar, _Rows, _Cols>> {
  public:
    Matrix() {
        static_assert(_Rows * _Cols * sizeof(_Scalar) <= 2048,
                      "Limit for fixed size matrix is 2 KiB. "
                      "Use dynamic size for larger matrices.");
    }

    Matrix(std::array<_Scalar, _Rows * _Cols> data) : m_data(std::move(data)) {
        static_assert(_Rows * _Cols * sizeof(_Scalar) <= 2048,
                      "Limit for fixed size matrix is 2 KiB. "
                      "Use dynamic size for larger matrices.");
        if (data.size() != _Rows * _Cols) {
            throw std::invalid_argument(
                "data must have the same size as matrix.");
        }
    }

    static constexpr size_t rows() { return _Rows; }

    static constexpr size_t cols() { return _Cols; }

  private:
    friend Matrix_base<Matrix>;

    std::array<_Scalar, _Rows * _Cols> m_data;
    static constexpr size_t m_rows = _Rows;
    static constexpr size_t m_cols = _Cols;
};

/// Dynamic size matrix
template <typename _Scalar>
class Matrix<_Scalar, Dynamic, Dynamic>
    : public Matrix_base<Matrix<_Scalar, Dynamic, Dynamic>> {
  public:
    Matrix() = default;

    Matrix(size_t rows, size_t cols)
        : m_rows(rows), m_cols(cols), m_data(rows * cols) {}

    Matrix(size_t rows, size_t cols, std::vector<_Scalar> data)
        : m_rows(rows), m_cols(cols), m_data(std::move(data)) {
        assert(data.size() == rows * cols &&
               "data must have the same size as matrix.");
    }

    void resize(size_t rows, size_t cols) {
        m_rows = rows;
        m_cols = cols;
        m_data.resize(rows * cols);
    }

    size_t rows() { return m_rows; }

    size_t cols() { return m_cols; }

  private:
    friend Matrix_base<Matrix>;

    std::vector<_Scalar> m_data;
    size_t m_rows = 0;
    size_t m_cols = 0;
};

} // namespace Lin_alg
#endif
