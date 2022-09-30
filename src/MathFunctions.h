#ifndef LINALG_MATHFUNCTIONS_H
#define LINALG_MATHFUNCTIONS_H

#include "Base.h"
#include "MatrixBase.h"
#include "Timer.h"

#include <concepts>
#include <ranges>

namespace LinAlg {

template <typename _Lhs, typename _Rhs,
          int _Rows = traits<_Lhs>::Rows_at_compile_time,
          int _Cols = traits<_Rhs>::Cols_at_compile_time,
          typename Scalar = typename traits<_Lhs>::Scalar>
Matrix<Scalar, _Rows, _Cols> product(const _Lhs& lhs, const _Rhs& rhs) {

    LA_FUNC_TIMER;

    Matrix<Scalar, _Rows, _Cols> result;
    if constexpr (_Rows == Dynamic || _Cols == Dynamic) {
        result.resize(lhs.rows(), rhs.cols());
    }

    for (size_t row = 0; row < lhs.rows(); ++row) {
        for (size_t col = 0; col < rhs.cols(); ++col) {
            result(row, col) = 0;
            for (size_t k = 0; k < lhs.cols(); ++k) {
                result(row, col) += lhs(row, k) * rhs(k, col);
            }
        }
    }
    return result;
}

// LU-Decomposition
// The LU-Decomposition of a matrix A is given by the equation A*inv(P) = L*U
// Where P is the (column) permutation matrix and L and U are lower and upper
// triangular matrices. Since L can be represented as a unit triangular matrix
// (i.e. all diagonal elements are 1), the two matrices L and U are stored in
// the matrix LU without loss of information.
template <typename T> class LUDecomposition {

    using Scalar = typename traits<T>::Scalar;

    // Consolidated LU-Matrix
    T LU;

    // Column order permutation
    std::vector<size_t> permutation;

    // Number of column permutations (needed for determinant calculation)
    size_t perm_count = 0;

    // Constructor performing the LU-decomposition on matrix mat
    explicit LUDecomposition(const T& mat) : permutation(mat.cols()) {

        LA_FUNC_TIMER;

        // Check if square
        if (mat.rows() != mat.cols()) {
            throw std::domain_error(
                "Matrix must be square to perform LU-decomposition.");
        }

        // Column order array
        std::iota(permutation.begin(), permutation.end(), 0);

        // Tolerance TODO: Consider to psas as parameter rather than hard-code
        static constexpr double tol = 1e-6;

        for (size_t col = 0; col < mat.cols(); ++col) {
            // Check if we need to pivot
            size_t j = col;
            while (std::abs(mat(col, j)) < tol && j < mat.cols()) {
                ++j;
                if (j >= mat.cols()) {
                    throw std::domain_error("Failure: Matrix is degenerate.");
                }
            }
            if (j != col) {
                // pivot columns
                std::swap(permutation[j], permutation[col]);
                // swap_cols(LU, j, col);
                ++perm_count;
            }

            for (size_t row = 0; row < mat.rows(); ++row) {
                LU(row, col) = mat(row, permutation[col]);
                for (size_t k = 0; k < row && k < col; ++k) {
                    LU(row, col) -= LU(row, k) * LU(k, col);
                }
                if (row > col) {
                    LU(row, col) /= LU(col, col);
                }
            }
        }
    }

    template <std::ranges::view Column> Column solve_col(Column b) {

        if (LU.cols() != b.size()) {
            throw std::range_error("Dimensions do not match.");
        }

        Column intermediate(b.size());

        for (size_t row = 0; row < LU.rows(); ++row) {
            intermediate[row] = b[row];
            for (size_t k = 0; k < row; ++k) {
                intermediate[row] -= LU(row, k) * intermediate[k];
            }
        }

        Column result(b.size());

        for (size_t row = LU.rows(); row > 0; --row) {
            result[row - 1] = intermediate[row - 1];
            for (size_t k = LU.rows(); k > row; --k) {
                result[row - 1] -= LU(row - 1, k - 1) * result[k - 1];
            }
            result[row - 1] /= LU(row - 1, row - 1);
        }
        if (perm_count) {
            reorder(result, permutation);
        }

        return result;
    }

    // TODO: Unsafe. Maybe replace with span or constrain "Container" type
    // otherwise.
    template<std::ranges::viewable_range OtherT>
    std::ranges::subrange solve(OtherT B) {

        // TODO: To apply solve_col, columns of B have to be subscriptable
    }

    double determinant() {
        // Switch sign if odd number of permutations
        double det = 1.0 - 2.0 * (perm_count % 2);
        for (size_t i = 0; i < LU.cols(); ++i) {
            det *= LU(i, i);
        }
        return det;
    }
};

/*
/// Solve AX = B for X
template <typename T, typename OtherT,
          typename Scalar = typename traits<T>::Scalar>
OtherT solve(const T& A, const OtherT& b) {

    if (A.cols() != b.rows() || b.cols() != 1) {
        throw std::range_error("Dimensions do not match.");
    }

    const auto& [LU, permutation, perm_count] = LUDecomposition(A);

    OtherT intermediate;

    for (size_t row = 0; row < A.rows(); ++row) {
        intermediate[row] = b[row];
        for (size_t k = 0; k < row; ++k) {
            intermediate[row] -= LU(row, k) * intermediate[k];
        }
    }

    OtherT result;

    for (size_t row = A.rows(); row > 0; --row) {
        result[row - 1] = intermediate[row - 1];
        for (size_t k = A.rows(); k > row; --k) {
            result[row - 1] -= LU(row - 1, k - 1) * result[k - 1];
        }
        result[row - 1] /= LU(row - 1, row - 1);
    }
    if (perm_count) {
        reorder(result, permutation);
    }

    return result;
}
*/

} // namespace Lin_alg
#endif
