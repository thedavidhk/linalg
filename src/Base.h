#ifndef LINALG_BASE_H
#define LINALG_BASE_H

namespace LinAlg {

template<typename T> struct traits;
template<typename T> struct traits<const T> : traits<T> {};

constexpr int Dynamic = -1;

template<typename T>
constexpr int Fixed_size_limit = 100000 / sizeof(T);

template<typename _Scalar, int _Rows, int _Cols>
class Matrix;
template<typename Derived> class Matrix_base;

template <typename _Scalar>
using MatrixX = Matrix<_Scalar, Dynamic, Dynamic>;

template <typename _Scalar, int size>
using Vector = Matrix<_Scalar, size, 1>;
template <typename _Scalar>
using VectorX = Matrix<_Scalar, Dynamic, 1>;

template <typename _Scalar, int size>
using Row_vector = Matrix<_Scalar, 1, size>;
template <typename _Scalar>
using Row_vectorX = Matrix<_Scalar, 1, Dynamic>;

// Double types

template <int _Rows, int _Cols>
using Matrixd = Matrix<double, _Rows, _Cols>;

template <int _Rows> using Vectord = Matrix<double, _Rows, 1>;

template <int _Cols> using Row_vectord = Matrix<double, 1, _Cols>;

using Matrix_xd = Matrix<double, Dynamic, Dynamic>;

using Vector_xd = Matrix<double, Dynamic, 1>;

using Row_vector_xd = Matrix<double, 1, Dynamic>;

// Float types

template <int _Rows, int _Cols>
using MatrixF = Matrix<float, _Rows, _Cols>;

template <int _Rows> using VectorF = Matrix<float, _Rows, 1>;

template <int _Cols> using Row_vectorF = Matrix<float, 1, _Cols>;

using Matrix_xf = Matrix<float, Dynamic, Dynamic>;

using Vector_xf = Matrix<float, Dynamic, 1>;

using Row_vector_xf = Matrix<float, 1, Dynamic>;


// Integer types

template <int _Rows, int _Cols> using MatrixI = Matrix<int, _Rows, _Cols>;

template <int _Rows> using VectorI = Matrix<float, _Rows, 1>;

template <int _Cols> using Row_vectorI = Matrix<float, 1, _Cols>;

using Matrix_xi = Matrix<int, Dynamic, Dynamic>;

using Vector_xi = Matrix<int, Dynamic, 1>;

using Row_vector_xi = Matrix<int, 1, Dynamic>;


}

#endif
