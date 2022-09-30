#ifndef LINALG_MATRIXBASE_H
#define LINALG_MATRIXBASE_H

#include "Base.h"
#include "Product.h"
#include "Timer.h"

#include <random>

namespace LinAlg {

template<typename T>
concept Matrix_type = traits<T>::is_matrix;

template <typename Derived> class Matrix_base {

    using Scalar = typename traits<Derived>::Scalar;

  public:
    // Get element by single index
    Scalar& operator[](size_t index);
    const Scalar& operator[](size_t index) const;

    // Get element by row and column index
    Scalar& operator()(size_t row, size_t col);
    const Scalar& operator()(size_t row, size_t col) const;

    // Element wise addition, subtraction, and multiplication
    // These operations require both matrices to be of the exact same type,
    // i.e. to have the same dimensions and the same scalar type.
    Derived& operator+=(const Derived& rhs);
    Derived& operator-=(const Derived& rhs);
    Derived& operator*=(const Derived& rhs);

    // (Element wise) scalar multiplication
    Derived& operator*=(Scalar rhs);

    // Fill matrix with random numbers
    void randomize_uniform(Scalar lower, Scalar upper);

    // Pretty print matrix
    void print(int precision = 3) const;

    // Get dimensions of matrix
    constexpr size_t rows() const { return derived().rows(); }
    constexpr size_t cols() const { return derived().cols(); }
    constexpr size_t size() const { return rows() * cols(); }

    // Return iterators to the underlying data container
    auto begin();
    const auto begin() const;
    auto end();
    const auto end() const;

  private:
    Matrix_base(){};
    friend Derived;

    // Returns a reference to the derived instance
    Derived& derived();
    const Derived& derived() const;

    // Converts row and column indices to single data index in a column major
    // fashion
    inline size_t index(size_t row, size_t col) const {
        return row + derived().m_rows * col;
    }
};

// Helper Functions

// Swap two columns of matrix
template <typename Mat>
Mat& swap_cols(Mat& mat, size_t col_index1, size_t col_index2) {
    if (col_index1 > mat.cols() - 1 || col_index2 > mat.cols() - 1) {
        throw std::range_error("Indices out of range.");
    }
    std::swap_ranges(mat.begin() + col_index1 * mat.rows(),
                     mat.begin() + (col_index1 + 1) * mat.rows(),
                     mat.begin() + col_index2 * mat.rows());
    return mat;
}

// Reorder container based on container of indices
template <typename Container, typename Index_container>
Container& reorder(Container& container, Index_container indices) {
    assert(indices.size() >= container.size() &&
           "Index must have at least the same size as container.");
    for (size_t i = 0; i < container.rows() * container.cols(); i++) {
        auto current = i;
        while (i != indices[current]) {
            auto next = indices[current];
            std::swap(container[current], container[next]);
            indices[current] = current;
            current = next;
        }
        indices[current] = current;
    }
    return container;
}

// Implementations

template <typename Derived> auto Matrix_base<Derived>::begin() {
    return derived().m_data.begin();
}

template <typename Derived> const auto Matrix_base<Derived>::begin() const {
    return derived().m_data.begin();
}

template <typename Derived> auto Matrix_base<Derived>::end() {
    return derived().m_data.end();
}

template <typename Derived> const auto Matrix_base<Derived>::end() const {
    return derived().m_data.end();
}

// Get the derived Matrix object
template <typename Derived> Derived& Matrix_base<Derived>::derived() {
    return *static_cast<Derived*>(this);
}

// Get the constant derived Matrix object
template <typename Derived>
const Derived& Matrix_base<Derived>::derived() const {
    return *static_cast<const Derived*>(this);
}

template <typename Derived>
typename Matrix_base<Derived>::Scalar&
Matrix_base<Derived>::operator()(size_t row, size_t col) {
    return derived().m_data[index(row, col)];
}

template <typename Derived>
const typename Matrix_base<Derived>::Scalar&
Matrix_base<Derived>::operator()(size_t row, size_t col) const {
    return derived().m_data[index(row, col)];
}

template <typename Derived>
typename Matrix_base<Derived>::Scalar&
Matrix_base<Derived>::operator[](size_t index) {
    return derived().m_data[index];
}

template <typename Derived>
const typename Matrix_base<Derived>::Scalar&
Matrix_base<Derived>::operator[](size_t index) const {
    return derived().m_data[index];
}

template <typename Derived>
Derived& Matrix_base<Derived>::operator+=(const Derived& rhs) {
    for (size_t index = 0; index < derived().m_rows * derived().m_cols;
         ++index) {
        derived()[index] += rhs(index);
    }
    return derived();
}

template <typename T>
T operator+(const Matrix_base<T>& lhs, const Matrix_base<T>& rhs) {
    lhs += rhs;
}

template <typename Derived>
Derived& Matrix_base<Derived>::operator-=(const Derived& rhs) {
    for (size_t index = 0; index < derived().m_rows * derived().m_cols;
         ++index) {
        derived()[index] -= rhs(index);
    }
    return derived();
}

template <typename T>
T operator-(const Matrix_base<T>& lhs, const Matrix_base<T>& rhs) {
    lhs -= rhs;
}

template <typename Derived>
Derived& Matrix_base<Derived>::operator*=(Scalar scalar) {
    for (size_t index = 0; index < derived().m_rows * derived().m_cols;
         ++index) {
        derived()[index] *= scalar;
    }
    return derived();
}

template <typename Derived>
Derived& Matrix_base<Derived>::operator*=(const Derived& rhs) {
    for (size_t index = 0; index < derived().m_rows * derived().m_cols;
         ++index) {
        derived()[index] *= rhs(index);
    }
    return *this;
}

template <typename Derived>
void Matrix_base<Derived>::randomize_uniform(Scalar lower, Scalar upper) {

    static_assert(std::is_arithmetic<Scalar>::value,
                  "Type needs to be arithmetic for random number generation.");

    typedef typename std::conditional<
        std::is_integral<Scalar>::value, std::uniform_int_distribution<Scalar>,
        std::uniform_real_distribution<Scalar>>::type Distribution;

    Distribution dist(lower, upper);

    std::default_random_engine gen;

    for (Scalar elem : derived().m_data) {
        elem = dist(gen);
    }
}

template <typename Derived>
void Matrix_base<Derived>::print(int precision) const {

    int width = 4 + precision;

    std::cout << std::fixed << std::setprecision(precision);

    for (size_t row = 0; row < derived().rows(); ++row) {
        for (size_t col = 0; col < derived().cols(); ++col) {
            std::cout << std::setw(width) << derived()(row, col);
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

} // namespace Lin_alg

#endif
