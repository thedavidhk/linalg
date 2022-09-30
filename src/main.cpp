#include <algorithm>
#include <iostream>

#include "LinAlg.h"
#include "Timer.h"

int main() {
    using namespace LinAlg;

    LA_FUNC_TIMER;

    Matrix<double, 3, 3> A({0, 2, 0, 1, 1, 2, 2, 0, 1});
    A.print();

    Matrix<double, 3, 1> b({1, 2, 3});
    b.print();

    std::cout << "Solve Ax = b:\n";

    solve(A, b).print(4);

    // A.swap_cols(0, 2);
    // A.print(4);

    // Matrix<double, 3, 3> B = LUDecompose(A);
    // B.print();

    /*
    static constexpr size_t n = 100;

    Matrix<double> A(n, n);
    A.randomize_uniform(0, 2);

    Matrix<double> B(n, n);
    B.randomize_uniform(0, 2);

    product(A, B);

    MatrixX<double> C(n, n);
    C.randomize_uniform(0, 2);

    MatrixX<double> D(n, n);
    D.randomize_uniform(0, 2);

    product(C, D);
    */
}
