/**
 * This simple program uses LU decomposition to solve a linear system Ax = b.
 */

#include "matlib.h"

#include <iostream>

using namespace matlib;

int main() {
    Matrix A {
        {0.000001, 3000000, 2000000},
        { 1000000, 2000000, 3000000},
        { 2000000, 1000000, 2000000}
    };
    Matrix b {
        {12000000.000001},
        {14000000},
        {10000000}
    };

    // Perform LU decomposition of matrix A
    LU lu {A};
    
    // Solve Ax = b
    Matrix x = lu.solve(b);

    // Alternative:
    // Matrix x = solve(A, b);

    std::cout << "Solution is:" << std::endl;
    std::cout << lu.solve(b) << std::endl;

    return 0;
}
