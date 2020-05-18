#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <string>
#include <initializer_list>

/**
 * Offers useful classes and functions for work with matrices. Main class is
 * Matrix class that represents a real valued matrix. Binary operators like 
 * matrix addition, matrix multiplication, scalar multiplication and similar 
 * can be used as operators.
 * 
 * For solving linear systems and computing inverses and determinants class LU 
 * can be used.
 * 
 * @author Antun Magdic
 */
namespace matlib {

/**
 * Real valued matrix with basic algebraic operations. double is used for
 * storing values. Indexing starts at 1.
 */
class Matrix {
public:

    /**
     * Creates a new Matrix with given number of rows and columns. Values 
     * of elements aren't initialized to any value. If either given number of 
     * rows or given number of columns is negative or 0, behaviour is 
     * undefined.
     * 
     * @param rows number of rows of a matrix to create
     * @param cols number of columns of a matrix to create
     */
    Matrix(int rows, int cols);

    /**
     * Creates a new Matrix from given data. Elements of the matrix to create 
     * are stored at data row by row. Matrix has rows rows and cols columns.
     * If use_data is set to true, given memory will be used by Matrix without 
     * copying. By default data is copied for use by Matrix.
     * If memory at data isn't at least rows*cols big, behaviour is undefined.
     * 
     * @param data pointer to matrix elements stored in memory row by row
     * @param rows number of rows of matrix to create
     * @param cols number of columns of matrix to create
     * @param use_data if set to true values at data will not be copied, but 
     *        given memory will be used
     */
    Matrix(double* data, int rows, int cols, bool use_data=false);

    /**
     * Creates a new Matrix from given data. Pointers to rows of matrix to 
     * create are stored at row_pointers. Matrix has rows rows and cols 
     * columns.
     * If row_pointers isn't at least rows big or if any row_pointers[i] isn't
     * at least cols big, behaviour is undefined.
     * 
     * @param row_pointers pointer to the array from which every element points
     *        to an array of values in a single row
     * @param rows number of rows
     * @param cols number of columns
     */
    Matrix(double** row_pointers, int rows, int cols);

    /**
     * Creates a new Matrix from given data.
     * 
     * @param data vector of rows
     * 
     * @throws std::logic_error if all row sizes aren't the same or if some 
     *         rows are empty or if there are zero rows
     */
    Matrix(const std::vector<std::vector<double>>& data);

    /**
     * Parses the given string str and produces a Matrix. Elements in rows are 
     * separated by tabs or spaces and rows are separated by newline 
     * characters. No empty lines are allowed between matrix rows.
     * In the following example matrix has 3 rows and 2 columns:
     * 200 3.5
     * -1.3333 -512
     * 3.14 2.71
     * 
     * @param str string with matrix elements
     * 
     * @throws std::invalid_argument if the given string isn't formatted 
     *         correctly
     */
    Matrix(const std::string& str);

    /**
     * Copy constructor.
     */
    Matrix(const Matrix& m);

    /**
     * Move constructor.
     */
    Matrix(Matrix&& m);

    /**
     * Defined so Matrix elements can be easily initialized using initializer 
     * list.
     * 
     * Example:
     * Matrix m {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
     * 
     * @throws std::logic_error if all row sizes aren't the same
     */
    Matrix(std::initializer_list<std::initializer_list<double>> list);

    /**
     * Destructor.
     */
    ~Matrix();

    /**
     * Makes a copy of m.
     */
    Matrix& operator=(const Matrix& m);

    /**
     * Number of rows of this Matrix.
     * 
     * @return number of rows of this Matrix
     */
    int row_dim() const;

    /**
     * Number of columns of this Matrix.
     * 
     * @return number of columns of this Matrix
     */
    int col_dim() const;
    
    /**
     * Returns a reference to matrix element in ith row and jth column. 
     * Returned reference can be used as both rvalue and lvalue reference.
     * Both row indices and column indices start at 1.
     * No checking of validity of i and j is performed. If i and j aren't valid
     * for this matrix, behaviour is undefined.
     * 
     * @param i row index of element to return
     * @param j column index of element to return
     * @return reference to value in ith row and jth column
     */
    inline double& operator()(int i, int j) const {
        if (i > rows || j > cols || i < 1 || j < 1) {
            throw std::logic_error(std::to_string(i) + ", " + std::to_string(j));
        }
        return *(data + (i-1)*cols + j-1);
    }

    /**
     * Used for getting submatrices of a matrix.
     */
    struct slice {
        int from;
        int to;
    };
    
    /**
     * Returns a submatrix of this Matrix. Elements of a new matrix are 
     * specified by rows and cols slices. Both rows.from and rows.to are 
     * included as well as cols.from and cols.to.
     * Valid values for rows slice are from 1 to row dimension of this matrix 
     * (both included). Valid values for cols slice are from 1 to column 
     * dimension of this matrix (both included). rows.from should be less than 
     * or equal to rows.to and cols.from should be less than or equal to 
     * cols.to. Validity of given arguments isn't checked and if constraints 
     * given here aren't satisfied, the behaviour is undefined.
     * Returned Matrix cannot be used as lvalue to change elements of original 
     * matrix.
     * 
     * @param rows new matrix will contain rows from rows.from to rows.to (both 
     *        included)
     * @param cols new matric will contain columns from cols.from to cols.to 
     *        (both included)
     */
    Matrix operator()(slice rows, slice cols) const;

    /**
     * Adds Matrix rhs to this one.
     * 
     * @param rhs Matrix to add to this one
     * 
     * @throws dimension_mismatch_error if dimensions of this Matrix and rhs 
     *         aren't the same
     */
    Matrix& operator+=(const Matrix& rhs);

    /**
     * Subtracts Matrix rhs from this one.
     * 
     * @param rhs Matrix to subtract from this one
     * 
     * @throws dimension_mismatch_error if dimensions of this Matrix and rhs 
     *         aren't the same
     */
    Matrix& operator-=(const Matrix& rhs);

    /**
     * Adds given value to all elements of this Matrix.
     * 
     * @param value value to add to all elements of this Matrix
     */
    Matrix& operator+=(double value);

    /**
     * Subtracts given value from all elements of this Matrix.
     * 
     * @param value value to subtract from all elements of this Matrix
     */
    Matrix& operator-=(double value);

    /**
     * Multiplies all elements of this Matrix by given value.
     * 
     * @param value value with which all elements of this Matrix are multiplied
     */
    Matrix& operator*=(double value);

    /**
     * Divides all elements of this Matrix by given value.
     * 
     * @param value value by which all elements of this Matrix are divided
     */
    Matrix& operator/=(double value);

    /**
     * Returns a new Matrix that is a transpose of this one.
     * 
     * @return transpose of this Matrix
     */
    Matrix transpose();

    /**
     * Returns an identity Matrix with n rows and n columns.
     * 
     * @param n number of rows and columns of matrix to return
     * @return identity matrix with n rows and n columns
     * 
     * @throws std::logic_error if n is less than 1
     */
    static Matrix identity(int n);

    /**
     * Returns a Matrix of zeros withm rows and n columns.
     * 
     * @param m number of rows
     * @param n number of columns
     * @return matrix of zeros with m rows and n cols
     * 
     * @throws std::logic_error if m or n is less than 1
     */
    static Matrix zeros(int m, int n);

    /**
     * Loads a Matrix from a file given by its filestream ifs. File must be 
     * formatted so that a matrix can be initialized using a string in file.
     * 
     * @param ifs opened input file stream
     * @return Matrix read from given file
     * 
     * @throws std::invalid_argument if the given string isn't formatted 
     *         correctly
     */
    static Matrix from_file(std::ifstream& ifs);

    /**
     * Prints this Matrix to given file in a way that can be used as input to
     * Matrix::from_file.
     * 
     * @param ofs opened output file stream
     */
    void to_file(std::ofstream& ofs);

    /**
     * Checks if given matrices are equal. Matrices are equal if they have
     * equal dimensions and if all values are equals. Two values a and b are 
     * considered if abs(a - b) <= eps.
     * 
     * @param lhs
     * @param rhs
     * @param eps used to compare values (default is 1e-5)
     */
    static bool equals(const Matrix& lhs, const Matrix& rhs, double eps=1e-5);

    /**
     * Returns string representation of this Matrix. Returned string can be 
     * used as valid input to Matrix::parse.
     * 
     * @param precision number of digits after decimal point (default is 3)
     * @return string representation of this Matrix
     */
    std::string to_string(int precision=3) const;

private:

    /** 
     * Pointer to the array used to store matrix elements.
     * Element in the i-th row and j-th column can be retrieved like this:
     *      elem_ij = *(data + (i-1) * cols + j-1)
     * (Indices start at 1)
     */
    double* data;

    /** number of rows */
    int rows;

    /** number of columns */
    int cols;

    /** Parses the given string and produces a vector of rows */
    static std::vector<std::vector<double>> parse(const std::string& str);

};

/**
 * Adds matrices lhs and rhs together and produces a new Matrix as the result 
 * of addition.
 * 
 * @param lhs first operand
 * @param rhs second operand
 * @return new Matrix that is the result of addition of lhs and rhs
 * 
 * @throws dimension_mismatch_error if dimensions of lhs and rhs aren't the 
 *         same
 */
Matrix operator+(const Matrix& lhs, const Matrix& rhs);

/**
 * Subtracts Matrix rhs from Matrix lhs and produces a new Matrix as the result
 * of subtraction.
 * 
 * @param lhs first operand
 * @param rhs second operand
 * @return new Matrix that is the result of subtraction of rhs from lhs
 * 
 * @throws dimension_mismatch_error if dimensions of lhs and rhs aren't the
 *         same
 */
Matrix operator-(const Matrix& lhs, const Matrix& rhs);

/**
 * Multiplies all elements of rhs by given value lhs and produces a new Matrix
 * as the result of scalar multiplication.
 * 
 * @param lhs scalar
 * @param rhs matrix
 * @return new Matrix that is the result of scalar multiplication of rhs by lhs
 */
Matrix operator*(double lhs, const Matrix& rhs);

/**
 * Divides all elements of lhs by given value rhs and produces a new Matrix
 * as the result of scalar division.
 * 
 * @param lhs matrix
 * @param rhs scalar
 * @return new Matrix that is the result of scalar division of lhs by rhs
 */
Matrix operator*(const Matrix& lhs, double rhs);

/**
 * Mulitplies matrices lhs and rhs and produces a new Matrix as the result of 
 * matrix multiplication.
 * 
 * @param lhs first operand
 * @param rhs second operand
 * @return new Matrix that is the result of matrix multiplication of lhs and 
 *         rhs
 * 
 * @throws dimension_mismatch_error if column dimension of lhs isn't the same 
 *         as row dimension of rhs
 */
Matrix operator*(const Matrix& lhs, const Matrix& rhs);

/**
 * Equivalent to Matrix::equals(lhs, rhs).
 */
bool operator==(const Matrix& lhs, const Matrix& rhs);

/**
 * Stacks matrices m1 and m2 horizontally. m1 and m2 must have the same row
 * dimension. Column dimension of the new matrix will be the sum of column
 * dimensions of m1 and m2.
 * If matrices m1 and m2 don't have the same row dimension, the behaviour is 
 * undefined.
 * 
 * @param m1 first matrix
 * @param m2 second matrix
 * @return Matrix gotten from appending columns of m2 to m1
 * 
 * @throws dimension_mismatch_error if row dimensions of m1 and m2 aren't the 
 *         same
 */
Matrix cbind(const Matrix& m1, const Matrix& m2);

/**
 * Stacks matrices m1 and m2 vertically. m1 and m2 must have the same column 
 * dimension. Row dimension of the new matrix will be the sum of row dimensions
 * of m1 and m2.
 * 
 * @param m1 first matrix
 * @param m2 second matrix
 * @return Matrix gotten from appending rows of m2 to m1
 * 
 * @throws dimension_mismatch_error if column dimensions of m1 and m2 aren't
 *         the same
 */
Matrix rbind(const Matrix& m1, const Matrix& m2);

/**
 * Solves the system AX = B in matrix form. A must be square non singular 
 * matrix and B must have the same number of rows as A.
 * 
 * @param a A
 * @param b B
 * @param eps used to check if A is singular
 * @return X
 * 
 * @throws dimension_mismatch_error if B doesn't have the same number of rows 
 *         as A
 * @throws singular_matrix_error if A is singular
 */
Matrix solve(const Matrix& a, const Matrix& b, double eps=1E-8);

/**
 * Computes the inverse of Matrix m. m must be square non singular matrix.
 * 
 * @param m
 * @param eps used to check if A is singular
 * @return inverse of m
 * 
 * @throws singular_matrix_error if m is singular
 */
Matrix inv(const Matrix& m, double eps=1E-8);

/**
 * Computes the determinant of Matrix m. m must be square matrix.
 * 
 * @param m
 * @return determinant of m
 */
double det(const Matrix& m);

/**
 * LU (or LUP) decomposition of a matrix. This class can be used to decompose 
 * a matrix to a product of lower triangular and upper triangular matrices.
 * 
 * Using this class inverse and determinant of a matrix can be calculated 
 * easily and efficiently.
 * 
 * This class can be also used to solve linear systems of n equations with n
 * unknowns and matrix equations of form AX = B, where A, X and B are matrices.
 * 
 * LU (or LUP) decomposition can only be performed on square nonsingular 
 * matrices.
 * 
 * In the case of LU decomposition matrix A is decomposed to a product LU, 
 * where L is a lower triangular matrix with ones on its main diagonal and U is
 * an upper triangular matrix.
 * 
 * In the case of LUP decomposition matrix PA is decomposed to a product LU,
 * where L is a lower triangular matrix with ones on its main diagonal and U is
 * an upper triangular matrix. Matrix P is a permutation matrix and is used to 
 * avoid small numbers on main diagonal while performing the decomposition.
 */
class LU {
public:

    /**
     * Creates an object that stores information about LU decomposition of 
     * matrix m. If exchange_rows is true LUP decomposition is performed, 
     * otherwise LU decomposition is performed.
     * 
     * @param m matrix to decompose
     * @param exchange_rows if true LUP is performed, otherwise LU is performed
     * @param eps used to check if matrix is singular
     * 
     * @throws matlib_error if matrix m isn't square
     * @throws singular_matrix_error if matrix is singular or if exchange_rows
     *         is true and LU decomposition cannot be performed because of a 
     *         zero (or any number whose absolute value is less than eps) on
     *         main diagonal
     */
    LU(const Matrix& m, bool exchange_rows=true, double eps=1E-8);

    /**
     * Desturctor.
     */
    ~LU();

    /**
     * Returns L matrix.
     * 
     * @return L
     */
    Matrix l();

    /**
     * Returns U matrix.
     * 
     * @return U
     */
    Matrix u();

    /**
     * Solves the system AX = B in matrix form. B must have the same number of
     * rows as A. A is the matrix that was used to construct this LU object.
     * 
     * @param m B
     * @param eps used to check for singular matrix
     * @return X
     * 
     * @throws dimension_mismatch_error if B doesn't have the same number of rows 
     *         as A
     * @throws singular_matrix_error if LU is singular
     */
    Matrix solve(const Matrix& m, double eps=1E-8);

    /**
     * Computes the inverse of the matrix that was used to construct this LU 
     * object.
     * 
     * @param eps used to check for singular matrix
     * @return inverse of the matrix used to construct this LU object
     * 
     * @throws singular_matrix_error if LU is singular
     */
    Matrix inverse(double eps=1E-8);

    /**
     * Computes the determinant of the matrix that was used to construct this 
     * LU object.
     */
    double determinant();
    
private:

    /** L and U matrices stored in one Matrix */
    Matrix lu;

    /** 
     * Used to store permutation of rows. Value at ith index is the position 
     * of row i from the original matrix. Value at index 0 is the number of row
     * exchanges made during the decomposition process. If this is an LU 
     * decomposition (not LUP) p is nullptr.
     */
    int *p = nullptr;

    /** Performs LU decomposition of Matrix lu. */
    void lu_decomposition(double eps);

    /** Performs LUP decomposition of Matrix lu. */
    void lup_decomposition(double eps);

    /** Performs forward substitution. */
    Matrix forward_subs(const Matrix& b);

    /** Performs backward substitution. */
    Matrix backward_subs(const Matrix& y, double eps=1E-8);

};

/**
 * Used to indicate an error that occurred while performing some of the matrix 
 * operations.
 */
class matlib_error : public std::logic_error {
public:
    explicit matlib_error(const std::string& str) :std::logic_error {str} { }
    explicit matlib_error(const char* str) :std::logic_error {str} { }
};

/**
 * Used to indicate invalid dimensions of operands or arguments.
 */
class dimension_mismatch_error : public matlib_error {
public:
    explicit dimension_mismatch_error(const std::string& str) 
    :matlib_error {str} { }
    explicit dimension_mismatch_error(const char* str) 
    :matlib_error {str} { }
    explicit dimension_mismatch_error(const Matrix& m1, const Matrix& m2)
    :matlib_error {error_string(m1, m2)} { }
private:
    /** returns error message for dimension mismatch of m1 and m2 */
    static std::string error_string(const Matrix& m1, const Matrix& m2);
};

/**
 * Used to indicate a singular matrix that has been encountered in operation
 * which is not defined for singular matrices.
 */
class singular_matrix_error : public matlib_error {
public:
    explicit singular_matrix_error(const std::string& str) 
    :matlib_error {str} { }
    explicit singular_matrix_error(const char* str)
    :matlib_error {str} { }
};

} // namespace matlib

std::ostream& operator<<(std::ostream& out, matlib::Matrix& m);

std::ostream& operator<<(std::ostream& out, matlib::Matrix&& m);

std::ostream& operator<<(std::ostream& out, matlib::Matrix::slice& s);

std::ostream& operator<<(std::ostream& out, matlib::LU& lu);

#endif // MATRIX_H
