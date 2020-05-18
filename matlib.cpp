#include "matlib.h"

#include <cstring>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <string.h>

namespace matlib {

// UTILITY FUNCTION HEADERS
static bool dim_eq(const Matrix& m1, const Matrix& m2);
static std::string dim_str(const Matrix& m);
static void singular_error(double value);

Matrix::Matrix(int rows, int cols) {
    this->rows = rows;
    this->cols = cols;
    data = new double[rows * cols];
}

Matrix::Matrix(double* data, int rows, int cols, bool use_data) {
    this->rows = rows;
    this->cols = cols;
    if (use_data) {
        this->data = data;
    } else {
        this->data = new double[rows * cols];
        memcpy(this->data, data, rows*cols * sizeof(double));
    }
}

Matrix::Matrix(double** row_pointers, int rows, int cols) {
    this->rows = rows;
    this->cols = cols;
    data = new double[rows * cols];
    for (int i = 0; i < rows; i++) {
        memcpy(this->data + i * cols, row_pointers[i], cols * sizeof(double));
    }
}

Matrix::Matrix(const std::vector<std::vector<double>>& data) {
    rows = data.size();
    if (rows == 0) {
        throw std::logic_error("Data contains no rows");
    }
    cols = data[0].size();
    if (cols == 0) {
        throw std::logic_error("Column dimension 0");
    }
    this->data = new double[rows * cols];
    int i = 1;
    for (std::vector<double> row : data) {
        if (row.size() != cols) {
            throw std::logic_error(
                "Different column dimensions in different rows, row: " + 
                std::to_string(i));
        }
        int j = 1;
        for (double d : row) {
            this->operator()(i, j) = d;
            j++;
        }
        i++;
    }
}

Matrix::Matrix(const std::string& str) :Matrix {parse(str)} { }

Matrix::Matrix(const Matrix& m) {
    rows = m.rows;
    cols = m.cols;
    data = new double[rows * cols];
    memcpy(this->data, m.data, rows*cols * sizeof(double));
}

Matrix::Matrix(Matrix&& m) {
    rows = m.rows;
    cols = m.cols;
    data = m.data;
    m.data = nullptr;
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) {
    rows = list.size();
    cols = -1;
    int i = 0;
    for (auto& row : list) {
        if (cols == -1) {
            cols = row.size();
            data = new double[rows * cols];
        }
        if (cols != -1 && row.size() != cols) {
            throw std::logic_error(
                "Variable lengths of rows in initializer list");
        }
        for (double value : row) {
            data[i] = value;
            i++;
        }
    }
}

Matrix::~Matrix() {
    if (data != nullptr) {
        // data can be nullptr after move has occurred
        delete[] data;
    }
}

Matrix& Matrix::operator=(const Matrix& m) {
    if (rows * cols < m.rows * m.cols) {
        delete[] data;
        data = new double[m.rows * m.cols];
    }
    rows = m.rows;
    cols = m.cols;
    memcpy(data, m.data, rows*cols * sizeof(double));
    return *this;
}

int Matrix::row_dim() const {
    return rows;
}

int Matrix::col_dim() const {
    return cols;
}

Matrix Matrix::operator()(Matrix::slice rows, Matrix::slice cols) const {
    Matrix m {rows.to - rows.from + 1, cols.to - cols.from + 1};
    for (int i = rows.from; i <= rows.to; i++) {
        for (int j = cols.from; j <= cols.to; j++) {
            m(i - rows.from + 1, j - cols.from + 1) = this->operator()(i, j);
        }
    }
    return m;
}

Matrix& Matrix::operator+=(const Matrix& rhs) {
    if (!dim_eq(*this, rhs)) {
        throw dimension_mismatch_error(*this, rhs);
    }
    for (int i = 0; i < rows*cols; i++) {
        data[i] += rhs.data[i];
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
    if (!dim_eq(*this, rhs)) {
        throw dimension_mismatch_error(*this, rhs);
    }
    for (int i = 0; i < rows*cols; i++) {
        data[i] -= rhs.data[i];
    }
    return *this;
}

Matrix& Matrix::operator+=(double value) {
    for (int i = 0; i < rows*cols; i++) {
        data[i] += value;
    }
    return *this;
}

Matrix& Matrix::operator-=(double value) {
    for (int i = 0; i < rows*cols; i++) {
        data[i] -= value;
    }
    return *this;
}

Matrix& Matrix::operator*=(double value) {
    for (int i = 0; i < rows*cols; i++) {
        data[i] *= value;
    }
    return *this;
}

Matrix& Matrix::operator/=(double value) {
    for (int i = 0; i < rows*cols; i++) {
        data[i] /= value;
    }
    return *this;
}

Matrix Matrix::transpose() {
    Matrix res {cols, rows};
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            res(j, i) = this->operator()(i, j);
        }
    }
    return res;
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
    Matrix res {lhs};
    res += rhs;
    return res;
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs) {
    Matrix res {lhs};
    res -= lhs;
    return res;
}

Matrix operator*(double lhs, const Matrix& rhs) {
    Matrix res {rhs};
    res *= lhs;
    return res;
}

Matrix operator/(const Matrix& lhs, double rhs) {
    Matrix res {lhs};
    res /= rhs;
    return res;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.col_dim() != rhs.row_dim()) {
        throw dimension_mismatch_error(lhs, rhs);
    }
    Matrix res {lhs.row_dim(), rhs.col_dim()};
    for (int i = 1; i <= res.row_dim(); i++) {
        for (int j = 1; j <= res.col_dim(); j++) {
            double sum = 0;
            for (int k = 1; k <= lhs.col_dim(); k++) {
                sum += lhs(i, k) * rhs(k, j);
            }
            res(i, j) = sum;
        }
    }
    return res;
}

bool Matrix::equals(const Matrix& lhs, const Matrix& rhs, double eps) {
    if (!dim_eq(lhs, rhs)) return false;
    for (int i = 0; i < lhs.rows * lhs.cols; i++) {
        if (fabs(lhs.data[i] - rhs.data[i]) > eps) return false;
    }
    return true;
}

bool operator==(const Matrix& lhs, const Matrix& rhs) {
    return Matrix::equals(lhs, rhs);
}

Matrix cbind(const Matrix& m1, const Matrix& m2) {
    int m1r = m1.row_dim();
    int m2r = m2.row_dim();
    int m1c = m1.col_dim();
    int m2c = m2.col_dim();
    if (m1r != m2r) {
        throw dimension_mismatch_error(
            "Row dimension of m1: " + std::to_string(m1r) + ", " +
            "row dimension of m2: " + std::to_string(m2r));
    }

    Matrix m {m1r, m1c + m2c};

    // Copy m1 into m
    for (int i = 1; i <= m1r; i++) {
        for (int j = 1; j <= m1c; j++) {
            m(i, j) = m1(i, j);
        }
    }

    // Append m2
    for (int i = 1; i <= m2r; i++) {
        for (int j = 1; j <= m2c; j++) {
            m(i, m1c + j) = m2(i, j);
        }
    }

    return m;
}

Matrix rbind(const Matrix& m1, const Matrix& m2) {
    int m1r = m1.row_dim();
    int m2r = m2.row_dim();
    int m1c = m1.col_dim();
    int m2c = m2.col_dim();
    if (m1c != m2c) {
        throw dimension_mismatch_error(
            "Column dimension of m1: " + std::to_string(m1c) + ", " +
            "column dimension of m2: " + std::to_string(m2c));
    }

    Matrix m {m1r + m2r, m1c};

    // Copy m1 into m
    for (int i = 1; i <= m1r; i++) {
        for (int j = 1; j <= m1c; j++) {
            m(i, j) = m1(i, j);
        }
    }

    // Append m2
    for (int i = 1; i <= m2r; i++) {
        for (int j = 1; j <= m2c; j++) {
            m(m1r + i, j) = m2(i, j);
        }
    }

    return m;
}

// Returns a vector of doubles find in given string separated by spaces or 
// tabs.
// Throws std::invalid_argument if numbers aren't formatted correctly.
static std::vector<double> parse_row(const std::string& str) {
    std::vector<double> row;
    std::istringstream iss {str};
    std::string s;
    while (!iss.eof()) {
        iss >> s;
        try {
            row.push_back(std::stod(s));
        } catch (std::out_of_range) {
            throw std::invalid_argument("Value to large for double: " + s);
        } catch (std::invalid_argument) {
            throw std::invalid_argument("Invalid format: " + s);
        }
    }
    return row;
}

// True if str consists only of whitespaces, false otherwise
static bool is_whitespace(const std::string& str) {
    for (auto c : str) {
        if (c == ' ') continue;
        if (c == '\t') continue;
        if (c == '\n') continue;
        if (c == '\r') continue;
        return false;
    }
    return true;
}

std::vector<std::vector<double>> Matrix::parse(const std::string& str) {
    std::vector<std::vector<double>> rows;
    int current_index = 0;
    int cols = -1;
    while (true) {
        int newline = str.find('\n', current_index);
        newline = newline == std::string::npos ? str.length() : newline;

        if (newline <= current_index) break;

        std::string strrow = str.substr(
            current_index, newline - current_index);
        current_index = newline + 1;

        if (is_whitespace(strrow)) {
            if (cols == -1) {
                // Skip empty lines in the beginning
                continue;
            } else {
                // End of input if empty line
                break;
            }
        }

        std::vector<double> row = parse_row(strrow);
        if (cols != row.size() && cols != -1) {
            throw std::invalid_argument("Column lengths don't match: " + 
                std::to_string(cols) + ", " + std::to_string(row.size()));
        }
        cols = row.size();
        rows.push_back(row);
    }
    return rows;
}

Matrix Matrix::identity(int n) {
    if (n < 1) {
        throw std::logic_error(
            "Requested matrix dimension: " + std::to_string(n));
    }
    Matrix id {Matrix::zeros(n, n)};
    for (int i = 1; i <= n; i++) {
        id(i, i) = 1;
    }
    return id;
}

Matrix Matrix::zeros(int m, int n) {
    if (m < 1) {
        throw std::logic_error(
            "Requested matrix with " + std::to_string(m) + " rows");
    }
    if (n < 1) {
        throw std::logic_error(
            "Requested matrix with " + std::to_string(n) + " columns");
    }
    Matrix z {m, n};
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            z(i, j) = 0;
        }
    }
    return z;
}

Matrix Matrix::from_file(std::ifstream& ifs) {
    std::stringstream ss;
    ss << ifs.rdbuf();
    return Matrix {ss.str()};
}

void Matrix::to_file(std::ofstream& ofs) {
    ofs << to_string() << '\n';
}

// Returns the number of digits preceding the decimal point including minus if
// given number d is negative.
static int count_digits(double d) {
    if (isnan(d)) return 3;
    if (1.0/0 == d) return 3;
    if (-1.0/0 == d) return 4;
    int count = 0;
    while (d >= 1) {
        count++;
        d /= 10;
    }
    count = count == 0 ? 1 : count; // 0 is the first digit
    count = d < 0 ? count + 1 : count;
    return count;
}

std::string Matrix::to_string(int precision) const {
    std::string result;
    std::string fmt {'%'};
    int max_digits = 0;
    for (int i = 0; i < rows*cols; i++) {
        int d = count_digits(data[i]);
        max_digits = max_digits >= d ? max_digits : d;
    }
    fmt += std::to_string(1 + max_digits + precision) + '.' + 
        std::to_string(precision) + 'f';
    const char* cfmt = fmt.c_str();
    std::unique_ptr<char> buff {new char[max_digits + precision + 2]};
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            sprintf(buff.get(), cfmt, this->operator()(i, j));
            result += buff.get();
            result += ' ';
            result += ' ';
        }
        result.pop_back();
        result += '\n';
    }
    result.pop_back();
    return result;
}

Matrix solve(const Matrix& a, const Matrix& b, double eps) {
    LU lu {a, true, eps};
    return lu.solve(b, eps);
}

Matrix inv(const Matrix& m, double eps) {
    LU lu {m, true, eps};
    return lu.inverse(eps);
}

double det(const Matrix& m) {
    try {
        LU lu {m, true, 0};
        return lu.determinant();
    } catch (singular_matrix_error) {
        return 0;
    }
}

//
// LU decomposition
//

LU::LU(const Matrix& m, bool exchange_rows, double eps) 
:lu {m} {
    if (m.row_dim() != m.col_dim()) {
        throw matlib_error(
            "LU decomposition is allowed only for square matrices");
    }
    if (exchange_rows) {
        // +1 size for easier indexing
        p = new int[m.row_dim()+1];
        lup_decomposition(eps);
    } else {
        lu_decomposition(eps);
    }
}

LU::~LU() {
    if (p != nullptr) delete[] p;
}

void LU::lu_decomposition(double eps) {
    int n = lu.row_dim();
    // iterate over elements on main diagonal
    for (int d = 1; d <= n-1; d++) {
        // iterate over rows
        for (int i = d+1; i <= n; i++) {
            if (abs(lu(d, d)) < eps) singular_error(lu(d, d));
            lu(i, d) /= lu(d, d);
            // iterate over elements of row i
            for (int j = d+1; j <= n; j++) {
                lu(i, j) -= lu(i, d) * lu(d, j);
            }
        }
    }
}

static void swap(int& m, int& n) {
    int tmp = m;
    m = n;
    n = tmp;
}

static void swap(double& x, double& y) {
    double tmp = x;
    x = y;
    y = tmp;
}

void LU::lup_decomposition(double eps) {
    int n = lu.row_dim();
    // initialize p[0] to 0 for counting row exchanges
    for (int i = 0; i <= n; i++) p[i] = i;
    for (int d = 1; d <= n-1; d++) {
        int pivot = d;
        for (int i = d+1; i <= n; i++) {
            if (abs(lu(i, d)) > abs(lu(pivot, d))) {
                pivot = i;
            }
        }
        swap(p[d], p[pivot]);
        if (d != pivot) {
            p[0]++;
            for (int j = 1; j <= n; j++) {
                swap(lu(d, j), lu(pivot, j));
            }
        }

        for (int i = d+1; i <= n; i++) {
            if (abs(lu(d, d)) < eps) singular_error(lu(d, d));
            lu(i, d) /= lu(d, d);
            for (int j = d+1; j <= n; j++) {
                lu(i, j) -= lu(i, d) * lu(d, j);
            }
        }
    }
}

Matrix LU::forward_subs(const Matrix& b) {
    Matrix y {b};
    int cols = b.col_dim();
    int n = lu.row_dim();
    for (int col = 1; col <= cols; col++) {
        for (int j = 1; j <= n-1; j++) {
            for (int i = j+1; i <= n; i++) {
                y(i, col) -= lu(i, j) * y(j, col);
            }
        }
    }
    return y;
}

Matrix LU::backward_subs(const Matrix& y, double eps) {
    Matrix x {y};
    int cols = y.col_dim();
    int n = lu.row_dim();
    for (int col = 1; col <= cols; col++) {
        for (int j = n; j >= 1; j--) {
            if (abs(lu(j, j)) < eps) singular_error(lu(j, j));
            x(j, col) /= lu(j, j);
            for (int i = 1; i <= j-1; i++) {
                x(i, col) -= lu(i, j) * x(j, col);
            }
        }
    }
    return x;
}

Matrix LU::l() {
    Matrix l {lu};
    int n = lu.row_dim();
    for (int i = 1; i <= n; i++) l(i, i) = 1;
    for (int i = 1; i <= n-1; i++) {
        for (int j = i+1; j <= n; j++) {
            l(i, j) = 0;
        }
    }
    return l;
}

Matrix LU::u() {
    Matrix u {lu};
    int n = lu.row_dim();
    for (int i = 2; i <= n; i++) {
        for (int j = 1; j <= i-1; j++) {
            u(i, j) = 0;
        }
    }
    return u;
}

// Returns a row permutation of m. Row i of original matrix will bi row p[i] of
// the returned matrix.
static Matrix permute(const Matrix& m, const int* p) {
    int rows = m.row_dim();
    int cols = m.col_dim();
    Matrix permuted {rows, cols};
    for (int i = 1; i <= rows; i++) {
        int row_loc = p[i];
        for (int j = 1; j <= cols; j++) {
            permuted(i, j) = m(row_loc, j);
        }
    }
    return permuted;
}

Matrix LU::solve(const Matrix& m, double eps) {
    if (m.row_dim() != lu.row_dim()) {
        throw dimension_mismatch_error(
            "Cannot solve a system of matrix dimension " + dim_str(lu) + 
            " for a matrix with " + std::to_string(m.row_dim()) + "rows");
    }
    if (p == nullptr) {
        // lu decomposition
        Matrix y {forward_subs(m)};
        Matrix x {backward_subs(y, eps)};
        return x;
    } else {
        // lup decomposition
        Matrix mp = std::move(permute(m, p));
        Matrix y {forward_subs(mp)};
        Matrix x {backward_subs(y, eps)};
        return x;
    }
}

Matrix LU::inverse(double eps) {
    return solve(Matrix::identity(lu.row_dim()), eps);
}

double LU::determinant() {
    int n = lu.row_dim();
    double d = 1;
    for (int i = 1; i <= n; i++) {
        d *= lu(i, i);
    }
    if (p != nullptr) {
        d *= 1 - 2 * (p[0] % 2);
    }
    return d;
}

std::string dimension_mismatch_error::error_string
(const Matrix& m1, const Matrix& m2) {
    return "Dimensions: (" + dim_str(m1) + "), (" + dim_str(m2) + ")";
}

//
// UTILITY FUNCTIONS
//

// True if m1 and m2 have same dimensions, false otherwise.
static bool dim_eq(const Matrix& m1, const Matrix& m2) {
    return m1.row_dim() == m2.row_dim() && m1.col_dim() == m2.col_dim();
}

// Returns the string of form "row_dim, col_dim".
static std::string dim_str(const Matrix& m) {
    return std::to_string(m.row_dim()) + ", " + std::to_string(m.col_dim());
}

// Throws singular_matrix_error with message "Value near zero encountered: 
// ${value}".
static void singular_error(double value) {
    std::string error_message {"Value near zero encountered: "};
    char buff[20];
    sprintf(buff, "%.3e", value);
    error_message += buff;
    throw singular_matrix_error(error_message);
}

} // namespace matlib

std::ostream& operator<<(std::ostream& out, matlib::Matrix& m) {
    return out << m.to_string();
}

std::ostream& operator<<(std::ostream& out, matlib::Matrix&& m) {
    return out << m;
}

std::ostream& operator<<(std::ostream& out, matlib::Matrix::slice& s) {
    return out << '[' << s.from << ',' << ' ' << s.to << ']';
}

std::ostream& operator<<(std::ostream& out, matlib::LU& lu) {
    matlib::Matrix l {lu.l()};
    matlib::Matrix u {lu.u()};
    out << l << std::endl << std::endl;
    out << u;
    return out;
}
