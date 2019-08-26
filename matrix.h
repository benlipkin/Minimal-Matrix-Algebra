#ifndef MATRIX_H
#define MATRIX_H

class Matrix {
private:
    double **m_array;
    int m_rows;
    int m_cols;
public:
    Matrix(int rows, int cols);
    Matrix(const Matrix& src);
    ~Matrix();
    Matrix& operator=(const Matrix& src);
    double& operator()(int row, int col);
    double operator()(int row, int col) const;
    friend Matrix operator+(const Matrix &mat1, const Matrix &mat2);
    friend Matrix operator-(const Matrix &mat1, const Matrix &mat2);
    friend Matrix operator*(const Matrix &mat1, const Matrix &mat2);
    friend Matrix operator*(const Matrix &mat, double scl);
    friend std::ostream& operator<<(std::ostream &out, const Matrix &mat);
    int size(int dim) const;
    Matrix get(int rmin, int rmax, int cmin, int cmax) const;
    Matrix& set(Matrix mat, int rmin, int rmax, int cmin, int cmax);
    Matrix& fill(double val);
    Matrix& rand(double min, double max);
    Matrix& randn(double mean, double sigma);
    Matrix transpose();
};

#endif
