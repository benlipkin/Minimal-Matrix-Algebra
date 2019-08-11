#include <iostream>
#include <iomanip>
#include <random>
#include "matrix.h"

static std::default_random_engine prng;

Matrix::Matrix(int rows, int cols): m_rows(rows), m_cols(cols) {
    assert(rows>0 && cols>0);
    m_array=new double*[m_rows];
    for (int row=0; row<m_rows; row++) {
        m_array[row]=new double[m_cols];
    }
    this->fill(0);
}

Matrix::Matrix(const Matrix& src) {
    m_rows=src.m_rows;
    m_cols=src.m_cols;
    m_array=new double*[m_rows];
    for (int row=0; row<m_rows; row++) {
        m_array[row]=new double[m_cols];
        for (int col=0; col<m_cols; col++) {
            m_array[row][col]=src.m_array[row][col];
        }
    }
}

Matrix::~Matrix() {
    for (int row=0; row<m_rows; row++) {
        delete[] m_array[row];
        m_array[row]=nullptr;
    }
    delete[] m_array;
    m_array=nullptr;
}

Matrix& Matrix::operator=(const Matrix& src) {
    if (this==&src) {
        return *this;
    }
    for (int row=0; row<m_rows; row++) {
        delete[] m_array[row];
        m_array[row]=nullptr;
    }
    delete[] m_array;
    m_rows=src.m_rows;
    m_cols=src.m_cols;
    m_array=new double*[m_rows];
    for (int row=0; row<m_rows; row++) {
        m_array[row]=new double[m_cols];
        for (int col=0; col<m_cols; col++) {
            m_array[row][col]=src.m_array[row][col];
        }
    }
    return *this;
}

double& Matrix::operator()(int row, int col) {
    assert(row>=0 && row<m_rows && col>=0 && col<m_cols);
    return m_array[row][col];
}

double Matrix::operator()(int row, int col) const {
    assert(row>=0 && row<m_rows && col>=0 && col<m_cols);
    return m_array[row][col];
}

Matrix operator+(const Matrix &mat1, const Matrix &mat2) {
    assert(mat1.m_rows==mat2.m_rows && mat1.m_cols==mat2.m_cols);
    Matrix out(mat1.m_rows,mat1.m_cols);
    for (int row=0; row<mat1.m_rows; row++) {
        for (int col=0; col<mat1.m_cols; col++) {
            out.m_array[row][col]+=mat1.m_array[row][col];
            out.m_array[row][col]+=mat2.m_array[row][col];
        }
    }
    return out;
}

Matrix operator-(const Matrix &mat1, const Matrix &mat2) {
    assert(mat1.m_rows==mat2.m_rows && mat1.m_cols==mat2.m_cols);
    Matrix out(mat1.m_rows,mat1.m_cols);
    for (int row=0; row<mat1.m_rows; row++) {
        for (int col=0; col<mat1.m_cols; col++) {
            out.m_array[row][col]+=mat1.m_array[row][col];
            out.m_array[row][col]-=mat2.m_array[row][col];
        }
    }
    return out;
}

Matrix operator*(const Matrix &mat1, const Matrix &mat2) {
    assert(mat1.m_cols==mat2.m_rows);
    Matrix out(mat1.m_rows,mat2.m_cols);
    for (int row=0; row<mat1.m_rows; row++) {
        for (int col=0; col<mat2.m_cols; col++) {
            for (int idm=0; idm<mat1.m_cols; idm++) {
                out.m_array[row][col]+=mat1.m_array[row][idm]*mat2.m_array[idm][col];
            }
        }
    }
    return out;
}

Matrix operator*(const Matrix &mat, double scl) {
    Matrix out(mat.m_rows,mat.m_cols);
    for (int row=0; row<mat.m_rows; row++) {
        for (int col=0; col<mat.m_cols; col++) {
            out.m_array[row][col]=scl*mat.m_array[row][col];
        }
    }
    return out;
}

std::ostream& operator<<(std::ostream &out, const Matrix &mat) {
    std::cout << std::scientific << std::setprecision(2);
    for (int row=0; row<mat.m_rows; row++) {
        for (int col=0; col<mat.m_cols; col++) {
            std::cout << std::setw(10) << std::right <<  mat.m_array[row][col] << ' ';
        }
        std::cout << '\n';
    }
    return out;
}

int Matrix::size(int dim) const {
    assert(dim==0 || dim==1);
    return (dim) ? m_cols : m_rows;
}

Matrix Matrix::get(int rmin, int rmax, int cmin, int cmax) const {
    assert(rmin>=0 && cmin>=0 && rmax<=m_rows && cmax<=m_cols);
    Matrix out(rmax-rmin,cmax-cmin);
        for (int row=rmin; row<rmax; row++) {
            for (int col=cmin; col<cmax; col++) {
                out.m_array[row-rmin][col-cmin]=m_array[row][col];
            }
        }
    return out;
}

Matrix& Matrix::set(Matrix mat, int rmin, int rmax, int cmin, int cmax) {
    assert(mat.m_rows==rmax-rmin && mat.m_cols==cmax-cmin);
    for (int row=rmin; row<rmax; row++) {
        for (int col=cmin; col<cmax; col++) {
            m_array[row][col]=mat.m_array[row-rmin][col-cmin];
        }
    }
    return *this;
}

Matrix& Matrix::fill(double val) {
    for (int row=0; row<m_rows; row++) {
        for (int col=0; col<m_cols; col++) {
            m_array[row][col]=val;
        }
    }
    return *this;
}

Matrix& Matrix::rand(double min, double max) {
    std::uniform_real_distribution<double> unif(min,max);
    for (int row=0; row<m_rows; row++) {
        for (int col=0; col<m_cols; col++) {
            m_array[row][col]=unif(prng);
        }
    }
    return *this;
}

Matrix& Matrix::randn(double mu, double sigma) {
    std::normal_distribution<double> norm(mu,sigma);
    for (int row=0; row<m_rows; row++) {
        for (int col=0; col<m_cols; col++) {
            m_array[row][col]=norm(prng);
        }
    }
    return *this;
}

Matrix Matrix::transpose() {
    Matrix out(m_cols,m_rows);
    for (int row=0; row<m_rows; row++) {
        for (int col=0; col<m_cols; col++) {
            out.m_array[col][row]=m_array[row][col];
        }
    }
    return out;
}
