/**
    NCG++
    gamma_matrix.hpp

    Purpose:
    Matrix class especially for
    gamma matrices

    @author Peter Labus
    @version 0.1
    10.05.2017
*/

#pragma once

#include <iostream>
#include <vector>
#include <complex>


class GammaMatrix
{
  public:

    GammaMatrix(void) = delete;
    explicit GammaMatrix(const int size);
    explicit GammaMatrix(std::initializer_list< std::complex<int> > const& list);

    std::complex<int> operator()(const int row, const int col) const { return M_[row*size_+col]; }
    std::complex<int>& operator()(const int row, const int col) { return M_[row*size_+col]; }
    GammaMatrix operator*=(GammaMatrix const& other);
    GammaMatrix operator*=(std::complex<int> c);

    friend std::ostream& operator<<(std::ostream& os, GammaMatrix const& A);

    int size(void) const { return size_; }

  private:

    const int size_;
    std::vector< std::complex<int>  > M_;

    void print_(void) const;
};

// These are addition, substraction, multiplication and outer (tensor) multiplication
GammaMatrix operator+(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix operator-(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix operator*(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix operator*(std::complex<int> c, GammaMatrix const& A);
GammaMatrix operator%(GammaMatrix const& A, GammaMatrix const& B);

GammaMatrix commutator(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix anticommutator(GammaMatrix const& A, GammaMatrix const& B);


struct PauliMatrices
{
  PauliMatrices(void);
  const GammaMatrix sigma1;
  const GammaMatrix sigma2;
  const GammaMatrix sigma3;
};

GammaMatrix Unity(const int d);
