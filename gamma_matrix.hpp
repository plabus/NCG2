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
#include <cstdint>
#include <cassert>


class GammaMatrix
{
  public:

    GammaMatrix(void) = delete;
    explicit GammaMatrix(const int size);

    friend std::ostream& operator<<(std::ostream& os, GammaMatrix const& A);

    std::complex<int8_t> operator()(const int row, const int col) const { return M_[row][col]; }
    std::complex<int8_t>& operator()(const int row, const int col) { return M_[row][col]; }
    int size(void) const { return size_; }

  private:

    const int size_;
    std::vector< std::vector< std::complex<int8_t> > > M_;

    void print_(void) const;
};

GammaMatrix operator+(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix operator-(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix operator*(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix commutator(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix anticommutator(GammaMatrix const& A, GammaMatrix const& B);


struct PauliMatrices
{
  PauliMatrices(void);
  GammaMatrix sigma1;
  GammaMatrix sigma2;
  GammaMatrix sigma3;
};
