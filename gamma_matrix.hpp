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

    explicit GammaMatrix(const int size)
      :
    size_(size),
    M_(size, std::vector< std::complex<int8_t> >(size, 0))
    {}


    std::complex<int8_t> operator()(const int row, const int col) const { return M_[row][col]; }
    std::complex<int8_t>& operator()(const int row, const int col) { return M_[row][col]; }
    int size(void) const { return size_; }
    void print(void) const
    {
      for(auto row = 0; row < M_.size(); ++row)
      {
        for(auto col = 0; col < M_.size(); ++col)
        {
          std::cout << " "
                    << static_cast<std::complex<int>>((*this)(row,col));
        }
        std::cout << std::endl;
      }
    }

  private:

    const int size_;
    std::vector< std::vector< std::complex<int8_t> > > M_;
};

GammaMatrix operator+(GammaMatrix const& A, GammaMatrix const& B)
{
  assert(A.size()==B.size() && "GammaMatrix Addition: ERROR: Matrices have different sizes");
  const int size = A.size();
  GammaMatrix C(size);

  for(auto i = 0; i < size; ++i)
    for(auto j = 0; j < size; ++j)
      C(i,j) = A(i,j) + B(i,j);

  return C;
}

GammaMatrix operator-(GammaMatrix const& A, GammaMatrix const& B)
{
  assert(A.size()==B.size() && "GammaMatrix Substraction: ERROR: Matrices have different sizes");
  const int size = A.size();
  GammaMatrix C(size);

  for(auto i = 0; i < size; ++i)
    for(auto j = 0; j < size; ++j)
      C(i,j) = A(i,j) - B(i,j);

  return C;
}

GammaMatrix operator*(GammaMatrix const& A, GammaMatrix const& B)
{
  assert(A.size()==B.size() && "GammaMatrix Multiplication: ERROR: Matrices have different sizes");
  const int size = A.size();
  GammaMatrix C(size);

  for(auto i = 0; i < size; ++i)
    for(auto j = 0; j < size; ++j)
      for(auto k = 0; k < size; ++k)
        C(i,j) += A(i,k) * B(k,j);

  return C;
}

GammaMatrix commutator(GammaMatrix const& A, GammaMatrix const& B)
{
  assert(A.size()==B.size() && "GammaMatrix Commutator: ERROR: Matrices have different sizes");
  return A*B-B*A;
}

GammaMatrix anticommutator(GammaMatrix const& A, GammaMatrix const& B)
{
  assert(A.size()==B.size() && "GammaMatrix Anticommutator: ERROR: Matrices have different sizes");
  return A*B+B*A;
}
