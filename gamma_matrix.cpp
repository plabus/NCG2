/**
    NCG++
    gamma_matrix.cpp

    Purpose:


    @author Peter Labus
    @version 0.1
    10.05.2017
*/

#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include "gamma_matrix.hpp"


GammaMatrix::GammaMatrix(const int size)
  :
size_(size),
M_(size, std::vector< std::complex<int> >(size, 0))
{}

std::ostream& operator<<(std::ostream& os, const GammaMatrix& A)
{
  A.print_();
  return os;
}

void GammaMatrix::print_(void) const
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



PauliMatrices::PauliMatrices(void)
  :
sigma1( GammaMatrix(2) ),
sigma2( GammaMatrix(2) ),
sigma3( GammaMatrix(2) )
{
  std::complex<int> I(0,1);

  sigma1(0,1) = 1;
  sigma1(1,0) = 1;

  sigma2(0,1) = -I;
  sigma2(1,0) = I;

  sigma3(0,0) = 1;
  sigma3(1,1) = -1;
}
