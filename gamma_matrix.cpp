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
#include <algorithm>
#include <cassert>
#include <cmath>
#include "gamma_matrix.hpp"


GammaMatrix::GammaMatrix(const int size)
  :
size_(size),
M_(size*size, 0)
{}

GammaMatrix::GammaMatrix(std::initializer_list< std::complex<int> > const& list)
  :
size_(sqrt(list.size())),
M_(list)
{}

std::ostream& operator<<(std::ostream& os, const GammaMatrix& A)
{
  A.print_();
  return os;
}

GammaMatrix& GammaMatrix::operator=(GammaMatrix const& other)
{
  assert(this->size()==other.size() && "GammaMatrix Copy Constructor: ERROR: Matrices have different sizes");

  if(&other == this) return *this;

  std::copy(other.M_.begin(), other.M_.end(), M_.begin());
  return *this;
}

void GammaMatrix::print_(void) const
{
  for(auto row = 0; row < size_; ++row)
  {
    for(auto col = 0; col < size_; ++col)
    {
      std::cout << " " << (*this)(row,col);
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

GammaMatrix operator*(std::complex<int> c, GammaMatrix const& A)
{
  const int size = A.size();
  GammaMatrix B(size);

  for(auto i = 0; i < size; ++i)
    for(auto j = 0; j < size; ++j)
        B(i,j) = c * A(i,j);

  return B;
}

GammaMatrix operator/(GammaMatrix const& A, int c)
{
  const int size = A.size();
  GammaMatrix B(size);

  for(auto i = 0; i < size; ++i)
    for(auto j = 0; j < size; ++j)
        B(i,j) = A(i,j) / c;

  return B;
}

GammaMatrix& GammaMatrix::operator*=(GammaMatrix const& other)
{
  assert(size_==other.size() && "GammaMatrix Multiplication Assignment: ERROR: Matrices have different sizes");
  *this = *this * other;
  return *this;
}

GammaMatrix& GammaMatrix::operator*=(std::complex<int> c)
{
  *this = c * (*this);
  return *this;
}

GammaMatrix operator%(GammaMatrix const& A, GammaMatrix const& B)
{
  const int size_A = A.size();
  const int size_B = B.size();
  GammaMatrix C(size_A * size_B);

  for(auto i = 0; i < size_A; ++i)
    for(auto j = 0; j < size_A; ++j)
      for(auto ii = 0; ii < size_B; ++ii)
        for(auto jj = 0; jj < size_B; ++jj)
          C(i*size_B+ii,j*size_B+jj) += A(i,j) * B(ii,jj);

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
sigma1(GammaMatrix(
      { {0,0}, {1,0},
        {1,0}, {0,0} }
      )),
sigma2(GammaMatrix(
      { {0,0}, {0,-1},
        {0,1}, {0,0} }
      )),
sigma3(GammaMatrix(
      { {1,0}, {0,0},
        {0,0}, {-1,0} }
      ))
{}

GammaMatrix Unity(const int d)
{
  auto ret = GammaMatrix(d);
  for(auto i = 0; i < d; ++i)
    ret(i,i) = 1;
  return ret;
}
