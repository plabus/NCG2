/**
    NCG++
    time_evolution_operator.hpp

    Purpose:
    An abstract class that provides an operator
    to give a new Markov Chain element for Dirac Matrix D
    given an old state.
    An implementation will depend on the algorithm used.

    @author Peter Labus
    @version 0.1
    09.05.2017
*/

#pragma once

template<typename FT>
class TimeEvolutionOperator
{
  public:

    virtual void operator()(DiracMatrix<FT> const& D_old, DiracMatrix<FT>& D_new) const;


  private:

    const Algorithm algorithm_;
};
