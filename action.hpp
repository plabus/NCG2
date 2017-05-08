/**
    NCG++
    action.hpp

    Purpose:
    Provide an interface and wrapper
    to unify the action parameters and
    an instance of a Dirac matrix.

    @author Peter Labus
    @version 0.1
    08.05.2017
*/

#pragma once

#include "precisions.hpp"

class Action
{
  public:

    Action(void) = delete;
    explicit Action(
      const FloatingPointPrecision prec,
      const ActionParameters action_params,
      const DiracMatrix& D
    ) :
      prec_(prec),
      action_params_(action_params),
      D_(D)
    {};

    virtual double operator()(void) const; //!< Calculates the action of DiracMatrix with given parameters.


  private:

    const FloatingPointPrecision prec_;
    const ActionParameters action_params_;
    DiracMatrix& D_;
};
