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


template<typename FT>
class Action
{
  public:

    Action(void) = delete;
    explicit Action(
      const ActionParameters action_params,
      const DiracMatrix<FT>& D
    ) :
      action_params_(action_params),
      D_(D)
    {};

    virtual double operator()(void) const; //!< Calculates the action of DiracMatrix with given parameters.


  private:

    const ActionParameters action_params_;
    DiracMatrix<FT>& D_;
};
