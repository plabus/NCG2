/**
    NCG++
    markov_chain.cpp

    Purpose:
    Generates Markov Chains for matrix models in
    non-commutative geometries.

    @author Peter Labus
    @version 0.1
    08.05.2017
*/

#include "markov_chain.hpp"

MarkovChain::MarkovChain(
    const Algorithm algorithm,
    const ModelParameters model_params,
    const ActionParameters action_params,
    RNG& rng,
    const IOHandler& io_handler
) :
    algorithm_(algorithm),
    model_params_(model_params),
    action_params_(action_params),
    rng_(rng),
    io_handler_(io_handler)
{
  // TODO:
  // Implement me!
}


MarkovChain::~MarkovChain(void)
{
  // TODO:
  // Implement me!
}


void MarkovChain::generate_configurations(int number)
{
  for( auto iter = 0; iter != number; ++iter )
  {
    get_next_configuration();
    // TODO: I/O Handling
    // if( write_out(iter) )
    // {
    //   io_handler_.write();
    // }
  }
}


void MarkovChain::get_next_configuration(void)
{
  T_(D_old_, D_new_);
  const double delta_S = S_old_ - S_new_;
  if( delta_S <= 0 || exp(-delta_S) > rng_.uniform(0, 1) ) // accept
  {
    // FIXME:
    // Can we have a swap here??
    D_old_ = D_new_;
  }
  else                                                     // reject
  {
    // Nothing to do
  }
}
