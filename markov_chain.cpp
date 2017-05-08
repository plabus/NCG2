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

#include "MarkovChain.hpp"

explicit MarkovChain::MarkovChain(SimulationParameters const& sim_params)
{
  // TODO:
  // Implement me!
}


~MarkovChain::MarkovChain(void)
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
    D_old = D_new;
  }
  else                                                     // reject
  {
    // Nothing to do
  }
}
