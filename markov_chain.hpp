/**
    NCG++
    markov_chain.hpp

    Purpose:
    Generates Markov Chains for matrix models in
    non-commutative geometries.

    @author Peter Labus
    @version 0.1
    08.05.2017
*/

#pragma once

class MarkovChain
{
  public:

    explicit MarkovChain(SimulationParams sparams);

    void generate_configuration(int number);


  private:

    Algorithm algorithm_;
    ModelParameters model_params_;
    ActionParameters action_params_;
    RNG rng_;
    IOHandler io_handler_;

    Dirac_Matrix D_old_;
    Dirac_Matrix D_new_;
    Action S_old_;
    Action S_new_;
    TimeEvolutionOperator T_;

    void get_next_configuration(void);
};
