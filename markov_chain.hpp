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

#include "algorithms.hpp"
#include "model_parameters.hpp"
#include "action_parameters.hpp"
#include "rng.hpp"

class MarkovChain
{
  public:

    MarkovChain(void) = delete;
    explicit MarkovChain(
        const Algorithm algorithm,
        const ModelParameters model_params,
        const ActionParameters action_params,
        RNG& rng,
        const IOHandler& io_handler
    );
    ~MarkovChain(void);

    void generate_configurations(int number);


  private:

    const Algorithm algorithm_;
    const ModelParameters model_params_;
    const ActionParameters action_params_;
    RNG& rng_;
    const IOHandler& io_handler_;

    Dirac_Matrix D_old_;
    Dirac_Matrix D_new_;
    Action S_old_;
    Action S_new_;
    const TimeEvolutionOperator T_;

    void get_next_configuration(void);
};
