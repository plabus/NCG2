
/**
    NCG++
    action_parameters.hpp

    Purpose:
    Structure all parameters needed for
    an action of the form:
      S = Tr ( g_2 * D^2 + g_4 * D^4 + ... )

    @author Peter Labus
    @version 0.1
    08.05.2017
*/

#pragma once

class ActionParameters
{
  public:

    explicit ActionParameters(
        double G2,
        double G4
    ) :
      g2_(G2),
      g4_(G4)
    {}

    double g2(void) const { return g2_; }
    double g4(void) const { return g4_; }


  private:

    const double g2_;
    const double g4_;
};
