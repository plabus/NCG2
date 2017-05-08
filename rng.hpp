/**
    NCG++
    rng.hpp

    Purpose:
    Provide an interface for a
    Random Number Generator with
    all the functionality needed.

    @author Peter Labus
    @version 0.1
    08.05.2017
*/


// FIXME:
// Make this a template over the floating point type
class RNG
{
  public:

    virtual double uniform(double lower, double uppper);
    virtual double signed_uniform(double symetric_boundary);
    virtual int uniform_int(int lower, int upper);
};
