#ifndef hhAnalysis_bbwwMEM_MEMIntegratorBase_h
#define hhAnalysis_bbwwMEM_MEMIntegratorBase_h

/** \class MEMIntegratorBase
 *
 * Base class for VEGAS and VAMP integration classes.
 *
 * \author Christian Veelken, NICPB Tallinn
 *
 */

#include <iostream>
#include <assert.h>

namespace mem
{
  class MEMIntegratorBase
  {
   public:
    MEMIntegratorBase() {}
    virtual ~MEMIntegratorBase() {}

    /// compute integral of function g 
    /// the points xl and xh represent the lower left and upper right corner of a Hypercube in d-dimensional integration space 
    typedef double (*gPtr_C)(double*, size_t, void*);
    virtual void integrate(gPtr_C, const double*, const double*, unsigned, double&, double&) { assert(0); }
    typedef double (*gPtr_Fortran)(double**, size_t, void**);
    virtual void integrate(gPtr_Fortran, const double*, const double*, unsigned, double&, double&) { assert(0); }

    virtual void print(std::ostream&) const {}
  };
}

#endif

