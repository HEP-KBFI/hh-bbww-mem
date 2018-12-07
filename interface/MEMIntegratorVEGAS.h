#ifndef hhAnalysis_bbwwMEM_MEMIntegratorVEGAS_h
#define hhAnalysis_bbwwMEM_MEMIntegratorVEGAS_h

/** \class MEMIntegratorVEGAS
 *
 * Interface to VEGAS integration algorithm.
 *
 * The VEGAS algorithm is documented in:
 *  [1] "A New Algorithm for Adaptive Multidimensional Integration",
 *      G.P. Lepage, J. Comput. Phys. 27 (1978) 192.
 *
 * \author Christian Veelken, NICPB Tallinn
 *
 */

#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorBase.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <iostream>

namespace mem
{
  class MEMIntegratorVEGAS : public MEMIntegratorBase
  {
   public:
    MEMIntegratorVEGAS(unsigned, unsigned, double, unsigned);
    ~MEMIntegratorVEGAS();

    void integrate(MEMIntegratorBase::gPtr_C, const double*, const double*, unsigned, double&, double&);

    void print(std::ostream&) const {}

   protected:
    void setIntegrand(MEMIntegratorBase::gPtr_C, const double*, const double*, unsigned);

    MEMIntegratorBase::gPtr_C integrand_;

    gsl_monte_function* vegasIntegrand_;
    gsl_monte_vegas_state* vegasWorkspace_;
    mutable gsl_rng* vegasRnd_;
    unsigned numCallsGridOpt_;
    unsigned numCallsIntEval_;
    double maxChi2_;
    unsigned maxIntEvalIter_;
    double precision_;
    unsigned numDimensions_;

    /// lower and upper boundary of integration region
    mutable double* xl_;
    mutable double* xu_;
  };
}

#endif

