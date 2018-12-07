#ifndef hhAnalysis_bbwwMEM_MEMIntegratorVAMP_h
#define hhAnalysis_bbwwMEM_MEMIntegratorVAMP_h

/** \class MEMIntegratorVAMP
 *
 * Interface to "Vegas AMPlified: Anisotropy, Multi-channel sampling and Parallelization" (VAMP) integration algorithm.
 *
 * The VAMP algorithm is documented in:
 *  [1] "Vegas revisited: Adaptive Monte Carlo integration beyond factorization",
 *      T. Ohl, J. Comput. Phys. 120 (1999) 13.
 *  [2] https://whizard.hepforge.org/vamp.pdf
 *
 * \author Christian Veelken, NICPB Tallinn
 *
 */

#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorBase.h"

#include <vector>
#include <string>
#include <iostream>

namespace mem
{
  class MEMIntegratorVAMP : public MEMIntegratorBase
  {
   public:
    MEMIntegratorVAMP(unsigned, unsigned);
    ~MEMIntegratorVAMP();

    void integrate(MEMIntegratorBase::gPtr_Fortran, const double*, const double*, unsigned, double&, double&);

    void print(std::ostream&) const {}

   protected:
    void setIntegrand(MEMIntegratorBase::gPtr_Fortran, const double*, const double*, unsigned);

    MEMIntegratorBase::gPtr_Fortran integrand_;

    unsigned numCallsGridOpt_;
    unsigned numCallsIntEval_;
    unsigned numDimensions_;

    /// lower and upper boundary of integration region
    mutable double* xl_;
    mutable double* xu_;
  };
}

#endif

