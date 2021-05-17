//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_sm_ttbar2WbWb_Wp2jj_Wn2vl_H
#define MG5_sm_ttbar2WbWb_Wp2jj_Wn2vl_H

#include <complex>
#include <vector>

#include "hhAnalysis/bbwwMEM/interface/mg5/parameters/Parameters_sm_ttbar2WbWb_Wp2jj_Wn2vl.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > t t~ WEIGHTED=2 @1
// *   Decay: t > w+ b WEIGHTED=2
// *     Decay: w+ > u d~ WEIGHTED=2
// *   Decay: t~ > w- b~ WEIGHTED=2
// *     Decay: w- > e- ve~ WEIGHTED=2
// Process: g g > t t~ WEIGHTED=2 @1
// *   Decay: t > w+ b WEIGHTED=2
// *     Decay: w+ > c s~ WEIGHTED=2
// *   Decay: t~ > w- b~ WEIGHTED=2
// *     Decay: w- > e- ve~ WEIGHTED=2
// Process: g g > t t~ WEIGHTED=2 @1
// *   Decay: t > w+ b WEIGHTED=2
// *     Decay: w+ > u d~ WEIGHTED=2
// *   Decay: t~ > w- b~ WEIGHTED=2
// *     Decay: w- > mu- vm~ WEIGHTED=2
// Process: g g > t t~ WEIGHTED=2 @1
// *   Decay: t > w+ b WEIGHTED=2
// *     Decay: w+ > c s~ WEIGHTED=2
// *   Decay: t~ > w- b~ WEIGHTED=2
// *     Decay: w- > mu- vm~ WEIGHTED=2
//--------------------------------------------------------------------------

class mg5_sm_ttbar2WbWb_Wp2jj_Wn2vl
{
  public:

    // Constructor.
    mg5_sm_ttbar2WbWb_Wp2jj_Wn2vl()
    {
      for(std::size_t i = 0; i < nprocesses; ++i)
      {
        jamp2[i] = nullptr;
      }
    }

    // Destructor.
    virtual ~mg5_sm_ttbar2WbWb_Wp2jj_Wn2vl()
    {
      for(std::size_t i = 0; i < nprocesses; ++i)
      {
        if(jamp2[i])
        {
          delete [] jamp2[i];
          jamp2[i] = nullptr;
        }
      }
    }

    // Initialize process.
    virtual void
    initProc(const std::string & param_card_name);

    // Calculate flavour-independent parts of cross section.
    virtual void
    sigmaKin();

    // Evaluate sigmaHat(sHat).
    virtual double
    sigmaHat();

    // Info on the subprocess.
    virtual std::string
    name() const
    {
      return "g g > u d~ b e- ve~ b~ (sm)";
    }

    const std::vector<double> &
    getMasses() const
    {
      return mME;
    }

    // Get and set momenta for matrix element evaluation
    std::vector<double *>
    getMomenta()
    {
      return p;
    }

    void
    setMomenta(std::vector<double *> & momenta)
    {
      p = momenta;
    }

    void
    setInitial(int inid1,
               int inid2)
    {
      id1 = inid1;
      id2 = inid2;
    }

    // Get matrix element vector
    const double *
    getMatrixElements() const
    {
      return matrix_element;
    }

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 8;
    static const int nprocesses = 1;

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void
    calculate_wavefunctions(const int perm[],
                            const int hel[]);

    static const int nwavefuncs = 15;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 3;
    std::complex<double> amp[namplitudes];

    double matrix_1_gg_ttx_t_wpb_wp_udx_tx_wmbx_wm_emvex();

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses];

    // Color flows, used when selecting color
    double * jamp2[nprocesses];

    // Pointer to the model parameters
    Parameters_sm_ttbar2WbWb_Wp2jj_Wn2vl * pars;

    // vector with external particle masses
    std::vector<double> mME;

    // vector with momenta (to be changed each event)
    std::vector<double *> p;
    // Initial particle ids
    int id1, id2;

};


#endif  // MG5_sm_ttbar2WbWb_Wp2jj_Wn2vl_H
