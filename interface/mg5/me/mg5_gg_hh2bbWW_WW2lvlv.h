//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_BSM_gg_hh2bbWW_WW2lvlv_H
#define MG5_Sigma_BSM_gg_hh2bbWW_WW2lvlv_H

#include <complex>
#include <vector>

#include "hhAnalysis/bbwwMEM/interface/mg5/parameters/Parameters_BSM_gg_hh2bbWW_WW2lvlv.h"

using namespace std;

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > h h WEIGHTED=6 @1
// *   Decay: h > b b~ WEIGHTED=2
// *   Decay: h > w+ w- WEIGHTED=2
// *     Decay: w+ > e+ ve WEIGHTED=2
// *     Decay: w- > e- ve~ WEIGHTED=2
// Process: g g > h h WEIGHTED=6 @1
// *   Decay: h > b b~ WEIGHTED=2
// *   Decay: h > w+ w- WEIGHTED=2
// *     Decay: w+ > e+ ve WEIGHTED=2
// *     Decay: w- > mu- vm~ WEIGHTED=2
// Process: g g > h h WEIGHTED=6 @1
// *   Decay: h > b b~ WEIGHTED=2
// *   Decay: h > w+ w- WEIGHTED=2
// *     Decay: w+ > mu+ vm WEIGHTED=2
// *     Decay: w- > e- ve~ WEIGHTED=2
// Process: g g > h h WEIGHTED=6 @1
// *   Decay: h > b b~ WEIGHTED=2
// *   Decay: h > w+ w- WEIGHTED=2
// *     Decay: w+ > mu+ vm WEIGHTED=2
// *     Decay: w- > mu- vm~ WEIGHTED=2
//--------------------------------------------------------------------------

class mg5_BSM_gg_hh2bbWW_WW2lvlv
{
  public:

    // Constructor.
    mg5_BSM_gg_hh2bbWW_WW2lvlv()
    {
      for(std::size_t i = 0; i < nprocesses; ++i)
      {
        jamp2[i] = nullptr;
      }
    }

    // Destructor.
    virtual ~mg5_BSM_gg_hh2bbWW_WW2lvlv()
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
    sigmaKin(bool & firsttime);

    // Evaluate sigmaHat(sHat).
    virtual double
    sigmaHat();

    // Info on the subprocess.
    virtual std::string
    name() const
    {
      return "g g > b b~ e+ ve e- ve~ (BSM_gg_hh)";
    }

    // Set Higgs width
    virtual void
    setHiggsWidth(double higgsWidth);

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

    double kt;
    double kl;
    double c2g;
    double cg;
    double c2;

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void
    calculate_wavefunctions(const int perm[],
                            const int hel[]);

    static const int nwavefuncs = 13;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 2;
    std::complex<double> amp[namplitudes];

    double matrix_1_gg_hh_h_bbx_h_wpwm_wp_epve_wm_emvex();

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses];

    // Color flows, used when selecting color
    double * jamp2[nprocesses];

    // Pointer to the model parameters
    Parameters_BSM_gg_hh2bbWW_WW2lvlv * pars;

    // vector with external particle masses
    std::vector<double> mME;

    // vector with momenta (to be changed each event)
    std::vector < double * > p;
    // Initial particle ids
    int id1, id2;
    std::string param_card_name_;
};

#endif  // MG5_Sigma_BSM_gg_hh2bbWW_WW2lvlv_H
