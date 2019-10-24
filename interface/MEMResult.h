#ifndef hhAnalysis_bbwwMEM_MEMResult_h
#define hhAnalysis_bbwwMEM_MEMResult_h

#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TMath.h>

#include <string>
#include <vector>
#include <ostream>

namespace mem
{

class MEMPermutationBase
{
 public:
  MEMPermutationBase(double prob, double probErr)
    : prob_(prob)
    , probErr_(probErr)
  {}
  ~MEMPermutationBase()
  {}

 protected:
  double prob_;
  double probErr_;
};

class MEMbbwwPermutationDilepton : public MEMPermutationBase
{
 public:
  MEMbbwwPermutationDilepton(double prob, double probErr,
			     const MeasuredParticle* measuredBJet1, const MeasuredParticle* measuredBJet2,
			     const MeasuredParticle* measuredChargedLepton1, const MeasuredParticle* measuredChargedLepton2,
			     int chargedLeptonPermutation = kPermutationUndefined2L)
    : MEMPermutationBase(prob, probErr)
    , measuredBJet1_(measuredBJet1)
    , measuredBJet2_(measuredBJet2)
    , measuredChargedLepton1_(measuredChargedLepton1)
    , measuredChargedLepton2_(measuredChargedLepton2)
    , chargedLeptonPermutation_(chargedLeptonPermutation)
  {}
  ~MEMbbwwPermutationDilepton()
  {}

  const MeasuredParticle* getMeasuredBJet1() const
  {
    // in signal     hypothesis: "first" bottom quark (no special meaning)
    // in background hypothesis: bottom quark, originating from decay of top quark
    return measuredBJet1_;
  }
  const MeasuredParticle* getMeasuredBJet2() const 
  {
    // in signal     hypothesis: "second" bottom quark (no special meaning)
    // in background hypothesis: anti-bottom quark, originating decay of from anti-top quark
    return measuredBJet2_;
  }
  const MeasuredParticle* getMeasuredChargedLepton1() const
  {
    // in signal     hypothesis: lepton originating from decay of "on-shell" W boson
    // in background hypothesis: lepton of positive charge, originating from decay of top quark
    return measuredChargedLepton1_;
  }
  const MeasuredParticle* getMeasuredChargedLepton2() const
  {
    // in signal     hypothesis: lepton originating from decay of off-shell W boson
    // in background hypothesis: lepton of negative charge, originating decay of from anti-top quark
    return measuredChargedLepton2_;
  }

 protected:
  const MeasuredParticle* measuredBJet1_;
  const MeasuredParticle* measuredBJet2_;
  const MeasuredParticle* measuredChargedLepton1_;
  const MeasuredParticle* measuredChargedLepton2_;
  /// flag specific to signal hypothesis, 
  /// used to switch between associations of lepton+ and lepton- to on-shell and off-shell W bosons
  bool chargedLeptonPermutation_;
};

class MEMbbwwPermutationSingleLepton : public MEMPermutationBase
{
 public:
  MEMbbwwPermutationSingleLepton(double prob, double probErr,
				 const MeasuredParticle* measuredBJet1, const MeasuredParticle* measuredBJet2,
				 const MeasuredParticle* measuredChargedLepton, const MeasuredParticle* measuredHadWJet1, const MeasuredParticle* measuredHadWJet2,
				 int chargedLeptonPermutation = kPermutationUndefined1L)
    : MEMPermutationBase(prob, probErr)
    , measuredBJet1_(measuredBJet1)
    , measuredBJet2_(measuredBJet2)
    , measuredChargedLepton_(measuredChargedLepton)
    , measuredHadWJet1_(measuredHadWJet1)
    , measuredHadWJet2_(measuredHadWJet2)
    , chargedLeptonPermutation_(chargedLeptonPermutation)
  {}
  ~MEMbbwwPermutationSingleLepton()
  {}

  const MeasuredParticle* getMeasuredBJet1() const
  {
    // in signal     hypothesis: "first" bottom quark (no special meaning)
    // in background hypothesis: bottom quark, originating from decay of top quark
    return measuredBJet1_;
  }
  const MeasuredParticle* getMeasuredBJet2() const 
  {
    // in signal     hypothesis: lepton originating from decay of "on-shell" W boson
    // in background hypothesis: lepton of positive charge, originating from decay of top quark
    return measuredBJet2_;
  }
  const MeasuredParticle* getMeasuredChargedLepton() const
  {
    // in signal and background hypothesis: charged lepton (no special meaning)
    return measuredChargedLepton_;
  }
  const MeasuredParticle* getMeasuredHadWJet1() const
  {
    // in signal and background hypothesis: "first" jet originating from decay of W boson (no special meaning)
    return measuredHadWJet1_;
  }
  const MeasuredParticle* getMeasuredHadWJet2() const
  {
    // in signal and background hypothesis: "second" jet originating from decay of W boson (no special meaning)
    return measuredHadWJet2_;
  }

 protected:
  const MeasuredParticle* measuredBJet1_;
  const MeasuredParticle* measuredBJet2_;
  const MeasuredParticle* measuredChargedLepton_;
  const MeasuredParticle* measuredHadWJet1_;
  const MeasuredParticle* measuredHadWJet2_;
  /// flag specific to signal hypothesis, 
  /// used to switch between associations of charged lepton to on-shell and off-shell W bosons
  int chargedLeptonPermutation_;
};

}

template <class T_sig, class T_bkg>
class MEMResult
{
 public:
  MEMResult()
    : prob_signal_(0.)
    , probErr_signal_(0.)
    , prob_background_(0.)
    , probErr_background_(0.)
  {}
  ~MEMResult() 
  {}
  
  double getProb_signal() const
  { 
    return prob_signal_;
  }
  double getProbErr_signal() const
  { 
    return probErr_signal_;
  }
  double getProb_background() const
  {
    return prob_background_;
  }
  double getProbErr_background() const
  {
    return probErr_background_;
  }

  double getLikelihoodRatio() const
  {
    double prob_SplusB = prob_signal_ + prob_background_;    
    if ( prob_SplusB > 0. ) {
      return prob_signal_/prob_SplusB;
    } else {
      return 0.;
    }
  }
  double getLikelihoodRatioErr() const
  {
    double prob2_SplusB = mem::square(prob_signal_ + prob_background_);    
    if ( prob2_SplusB > 0. ) {
      return TMath::Sqrt(mem::square((prob_background_/prob2_SplusB)*probErr_signal_) + mem::square((prob_signal_/prob2_SplusB)*probErr_background_));
    } else {
      return 0.;
    }
  }

  double getScore() const
  {
    double memLR = this->getLikelihoodRatio();
    double memScore = TMath::Abs(TMath::Log(TMath::Max(1.e-15, TMath::Min(memLR, 1. - memLR)/0.5)));
    if ( memLR >= 0.5 ) memScore *= +1.;
    else memScore *= -1.;
    return memScore;
  }

  const std::vector<T_sig>& getPermutations_signal() const
  {
    return permutations_signal_;
  }
  const std::vector<T_bkg>& getPermutations_background() const
  {
    return permutations_background_;
  }
  
  friend class MEMbbwwAlgoDilepton;
  friend class MEMbbwwAlgoSingleLepton;

 protected:
  double prob_signal_;
  double probErr_signal_;
  double prob_background_;
  double probErr_background_;
  
  std::vector<T_sig> permutations_signal_;
  std::vector<T_bkg> permutations_background_;
};

typedef MEMResult<mem::MEMbbwwPermutationDilepton, mem::MEMbbwwPermutationDilepton> MEMbbwwResultDilepton;
typedef MEMResult<mem::MEMbbwwPermutationSingleLepton, mem::MEMbbwwPermutationSingleLepton> MEMbbwwResultSingleLepton;

std::ostream&
operator<<(std::ostream& stream, const MEMbbwwResultDilepton& result);
std::ostream&
operator<<(std::ostream& stream, const MEMbbwwResultSingleLepton& result);

#endif
