#ifndef hhAnalysis_bbwwMEM_MEMResult_h
#define hhAnalysis_bbwwMEM_MEMResult_h

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

 class MEMbbwwPermutationDilepton_sig : public MEMPermutationBase
{
 public:
  MEMbbwwPermutationDilepton_sig(double prob, double probErr,
				 const MeasuredParticle& measuredBJet1, const MeasuredParticle& measuredBJet2,
				 const MeasuredParticle& measuredChargedLeptonFromOnshellW, const MeasuredParticle& measuredChargedLeptonFromOffshellW)
    : MEMPermutationBase(prob, probErr)
    , measuredBJet1_(measuredBJet1)
    , measuredBJet2_(measuredBJet2)
    , measuredChargedLeptonFromOnshellW_(measuredChargedLeptonFromOnshellW)
    , measuredChargedLeptonFromOffshellW_(measuredChargedLeptonFromOffshellW)
  {}
  ~MEMbbwwPermutationDilepton_sig()
  {}

 protected:
  MeasuredParticle measuredBJet1_;
  MeasuredParticle measuredBJet2_;
  MeasuredParticle measuredChargedLeptonFromOnshellW_;
  MeasuredParticle measuredChargedLeptonFromOffshellW_;
};

class MEMbbwwPermutationDilepton_bkg : public MEMPermutationBase
{
 public:
  MEMbbwwPermutationDilepton_bkg(double prob, double probErr, 
				 const MeasuredParticle& measuredBJetFromTop, const MeasuredParticle& measuredChargedLeptonFromTop,
				 const MeasuredParticle& measuredBJetFromAntiTop, const MeasuredParticle& measuredChargedLeptonFromAntiTop)
    : MEMPermutationBase(prob, probErr)
    , measuredBJetFromTop_(measuredBJetFromTop)
    , measuredChargedLeptonFromTop_(measuredChargedLeptonFromTop)
    , measuredBJetFromAntiTop_(measuredBJetFromAntiTop)
    , measuredChargedLeptonFromAntiTop_(measuredChargedLeptonFromAntiTop)
  {}
  ~MEMbbwwPermutationDilepton_bkg()
  {}

 protected:
  MeasuredParticle measuredBJetFromTop_;
  MeasuredParticle measuredChargedLeptonFromTop_;
  MeasuredParticle measuredBJetFromAntiTop_;
  MeasuredParticle measuredChargedLeptonFromAntiTop_;
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
    double prob_SplusB = prob_signal_ + prob_background_;    
    if ( prob_SplusB > 0. ) {
      double prob2_SplusB = mem::square(prob_SplusB);
      return TMath::Sqrt(mem::square((prob_background_/prob2_SplusB)*probErr_signal_) + mem::square((prob_signal_/prob2_SplusB)*probErr_background_));
    } else {
      return 0.;
    }
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

 protected:
  double prob_signal_;
  double probErr_signal_;
  double prob_background_;
  double probErr_background_;
  
  std::vector<T_sig> permutations_signal_;
  std::vector<T_bkg> permutations_background_;
};

typedef MEMResult<mem::MEMbbwwPermutationDilepton_sig, mem::MEMbbwwPermutationDilepton_bkg> MEMbbwwResultDilepton;

std::ostream&
operator<<(std::ostream& stream, const MEMbbwwResultDilepton& result)
{
  stream << " probability for signal hypothesis = " << result.getProb_signal() << " +/- " << result.getProbErr_signal() << std::endl;
  stream << " probability for background hypothesis = " << result.getProb_background() << " +/- " << result.getProbErr_background() << std::endl;
  stream << "--> likelihood ratio = " << result.getLikelihoodRatio() << " +/- " << result.getLikelihoodRatioErr() << std::endl;
  return stream;
}

#endif
