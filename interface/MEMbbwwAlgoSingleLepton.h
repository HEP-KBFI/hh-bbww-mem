#ifndef hhAnalysis_bbwwMEM_MEMbbwwAlgoSingleLepton_h
#define hhAnalysis_bbwwMEM_MEMbbwwAlgoSingleLepton_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoBase.h"
#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandSingleLepton_signal.h"
#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandSingleLepton_background.h"
#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorBase.h"
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"
#include "hhAnalysis/bbwwMEM/interface/MEMResult.h"
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TBenchmark.h>
#include <TMatrixD.h>
#include <TMath.h>

#include <vector> 

class MEMbbwwAlgoSingleLepton : public MEMbbwwAlgoBase
{
 public:
  MEMbbwwAlgoSingleLepton(double, const std::string&, const std::string&, const std::string&, int = 0);
  ~MEMbbwwAlgoSingleLepton();

  /// set transfer functions for b-jets and MET
  void setBJet1TF(mem::BJetTF*);
  void setBJet2TF(mem::BJetTF*);
  void setHadRecoilTF(mem::HadRecoilTF*);

  /// fix (flag=true) mass of charged lepton plus neutrino or of the jet pair originating from the decay of the "on-shell" W boson to mW,
  /// or allow the mass to vary during the integration (flag=false)
  ///
  /// Note: flag has an effect on the likelihood of the HH->bbWW signal hypothesis only (not on the likelihood of the ttbar background hypothesis)
  void applyOnshellWmassConstraint_signal(bool);

  /// set maximum number of jet pairs considered for W->jj decay (always use limited number, in order to restrict computing time)
  void setMaxNumHadWJetPairs(int maxNumHadWJetPairs)
  {
    maxNumHadWJetPairs_ = maxNumHadWJetPairs;
  }

  /// set option for sorting jet pairs considered for W->jj decay:
  ///  - decreasing compatibility with W boson mass (jet pair most compatible with W boson mass is ranked first)
  ///  - increasing dR between jets (jet pair of minimal dR(jet1,jet2) is ranked first)
  ///  - decreasing pT of jet pair (jet pair of highest pT is ranked first)
  ///  - decreasing sum(pT) of the two jets (jet pair of maximal sum pT_jet1 + pT_jet2 is ranked first)
  enum { kSortHadWJetPairsByMass, kSortHadWJetPairsByDeltaR, kSortHadWJetPairsByPt, kSortHadWJetPairsByScalarPt };  
  void setSortHadJetPairOption(int sortHadJetPairOption)
  {
    sortHadJetPairOption_ = sortHadJetPairOption;
  }
  
  /// run integration 
  void integrate(const std::vector<mem::MeasuredParticle>&, double, double, const TMatrixD&);

  /// return probabilities for signal and background hypotheses
  MEMbbwwResultSingleLepton getResult() const { return result_; }

  class MeasuredHadWJetPair
  {
   public:
    MeasuredHadWJetPair(const mem::MeasuredParticle* measuredHadWJet1, const mem::MeasuredParticle* measuredHadWJet2)
      : jet1_(measuredHadWJet1)
      , jet2_(measuredHadWJet2)
    {
      if ( jet1_ )
      {
        assert(jet1_->type() == mem::MeasuredParticle::kHadWJet);
        p4_ += jet1_->p4();
      }
      if ( jet2_ )
      {
        assert(jet2_->type() == mem::MeasuredParticle::kHadWJet);
        p4_ += jet2_->p4();
      }
    }
    ~MeasuredHadWJetPair() {}
    const mem::MeasuredParticle* jet1() const { return jet1_; }
    const mem::MeasuredParticle* jet2() const { return jet2_; }
    const mem::LorentzVector& p4() const { return p4_; }
   private:
    const mem::MeasuredParticle* jet1_;
    const mem::MeasuredParticle* jet2_;
    mem::LorentzVector p4_; 
  };

 protected:
  /// set measured momenta of charged leptons, b-jets, and (light quark) jets originating from hadronic decay of W boson
  void setMeasuredParticles(const std::vector<mem::MeasuredParticle>&);

  /// pointers to integration classes for signal and background hypotheses
  mem::MEMbbwwIntegrandSingleLepton_signal* integrand_signal_;
  bool integrand_signal_applyOnshellWmassConstraint_;	
  mem::MEMbbwwIntegrandSingleLepton_background* integrand_background_;

  /// measured momenta of charged leptons, b-jets, and (light quark) jets originating from hadronic decay of W boson
  std::vector<mem::MeasuredParticle> measuredParticles_;
  const mem::MeasuredParticle* measuredChargedLepton_;
  std::vector<const mem::MeasuredParticle*> measuredHadWJets_; 
  std::vector<MeasuredHadWJetPair> measuredHadWJetPairs_; 
  const mem::MeasuredParticle* measuredLeadingBJet_; 
  const mem::MeasuredParticle* measuredSubleadingBJet_;

  /// maximum number of jet pairs considered for W->jj decay (always use limited number, in order to restrict computing time)
  int maxNumHadWJetPairs_;  

  /// option for sorting jet pairs considered for W->jj decay 
  int sortHadJetPairOption_;

  /// result of integration (probabilities for signal and background hypotheses)
  MEMbbwwResultSingleLepton result_;
};

#endif
