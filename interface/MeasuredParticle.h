#ifndef hhAnalysis_bbwwMEM_MeasuredTauLepton_h
#define hhAnalysis_bbwwMEM_MeasuredTauLepton_h

#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h" // LorentzVector, Vector

#include <string>

namespace mem
{
  
class MeasuredParticle
{
 public:
  /**
     \enum    MeasuredParticle::kType
     \brief   enumeration of all particle types
   */
  enum kType {
    kUndefinedType,
    kElectron,      /* < electron */ 
    kMuon,          /* < muon     */
    kBJet           /* < b-jet    */
  };

  MeasuredParticle();
  MeasuredParticle(int, double, double, double, double, int = 0);
  MeasuredParticle(const MeasuredParticle&);
  ~MeasuredParticle() {}
  
  /// return type of the particle
  int type() const { return type_; }
  std::string type_string() const;

  /// return measured pt of the particle in labframe
  double pt() const { return pt_; }
  /// return measured polar angle of the particle in labframe
  double theta() const { return theta_; }
  /// return measured pseudo-rapidity of the particle in labframe
  double eta() const { return eta_; }
  /// return measured azimuthal angle of the particle in labframe
  double phi() const { return phi_; }
  /// return mass
  double mass() const { return mass_; }    
    
  /// return measured energy of the particle in labframe
  double energy() const { return energy_; }
  /// return measured px of the particle in labframe
  double px() const { return px_; }
  /// return measured py of the particle in labframe
  double py() const { return py_; }
  /// return measured pz of the particle in labframe
  double pz() const { return pz_; }

  /// return the measured momentum of the particle in labframe
  double p() const { return p_; }

  /// return the measured charge of the particle
  int charge() const { return charge_; }    

  /// return the measured four-vector of the particle in the labframe
  const LorentzVector& p4() const { return p4_; }

  /// return the measured momentum vector of the particle in the labframe
  const Vector& p3() const { return p3_; }
    
  /// return auxiliary data-members to speed-up numerical computations
  double cosPhi_sinTheta() const { return cosPhi_sinTheta_; }
  double sinPhi_sinTheta() const { return sinPhi_sinTheta_; }
  double cosTheta() const { return cosTheta_; }

 protected:
  /// set measured momentum in all coordinates systems
  void initialize();

 private:
  /// decay type
  int type_;

  /// measured momentum in labframe (in polar coordinates)
  double pt_;
  double eta_;
  double phi_;
  double mass_;

  /// measured momentum in labframe (in cartesian coordinates)
  double energy_;
  double px_;
  double py_;
  double pz_;

  /// visible momentum in labframe (magnitude);
  double p_;  

  /// charge (electrons and muons only)
  int charge_;

  /// energy and momentum in labframe (four-vector)
  LorentzVector p4_;

  /// momentum in labframe 
  Vector p3_;

  /// auxiliary data-members to speed-up numerical computations
  double theta_;
  double cosPhi_sinTheta_;
  double sinPhi_sinTheta_;
  double cosTheta_;
};

}

#endif
