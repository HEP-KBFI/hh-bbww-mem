#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"

#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TMath.h>

using namespace mem;

MeasuredParticle::MeasuredParticle()
  : type_(kUndefinedType),
    pt_(0.),
    eta_(0.),
    phi_(0.),
    mass_(0.),
    charge_(0)
{
  initialize();
}

MeasuredParticle::MeasuredParticle(int type, double pt, double eta, double phi, double mass, int charge) 
  : type_(type), 
    pt_(roundToNdigits(pt)),    
    eta_(roundToNdigits(eta)),
    phi_(roundToNdigits(phi)),
    mass_(roundToNdigits(mass)),
    charge_(charge)
{
  //std::cout << "<MeasuredParticle>:" << std::endl;
  //std::cout << " Pt = " << pt_ << ", eta = " << eta_ << ", phi = " << phi_ << ", mass = " << mass_ << std::endl;
  if ( type_ == kElectron ) {
    mass_ = electronMass;
  } else if ( type_ == kMuon ) {
    mass_ = muonMass;
  } else if ( type_ == kBJet ) {
    mass_ = bottomQuarkMass;
  } 
  initialize();
  //std::cout << " En = " << energy_ << ", Px = " << px_ << ", Py = " << py_ << ", Pz = " << pz_ << std::endl;
}

MeasuredParticle::MeasuredParticle(const MeasuredParticle& measuredParticle)
  : type_(measuredParticle.type()), 
    pt_(measuredParticle.pt()),
    eta_(measuredParticle.eta()),
    phi_(measuredParticle.phi()),
    mass_(measuredParticle.mass()), 
    charge_(measuredParticle.charge())     
{
  initialize();
}

void MeasuredParticle::initialize()
{
  // CV: relations between pT and p, energy taken from http://en.wikipedia.org/wiki/Pseudorapidity
  p_  = pt_*TMath::CosH(eta_);
  px_ = pt_*TMath::Cos(phi_);
  py_ = pt_*TMath::Sin(phi_);
  pz_ = pt_*TMath::SinH(eta_);
  energy_ = TMath::Sqrt(p_*p_ + mass_*mass_);
  p4_ = LorentzVector(px_, py_, pz_, energy_);
  p3_ = Vector(px_, py_, pz_);
  theta_ = p4_.theta();
  cosPhi_sinTheta_ = TMath::Cos(phi_)*TMath::Sin(theta_);
  sinPhi_sinTheta_ = TMath::Sin(phi_)*TMath::Sin(theta_);
  cosTheta_ = TMath::Cos(theta_);
}

std::string MeasuredParticle::type_string() const
{
  if      ( type_ == kElectron ) return "electron";
  else if ( type_ == kMuon     ) return "muon";
  else if ( type_ == kBJet     ) return "b-jet";
  else assert(0);
}
