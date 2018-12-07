#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TMath.h>

namespace mem
{

double roundToNdigits(double x, int n)
{
  double tmp = TMath::Power(10., n);
  if ( x != 0. ) {
    tmp /= TMath::Power(10., TMath::Floor(TMath::Log10(TMath::Abs(x))));
  }
  double x_rounded = TMath::Nint(x*tmp)/tmp;
  //std::cout << "<roundToNdigits>: x = " << x << ", x_rounded = " << x_rounded << std::endl;
  return x_rounded;
}
 
void printLorentzVector(const std::string& label, const LorentzVector& particleP4)
{
  std::cout << " " << label << ":" 
	    << " Pt = " << particleP4.pt() << ", eta = " << particleP4.eta() << ", phi = " << particleP4.phi() << ", mass = " << particleP4.mass() << std::endl;
  std::cout << "(En = " << particleP4.energy() << ", Px = " << particleP4.px() << ", Py = " << particleP4.py() << ", Pz = " << particleP4.pz() << ")" << std::endl;
}

double compCosAngle(const LorentzVector& particle1P4, const LorentzVector& particle2P4)
{
  double particle1P_times_particle2P = particle1P4.P()*particle2P4.P();
  if ( particle1P_times_particle2P > 1.e-3 ) {
    return (particle1P4.px()*particle2P4.px() + particle1P4.py()*particle2P4.py() + particle1P4.pz()*particle2P4.pz())/particle1P_times_particle2P;
  } else {
    return 0.;
  }
}

double compCosAngle(const LorentzVector& particle1P4, double particle2Theta, double particle2Phi)
{
  double cosParticle2Theta = TMath::Cos(particle2Theta);
  double sinParticle2Theta = TMath::Sin(particle2Theta);
  double cosParticle2Phi = TMath::Cos(particle2Phi);
  double sinParticle2Phi = TMath::Sin(particle2Phi);
  double particle1P = particle1P4.P();
  if ( particle1P > 1.e-3 ) {
    return (particle1P4.px()*cosParticle2Phi*sinParticle2Theta + particle1P4.py()*sinParticle2Phi*sinParticle2Theta + particle1P4.pz()*cosParticle2Theta)/particle1P;
  } else {
    return 0.;
  }
}

std::vector<double> compBJet2En_Hbb(const LorentzVector& trueBJet1P4, const LorentzVector& measuredBJet2P4)
{
  std::cout << "<compBJet2En_Hbb>:" << std::endl;
  double delta_mH = 0.5*higgsBosonMass2 - bottomQuarkMass2;
  double a = trueBJet1P4.energy();
  double b = trueBJet1P4.P()*compCosAngle(trueBJet1P4, measuredBJet2P4);
  std::cout << "cosAngle(trueBJet1P4, measuredBJet2P4) = " << compCosAngle(trueBJet1P4, measuredBJet2P4) << ": b = " << b << std::endl;
  double a2_minus_b2 = a*a - b*b;
  if ( a2_minus_b2 > 1.e-3 ) {
    double term1 = a*delta_mH;
    double term2 = delta_mH*delta_mH - (a2_minus_b2)*bottomQuarkMass2;
    if ( term2 < 0. ) term2 = 0;
    double sqrt_term2 = TMath::Sqrt(term2);
    double term3 = b*sqrt_term2; // CV: taking the absolute value is not necessary, as we consider solutions with both signs anyway
    double trueBJet2En_solution1 = (term1 + term3)/a2_minus_b2;
    double trueBJet2En_solution2 = (term1 - term3)/a2_minus_b2;
    std::cout << "trueBJet2En: solution1 = " << trueBJet2En_solution1 << ", solution2 = " << trueBJet2En_solution2 << std::endl;
    //double trueBJet2En = TMath::Max(trueBJet2En_solution1, trueBJet2En_solution2);
    //double trueBJet2En = TMath::Min(trueBJet2En_solution1, trueBJet2En_solution2);
    //return trueBJet2En;
    std::vector<double> trueBJet2En;
    trueBJet2En.push_back(trueBJet2En_solution1);
    trueBJet2En.push_back(trueBJet2En_solution2);
    //return std::pair<double, double>std::vector<double>(trueBJet2En_solution1, trueBJet2En_solution2);
    return trueBJet2En;
  } else {
    //return 0.;
    std::vector<double> trueBJet2En;
    return trueBJet2En;
  }
}

double compNuEn_Wlnu(const LorentzVector& trueChargedLeptonP4, double trueNuTheta, double trueNuPhi)
{
  double denominator = 2.*trueChargedLeptonP4.energy()*(1. - compCosAngle(trueChargedLeptonP4, trueNuTheta, trueNuPhi));
  if ( denominator > 1.e-3 ) {
    double trueNuEn = wBosonMass2/denominator;
    return trueNuEn;
  } else {
    return 0.;
  }
}

double compNuStarEn_Hww(const LorentzVector& trueEllNuEllStarP4, double trueNuStarTheta, double trueNuStarPhi)
{
  double trueEllNuEllStarMass = trueEllNuEllStarP4.mass();
  double delta_mH = 0.5*(higgsBosonMass2 - trueEllNuEllStarMass*trueEllNuEllStarMass);
  if ( !(delta_mH > 0.) ) return 0.;
  double a = trueEllNuEllStarP4.energy();
  double b = trueEllNuEllStarP4.P()*compCosAngle(trueEllNuEllStarP4, trueNuStarTheta, trueNuStarPhi);
  double a_minus_b = a - b;
  if ( a_minus_b > 1.e-3 ) {
    double trueNuStarEn = delta_mH/a_minus_b;
    return trueNuStarEn;
  } else {
    return 0.;
  }
}

double compBJetEn_top(const LorentzVector& trueEllNuP4, const LorentzVector& measuredBJetP4)
{
  double delta_mt = 0.5*(topQuarkMass2 - (wBosonMass2 + bottomQuarkMass2));
  double a = trueEllNuP4.energy();  
  double b = trueEllNuP4.P()*compCosAngle(trueEllNuP4, measuredBJetP4);
  double a2_minus_b2 = a*a - b*b;
  if ( a2_minus_b2 > 1.e-3 ) {
    double term1 = a*delta_mt;
    double term2 = delta_mt*delta_mt - (a2_minus_b2)*bottomQuarkMass2;
    if ( term2 < 0. ) term2 = 0;
    double sqrt_term2 = TMath::Sqrt(term2);
    double term3 = b*sqrt_term2; // CV: taking the absolute value is not necessary, as we consider solutions with both signs anyway
    double trueBJetEn_solution1 = (term1 + term3)/a2_minus_b2;
    double trueBJetEn_solution2 = (term1 - term3)/a2_minus_b2;
    double trueBJetEn = TMath::Max(trueBJetEn_solution1, trueBJetEn_solution2);
    return trueBJetEn;
  } else {
    return 0.;
  }
}

LorentzVector buildLorentzVector(double energy, double theta, double phi)
{
  // special version of "buildLorentzVector" function for massless particles (neutrinos)
  double cosTheta = TMath::Cos(theta);
  double sinTheta = TMath::Sin(theta);
  double cosPhi = TMath::Cos(phi);
  double sinPhi = TMath::Sin(phi);
  double px = energy*cosPhi*sinTheta;
  double py = energy*sinPhi*sinTheta;
  double pz = energy*cosTheta;
  LorentzVector p4(px, py, pz, energy);
  return p4;
}

LorentzVector buildLorentzVector(double energy, double theta, double phi, double mass)
{
  // generic version of "buildLorentzVector" function for particles with non-zero mass
  double p2 = energy*energy - mass*mass;
  if ( p2 < 0. ) p2 = 0.;
  double p = TMath::Sqrt(p2);
  double cosTheta = TMath::Cos(theta);
  double sinTheta = TMath::Sin(theta);
  double cosPhi = TMath::Cos(phi);
  double sinPhi = TMath::Sin(phi);
  double px = p*cosPhi*sinTheta;
  double py = p*sinPhi*sinTheta;
  double pz = p*cosTheta;
  LorentzVector p4(px, py, pz, energy);
  return p4;
}

double compJacobiFactor_Hbb(const LorentzVector& trueBJet1P4, const LorentzVector& trueBJet2P4)
{
  double trueBJet1Beta = trueBJet1P4.P()/trueBJet1P4.energy();
  double trueBJet2Beta = trueBJet2P4.P()/trueBJet2P4.energy();
  if ( trueBJet2Beta > 1.e-3 ) {
    double inverse_jacobiFactor = trueBJet1P4.energy()*(1. - (trueBJet1Beta/trueBJet2Beta)*compCosAngle(trueBJet1P4, trueBJet2P4));
    if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
    return 1./inverse_jacobiFactor;
  } else {
    return 0.;
  }
}
 
double compJacobiFactor_Wlnu(const LorentzVector& trueChargedLeptonP4, const LorentzVector& trueNuP4)
{
  double inverse_jacobiFactor = 2.*trueChargedLeptonP4.energy()*(1. - compCosAngle(trueChargedLeptonP4, trueNuP4));
  if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
  return 1./inverse_jacobiFactor;
}

double compJacobiFactor_Hww(const LorentzVector& trueEllNuEllStarP4, const LorentzVector& trueNuStarP4)
{
  double trueEllNuEllStarBeta = trueEllNuEllStarP4.P()/trueEllNuEllStarP4.energy();
  double inverse_jacobiFactor = trueEllNuEllStarP4.energy()*(1. - trueEllNuEllStarBeta*compCosAngle(trueEllNuEllStarP4, trueNuStarP4));
  if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
  return 1./inverse_jacobiFactor;
}

double compJacobiFactor_top(const LorentzVector& trueEllNuP4, const LorentzVector& trueBJetP4)
{
  double trueEllNuBeta = trueEllNuP4.P()/trueEllNuP4.energy();
  double trueBJetBeta = trueBJetP4.P()/trueBJetP4.energy();
  if ( trueBJetBeta > 1.e-3 ) {
    double inverse_jacobiFactor = trueEllNuP4.energy()*(1. - (trueEllNuBeta/trueBJetBeta)*compCosAngle(trueEllNuP4, trueBJetP4));
    if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
    return 1./inverse_jacobiFactor;
  } else {
    return 0.;
  }
}

}
