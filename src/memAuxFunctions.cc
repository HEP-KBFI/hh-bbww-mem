#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TMath.h>

namespace mem
{

double 
roundToNdigits(double x, int n)
{
  double tmp = TMath::Power(10., n);
  if ( x != 0. ) 
  {
    tmp /= TMath::Power(10., TMath::Floor(TMath::Log10(TMath::Abs(x))));
  }
  double x_rounded = TMath::Nint(x*tmp)/tmp;
  return x_rounded;
}
 
void 
printLorentzVector(const std::string& label, const LorentzVector& particleP4)
{
  std::cout << label << ":" << std::endl;
  std::cout << " Pt = " << particleP4.pt() << ", eta = " << particleP4.eta() << ", phi = " << particleP4.phi() << ", mass = " << particleP4.mass() << std::endl;
  std::cout << "(En = " << particleP4.energy() << ", Px = " << particleP4.px() << ", Py = " << particleP4.py() << ", Pz = " << particleP4.pz() << ")" << std::endl;
}

void 
printLorentzVector(const std::string& label, const LorentzVector& trueParticleP4, const LorentzVector& recoParticleP4)
{
  std::cout << label << ":" << std::endl;
  std::cout << " true Pt = " << trueParticleP4.pt() << ", eta = " << trueParticleP4.eta() << ", phi = " << trueParticleP4.phi() << ", mass = " << trueParticleP4.mass() << std::endl;
  std::cout << "     (En = " << trueParticleP4.energy() << ", Px = " << trueParticleP4.px() << ", Py = " << trueParticleP4.py() << ", Pz = " << trueParticleP4.pz() << ")" << std::endl;
  std::cout << " rec. Pt = " << recoParticleP4.pt() << ", eta = " << recoParticleP4.eta() << ", phi = " << recoParticleP4.phi() << ", mass = " << recoParticleP4.mass() << std::endl;
  std::cout << "     (En = " << recoParticleP4.energy() << ", Px = " << recoParticleP4.px() << ", Py = " << recoParticleP4.py() << ", Pz = " << recoParticleP4.pz() << ")" << std::endl;
}

void 
printLorentzVector_NA(const std::string& label, const LorentzVector& trueParticleP4)
{
  std::cout << label << ":" << std::endl;
  std::cout << " true Pt = " << trueParticleP4.pt() << ", eta = " << trueParticleP4.eta() << ", phi = " << trueParticleP4.phi() << ", mass = " << trueParticleP4.mass() << std::endl;
  std::cout << "     (En = " << trueParticleP4.energy() << ", Px = " << trueParticleP4.px() << ", Py = " << trueParticleP4.py() << ", Pz = " << trueParticleP4.pz() << ")" << std::endl;
  std::cout << " rec. N/A" << std::endl;
}

double 
compCosAngle(double particle1Theta, double particle1Phi, double particle2Theta, double particle2Phi)
{
  double cosAngle = (TMath::Cos(particle1Phi)*TMath::Cos(particle2Phi) + TMath::Sin(particle1Phi)*TMath::Sin(particle2Phi))*TMath::Sin(particle1Theta)*TMath::Sin(particle2Theta) 
    + TMath::Cos(particle1Theta)*TMath::Cos(particle2Theta);
  return cosAngle;
}

std::vector<double> 
compBJetEn_Hbb(const LorentzVector& trueBJet1P4, double trueBJet2Theta, double trueBJet2Phi)
{
  double delta_mH = 0.5*higgsBosonMass2 - bottomQuarkMass2;
  double a = trueBJet1P4.energy();
  double b = trueBJet1P4.P()*compCosAngle(
    trueBJet1P4.theta(), trueBJet1P4.phi(), trueBJet2Theta, trueBJet2Phi);
  double a2_minus_b2 = a*a - b*b;
  if ( a2_minus_b2 > 1.e-3 ) 
  {
    double term1 = a*delta_mH;
    double term2 = delta_mH*delta_mH - (a2_minus_b2)*bottomQuarkMass2;
    if ( term2 < 0. ) term2 = 0;
    double sqrt_term2 = TMath::Sqrt(term2);
    double term3 = b*sqrt_term2; // CV: taking the absolute value is not necessary, as we consider solutions with both signs anyway
    double trueBJet2En_solution1 = (term1 + term3)/a2_minus_b2;
    double trueBJet2En_solution2 = (term1 - term3)/a2_minus_b2;
    std::vector<double> trueBJet2En;
    trueBJet2En.push_back(trueBJet2En_solution1);
    trueBJet2En.push_back(trueBJet2En_solution2);
    return trueBJet2En;
  } 
  else 
  {
    return std::vector<double>();
  }
}

double 
compNuEn_Wlnu(const LorentzVector& trueChargedLeptonP4, double trueNuTheta, double trueNuPhi, double q2W)
{
  double denominator = 2.*trueChargedLeptonP4.energy()*(1. - compCosAngle(
    trueChargedLeptonP4.theta(), trueChargedLeptonP4.phi(), trueNuTheta, trueNuPhi));
  if ( denominator > 1.e-3 ) 
  {
    double trueNuEn = q2W/denominator;
    return trueNuEn;
  } 
  else 
  {
    return 0.;
  }
}

double 
compNuStarEn_Hww(const LorentzVector& trueEllNuEllStarP4, double trueNuStarTheta, double trueNuStarPhi)
{
  double trueEllNuEllStarMass = trueEllNuEllStarP4.mass();
  double delta_mH = 0.5*(higgsBosonMass2 - trueEllNuEllStarMass*trueEllNuEllStarMass);
  if ( !(delta_mH > 0.) ) return 0.;
  double a = trueEllNuEllStarP4.energy();
  double b = trueEllNuEllStarP4.P()*compCosAngle(
    trueEllNuEllStarP4.theta(), trueEllNuEllStarP4.phi(), trueNuStarTheta, trueNuStarPhi);
  double a_minus_b = a - b;
  if ( a_minus_b > 1.e-3 ) 
  {
    double trueNuStarEn = delta_mH/a_minus_b;
    return trueNuStarEn;
  } 
  else 
  {
    return 0.;
  }
}

std::vector<double> 
compBJetEn_top(const LorentzVector& trueEllNuP4, double trueBJetTheta, double trueBJetPhi)
{
  double delta_mt = 0.5*(topQuarkMass2 - (wBosonMass2 + bottomQuarkMass2));
  double a = trueEllNuP4.energy();  
  double b = trueEllNuP4.P()*compCosAngle(
    trueEllNuP4.theta(), trueEllNuP4.phi(), trueBJetTheta, trueBJetPhi);
  double a2_minus_b2 = a*a - b*b;
  if ( a2_minus_b2 > 1.e-3 ) 
  {
    double term1 = a*delta_mt;
    double term2 = delta_mt*delta_mt - (a2_minus_b2)*bottomQuarkMass2;
    if ( term2 < 0. ) term2 = 0;
    double sqrt_term2 = TMath::Sqrt(term2);
    double term3 = b*sqrt_term2; // CV: taking the absolute value is not necessary, as we consider solutions with both signs anyway
    double trueBJetEn_solution1 = (term1 + term3)/a2_minus_b2;
    double trueBJetEn_solution2 = (term1 - term3)/a2_minus_b2;
    std::vector<double> trueBJetEn;
    trueBJetEn.push_back(trueBJetEn_solution1);
    trueBJetEn.push_back(trueBJetEn_solution2);
    return trueBJetEn;
  } 
  else 
  {
    return std::vector<double>();
  }
}

double 
compHadWJet2En_Wjj(const LorentzVector& trueHadWJet1P4, double trueHadWJet2Theta, double trueHadWJet2Phi, double q2W)
{
  double denominator = 2.*trueHadWJet1P4.energy()*(1. - compCosAngle(
    trueHadWJet1P4.theta(), trueHadWJet1P4.phi(), trueHadWJet2Theta, trueHadWJet2Phi));
  if ( denominator > 1.e-3 ) 
  {
    double trueHadWJet2En = q2W/denominator;
    return trueHadWJet2En;
  } 
  else 
  {
    return 0.;
  }
}

double 
compHadWJet2En_Hww(const LorentzVector& trueEllNuHadWJet1P4, double trueHadWJet2Theta, double trueHadWJet2Phi)
{
  double trueEllNuHadWJet1Mass = trueEllNuHadWJet1P4.mass();
  double delta_mH = 0.5*(higgsBosonMass2 - trueEllNuHadWJet1Mass*trueEllNuHadWJet1Mass);
  if ( !(delta_mH > 0.) ) return 0.;
  double a = trueEllNuHadWJet1P4.energy();
  double b = trueEllNuHadWJet1P4.P()*compCosAngle(
    trueEllNuHadWJet1P4.theta(), trueEllNuHadWJet1P4.phi(), trueHadWJet2Theta, trueHadWJet2Phi);
  double a_minus_b = a - b;
  if ( a_minus_b > 1.e-3 ) 
  {
    double trueHadWJet2En = delta_mH/a_minus_b;
    return trueHadWJet2En;
  } 
  else 
  {
    return 0.;
  }
}

LorentzVector 
findBJetEn_solution_Hbb(const std::vector<double>& trueBJetEn_solutions, double trueBJetTheta, double trueBJetPhi, const LorentzVector& trueOtherBJetP4, bool& trueBJetEn_foundSolution)
{
  LorentzVector trueBJetP4;
  trueBJetEn_foundSolution = false;
  double min_trueBJet_deltaMass = 1.e+3;
  for ( std::vector<double>::const_iterator trueBJetEn_solution = trueBJetEn_solutions.begin();
	trueBJetEn_solution != trueBJetEn_solutions.end(); ++trueBJetEn_solution ) 
  {
    if ( (*trueBJetEn_solution) > bottomQuarkMass ) 
    {
      LorentzVector trueBJetP4_solution = buildLorentzVector(*trueBJetEn_solution, trueBJetTheta, trueBJetPhi, bottomQuarkMass);
      double trueBJet_deltaMass = TMath::Abs((trueBJetP4_solution + trueOtherBJetP4).mass() - higgsBosonMass);
      if ( trueBJet_deltaMass < min_trueBJet_deltaMass ) 
      {
	trueBJetP4 = trueBJetP4_solution;
	trueBJetEn_foundSolution = true;
	min_trueBJet_deltaMass = trueBJet_deltaMass;
      }
    }
  }
  return trueBJetP4;
}

LorentzVector 
findBJetEn_solution_top(const std::vector<double>& trueBJetEn_solutions, double trueBJetTheta, double trueBJetPhi, const LorentzVector& trueWBosonP4, bool& trueBJetEn_foundSolution)
{
  LorentzVector trueBJetP4;
  trueBJetEn_foundSolution = false;
  double min_trueBJet_deltaMass = 1.e+3;
  for ( std::vector<double>::const_iterator trueBJetEn_solution = trueBJetEn_solutions.begin();
	trueBJetEn_solution != trueBJetEn_solutions.end(); ++trueBJetEn_solution ) 
  {
    if ( (*trueBJetEn_solution) > bottomQuarkMass ) 
    {
      LorentzVector trueBJetP4_solution = buildLorentzVector(*trueBJetEn_solution, trueBJetTheta, trueBJetPhi, bottomQuarkMass);
      double trueBJet_deltaMass = TMath::Abs((trueBJetP4_solution + trueWBosonP4).mass() - topQuarkMass);
      if ( trueBJet_deltaMass < min_trueBJet_deltaMass ) 
      {
	trueBJetP4 = trueBJetP4_solution;
	trueBJetEn_foundSolution = true;
	min_trueBJet_deltaMass = trueBJet_deltaMass;
      }
    }
  }
  return trueBJetP4;
}

LorentzVector 
buildLorentzVector(double energy, double theta, double phi)
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

LorentzVector 
buildLorentzVector(double energy, double theta, double phi, double mass)
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

double 
compJacobiFactor_Hbb(const LorentzVector& trueBJet1P4, const LorentzVector& trueBJet2P4)
{
  double trueBJet1Beta = trueBJet1P4.P()/trueBJet1P4.energy();
  double trueBJet2Beta = trueBJet2P4.P()/trueBJet2P4.energy();
  if ( trueBJet2Beta > 1.e-3 ) 
  {
    double inverse_jacobiFactor = trueBJet1P4.energy()*(1. - (trueBJet1Beta/trueBJet2Beta)*compCosAngle(
      trueBJet1P4.theta(), trueBJet1P4.phi(), trueBJet2P4.theta(), trueBJet2P4.phi()));
    if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
    return 1./inverse_jacobiFactor;
  } 
  else 
  {
    return 0.;
  }
}
 
double 
compJacobiFactor_Wlnu(const LorentzVector& trueChargedLeptonP4, const LorentzVector& trueNuP4)
{
  double inverse_jacobiFactor = 2.*trueChargedLeptonP4.energy()*(1. - compCosAngle(
    trueChargedLeptonP4.theta(), trueChargedLeptonP4.phi(), trueNuP4.theta(), trueNuP4.phi()));
  if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
  return 1./inverse_jacobiFactor;
}

double 
compJacobiFactor_Hww(const LorentzVector& trueEllNuEllStarP4, const LorentzVector& trueNuStarP4)
{
  double trueEllNuEllStarBeta = trueEllNuEllStarP4.P()/trueEllNuEllStarP4.energy();
  double inverse_jacobiFactor = trueEllNuEllStarP4.energy()*(1. - trueEllNuEllStarBeta*compCosAngle(
    trueEllNuEllStarP4.theta(), trueEllNuEllStarP4.phi(), trueNuStarP4.theta(), trueNuStarP4.phi()));
  if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
  return 1./inverse_jacobiFactor;
}

double 
compJacobiFactor_top(const LorentzVector& trueEllNuP4, const LorentzVector& trueBJetP4)
{
  double trueEllNuBeta = trueEllNuP4.P()/trueEllNuP4.energy();
  double trueBJetBeta = trueBJetP4.P()/trueBJetP4.energy();
  if ( trueBJetBeta > 1.e-3 ) 
  {
    double inverse_jacobiFactor = trueEllNuP4.energy()*(1. - (trueEllNuBeta/trueBJetBeta)*compCosAngle(
      trueEllNuP4.theta(), trueEllNuP4.phi(), trueBJetP4.theta(), trueBJetP4.phi()));
    if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
    return 1./inverse_jacobiFactor;
  } 
  else 
  {
    return 0.;
  }
}

double 
compJacobiFactor_Wjj(const LorentzVector& trueHadWJet1P4, const LorentzVector& trueHadWJet2P4)
{
  double inverse_jacobiFactor = 2.*trueHadWJet1P4.energy()*(1. - compCosAngle(
    trueHadWJet1P4.theta(), trueHadWJet1P4.phi(), trueHadWJet2P4.theta(), trueHadWJet2P4.phi()));
  if ( inverse_jacobiFactor < 1.e-3 ) inverse_jacobiFactor = 1.e-3;
  return 1./inverse_jacobiFactor;
}

}
