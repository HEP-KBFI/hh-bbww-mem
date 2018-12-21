//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "hhAnalysis/bbwwMEM/interface/mg5/helamps/HelAmps_BSM_gg_hh2bbWW_Wp2lv_Wn2jj.h"
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()

using namespace std; 

namespace MG5_BSM_gg_hh2bbWW_Wp2lv_Wn2jj
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], pow((pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2)), 0.5)); 
    if (pp == 0.0)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.0), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * pow(2.0 * p[0], 0.5), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
    if (pp == 0.000)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[abs(ip)]; 
      fo[3] = ip * nsf * sqm[abs(ip)]; 
      fo[4] = im * nsf * sqm[abs(im)]; 
      fo[5] = ip * sqm[abs(im)]; 
    }
    else
    {
      pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.00), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * pow(2.0 * p[0], 0.5); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}


void txxxxx(double p[4], double tmass, int nhel, int nst, complex<double>
    tc[18])
{
  complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  sqh = pow(0.5, 0.5); 
  sqs = pow(0.5/3, 0.5); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], pow(pt2 + p[3] * p[3], 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 

  ft[4][0] = complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if(nhel >= 0)
  {
    if(pp == 0)
    {
      ep[0] = complex<double> (0, 0); 
      ep[1] = complex<double> (-sqh, 0); 
      ep[2] = complex<double> (0, nst * sqh); 
      ep[3] = complex<double> (0, 0); 
    }
    else
    {
      ep[0] = complex<double> (0, 0); 
      ep[3] = complex<double> (pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = complex<double> (-sqh, 0); 
        ep[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }

  }

  // construct eps-
  if(nhel <= 0)
  {
    if(pp == 0)
    {
      em[0] = complex<double> (0, 0); 
      em[1] = complex<double> (sqh, 0); 
      em[2] = complex<double> (0, nst * sqh); 
      em[3] = complex<double> (0, 0); 
    }
    else
    {
      em[0] = complex<double> (0, 0); 
      em[3] = complex<double> (-pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = complex<double> (sqh, 0); 
        em[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }
  }

  // construct eps0
  if(fabs(nhel) <= 1)
  {
    if(pp == 0)
    {
      e0[0] = complex<double> (0, 0); 
      e0[1] = complex<double> (0, 0); 
      e0[2] = complex<double> (0, 0); 
      e0[3] = complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = complex<double> (pp/tmass, 0); 
      e0[3] = complex<double> (p[3] * emp, 0); 

      if(pt != 0)
      {
        e0[1] = complex<double> (p[1] * emp, 0); 
        e0[2] = complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = complex<double> (0, 0); 
        e0[2] = complex<double> (0, 0); 
      }
    }
  }

  if(nhel == 2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if(nhel == -2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if(tmass == 0)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = 0; 
    }
  }
  else if(tmass != 0)
  {
    if(nhel == 1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if(nhel == 0)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqs * (ep[i] * em[j] + em[i] * ep[j]
         + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if(nhel == -1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      throw_line_ext("impossible error?", TTHEXCEPTION_ERR_CODE_INVALID_HELICITY)
        << "Invalid helicity in txxxxx";
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for(j = 0; j < 4; j++ )
  {
    for(i = 0; i < 4; i++ )
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = pow(0.5, 0.5); 
  hel = double(nhel); 
  nsvahl = nsv * abs(hel); 
  pt2 = pow(p[1], 2) + pow(p[2], 2); 
  pp = min(p[0], pow(pt2 + pow(p[3], 2), 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = pow(pow(p[1], 2) + pow(p[2], 2), 0.5); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void SSS1_0(complex<double> S1[], complex<double> S2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  vertex = COUP * - cI * S3[2] * S2[2] * S1[2]; 
}


void VVS1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  complex<double> TMP3; 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP3 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP3; 
}


void FFV2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP0; 
  double P3[4]; 
  double OM3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP0 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F2[4] * F1[2] + F2[5] * F1[3] - P3[0] * OM3 * TMP0); 
  V3[3] = denom * - cI * (-F2[5] * F1[2] - F2[4] * F1[3] - P3[1] * OM3 * TMP0); 
  V3[4] = denom * - cI * (-cI * (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - P3[2]
      * OM3 * TMP0);
  V3[5] = denom * - cI * (F2[5] * F1[3] - F2[4] * F1[2] - P3[3] * OM3 * TMP0); 
}


void FFS1_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P3[4]; 
  complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * (+cI * (TMP1 + TMP2)); 
}


}  // end namespace $(namespace)s_BSM_gg_hh

