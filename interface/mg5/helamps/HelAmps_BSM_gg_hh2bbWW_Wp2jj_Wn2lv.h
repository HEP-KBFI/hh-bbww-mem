//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_BSM_gg_hh2bbWW_Wp2jj_Wn2lv_H
#define HelAmps_BSM_gg_hh2bbWW_Wp2jj_Wn2lv_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_BSM_gg_hh2bbWW_Wp2jj_Wn2lv
{
double Sgn(double e, double f); 

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void txxxxx(double p[4], double tmass, int nhel, int nst, std::complex<double>
    fi[18]);

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void SSS1_0(complex<double> S1[], complex<double> S2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);

void VVS1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

void FFV2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);

void FFS1_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

}  // end namespace MG5_BSM_hh2bbWW_Wp2jj_Wn2lv_hh

#endif  // HelAmps_BSM_gg_hh2bbWW_Wp2jj_Wn2lv_H
