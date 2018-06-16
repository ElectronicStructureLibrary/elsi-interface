/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Author: Lin Lin

This file is part of PEXSI. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.
 */
/// @file pole.cpp
/// @brief Implementation of the pole expansion subroutines.
/// @date 2011-11-15 Original version.
/// @date 2015-11-24 Add capability for updating the chemical potential.
#include "pexsi/pole.hpp"

namespace PEXSI{

using std::sqrt;
using std::exp;
using std::log;
using std::sin;
using std::cos;

// Fermi-Dirac distribution fd(z) = 2/(1+exp(beta*z))
Complex fd(Complex z, double beta, double mu){
  Complex val, ez;
  /* LL: VERY IMPORTANT TO AVOID OVERFLOW/UNDERFLOW! */
  if( z.real() >= 0 ){
    ez = exp(-beta*z);
    val = 2.0 * ez / (1.0 + ez);
  }
  else{
    ez = exp(beta* z);
    val = 2.0 / (1.0 + ez);
  }
  return val;
}

// Derivative of the Fermi-Dirac distribution with respect to mu
// (therefore there is no minus sign)
// fdDrvMu(z) = beta * 2.0 * exp(beta*z)/(1+exp(beta*z))^2
Complex fdDrvMu(Complex z, double beta, double mu){
  Complex val, ez;
  val = fd( z, beta, mu ) * beta;
  /* LL: VERY IMPORTANT TO AVOID OVERFLOW/UNDERFLOW! */
  if( z.real() >= 0 ){
    ez = exp(-beta*z);
    val *= 1.0 / ( 1.0 + ez );
  }
  else{
    ez = exp(beta* z);
    val *= ez / ( 1.0 + ez );
  }
  return val;
}

// Derivative of the Fermi-Dirac distribution with respect to T (1/beta)
// fdDrvT(z) = beta^2 * 2.0 * exp(beta*z) * z /(1+exp(beta*z))^2
Complex fdDrvT(Complex z, double beta, double mu){
  Complex val, ez;
  val = fd( z, beta, mu ) * beta * beta * z;
  /* LL: VERY IMPORTANT TO AVOID OVERFLOW/UNDERFLOW! */
  if( z.real() >= 0 ){
    ez = exp(-beta*z);
    val *= 1.0 / ( 1.0 + ez );
  }
  else{
    ez = exp(beta* z);
    val *= ez / ( 1.0 + ez );
  }
  return val;
}



// Energy function egy(z) = (z+mu) * fd(z)
Complex egy(Complex z, double beta, double mu){
  Complex val, ez;
  /* LL: VERY IMPORTANT TO AVOID OVERFLOW/UNDERFLOW! */
  if( z.real() >= 0 ){
    ez = exp(-beta*z);
    val = (z + mu) * 2.0 * ez / (1.0 + ez);
  }
  else{
    ez = exp(beta* z);
    val = (z + mu) * 2.0 / (1.0 + ez);
  }
  return val;
}

// Helmholtz free energy function hmz(z) = -2/beta*log(1+exp(-beta*z))
Complex hmz(Complex z, double beta, double mu){
  Complex val, ez;
  /* LL: VERY IMPORTANT TO AVOID OVERFLOW/UNDERFLOW! */
  if( z.real() >= 0 ){
    ez = exp(-beta*z);
    val = -2.0/beta*log(1.0+ez);
  }
  else{
    ez = exp(beta* z);
    val = 2.0*z - 2.0/beta*log(1.0+ez);
  }
  return val;
}


/*********************************************************************
  ELLIPKKP Computes the complete elliptic integral of the first kind,
  with complement.

  Input parameters:
  L

  Output parameters
  K, Kp

  K is the value of the complete elliptic integral of the first kind,
  evaluated at M=exp(-2*pi*L), 0 < L < Inf. Kp is the result for
  complementary parameter, which is useful when M < EPS.  Even when M
  < 1e-6, the built-in ELLIPKE of MATLAB can lose digits of accuracy
  for KP.

  Recall that the elliptic modulus k is related to the parameter
  M by M = k^2.

  ELLIPKKP uses the method of the arithmetic-geometric mean described
  in 17.6 of M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
  Functions," Dover, 1965.  Same method as in ELLIPKE, only
  interchanging 1 and 1-m to find KP.

  When m=exp(-2*pi*L) is extremely small, use O(m) approximations.

  Originally written by Toby Driscoll in 1999.

  Rewritten in C++ by
  Lin Lin
  Computer Research Division, Lawrence Berkeley National Lab
  Last modified:  10-24-2012
 *********************************************************************/
void ellipkkp(double* K, double* Kp, double *pL){
  double pi, m, a0, b0, s0, a1, b1, c1, w1;
  double mm, eps = 1e-15;
  int i1;
  double L = (*pL);

  pi = atan(1.0)*4.0;
  if( L == 0 ){
    fprintf(stderr, "L == 0. STOP!\n");
    return;
  }
  if( L > 10.0 ){
    *K = pi / 2.0;
    *Kp = pi * L + log(4.0);
    return;
  }

  m = exp(-2*pi*L);

  /* Calculate K */

  a0 = 1.0;
  b0 = sqrt(1.0 - m);
  s0 = m;
  i1 = 0; mm = 1.0;
  while( mm > eps ){
    a1 = (a0+b0) * 0.5;
    b1 = sqrt(a0*b0);
    c1 = (a0-b0)/2;
    i1++;
    w1 = pow(2.0, i1) * c1 * c1;
    mm = w1;
    s0 = s0 + w1;
    a0 = a1;
    b0 = b1;
  }

  *K = pi / (2.0 * a1);

  /* Calculate Kp */
  a0 = 1.0;
  b0 = sqrt(m);
  s0 = 1.0-m;
  i1 = 0;  mm = 1.0;
  while (mm > eps ){
    a1 = (a0+b0)/2.0;
    b1 = sqrt(a0 * b0);
    c1 = (a0-b0)/2.0;
    i1++;
    w1 = pow(2.0, i1) * c1 * c1;
    mm = w1;
    s0 = s0 + w1;
    a0 = a1;
    b0 = b1;
  }
  *Kp = pi / (2.0 * a1);
  return;
}

/*********************************************************************
  ELLIPJC Calculates Jacobi elliptic functions for complex argument.

  Input parameters:
  u, L

  Output parameters
  sn, cn, dn

  [sn,cn,dn] are the values of the Jacobi elliptic functions
  evaluated at complex argument U and parameter M=exp(-2*pi*L), 0 < L
  < Inf.  Recall that M = k^2, where k is the elliptic modulus.

  The entries of U are expected to lie within the rectangle |Re U| <
  K, 0 < Im U < Kp, where K,Kp are evaluated from ELLIPKKP.

  The built-in ELLIPJ of MATLAB can't handle compelx arguments, and
  standard transformations to handle this would require ELLIPJ
  called with parameter 1-M. When M < eps (or is even close),
  this can't be done accurately.

  The algorithm is the descending Landen transformation,
  described in L. Howell's PhD thesis from MIT. Additional
  formulas from Gradshteyn & Ryzhik, 5th ed., and Abramowitz
  & Stegun.

NOTE:  When calling this function outside *flag = 0
 *flag = 1 is only for recursive use.

 Originally written by Toby Driscoll in 1999.

 Rewritten in C++ by
 Lin Lin
 Computer Research Division, Lawrence Berkeley National Lab
 Last modified:  10-24-2012
 *********************************************************************/
void ellipjc(Complex* psn, Complex* pcn,
    Complex* pdn, Complex* pu,
    double *pL, int *flag) {

  double K, Kp;
  double L = (*pL);
  double m;
  double pi = std::atan(1.0)*4.0;
  double eps = 1e-15;
  double x, kappa, mu;
  int high;
  int ione = 1;

  /* C complex numbers */
  Complex sinu, cosu;
  Complex u;
  Complex snh, cnh, dnh, v;
  Complex sn1, cn1, dn1, denom;

  u = (*pu);

  /* Check and transform u in the upper half of the rectangle */
  high = 0;
  if( *flag == 0 ){
    ellipkkp(&K, &Kp, &L);
    if( u.imag() > Kp * 0.5 ){
      high = 1;
      u = Kp - u;
    }
    m = exp(-2.0*pi*L);
  }
  else{
    m = L;
  }


  /* Case 1 : m is already small */
  if( m < 4.0 * eps ){
    sinu = sin(u);
    cosu = cos(u);

    *psn = sinu + m/4.0 * (sinu * cosu - u) * cosu;

    *pcn = cosu + m/4.0 * (-sinu * cosu + u) * sinu;

    *pdn = 1.0 +  m/4.0 * (cosu * cosu - sinu * sinu - 1.0);

  }
  /* Case 2 : m is big, call recursive formula */
  else{
    if( m > 1e-3 )
      kappa = (1.0-sqrt(1.0-m))/(1.0+sqrt(1.0-m));
    else{
      x = m/4;
      kappa = x*(1.0+x*(2.0+x*(5.0+x*(14.0+x*(42.0+x*132)))));
    }
    mu = kappa * kappa;
    v  = u / (1.0 + kappa);
    /* Call ellipjc recursively */
    ellipjc(&sn1, &cn1, &dn1, &v, &mu, &ione);
    denom = 1.0 + kappa * sn1 * sn1;

    *psn = (1.0+kappa) * sn1 / denom;
    *pcn = cn1 * dn1 / denom;
    *pdn = (1.0-kappa*sn1*sn1) / denom;
  }


  if( high ){
    snh = (*psn);
    cnh = (*pcn);
    dnh = (*pdn);

    (*psn)  = -1.0/ (sqrt(m) * snh);
    (*pcn) = Z_I * dnh / (sqrt(m)*snh);
    (*pdn) = Z_I * cnh / snh;
  }
  return;
}



/// @brief Compute the pole expansion.
///
/// Generate the poles and weights for any function f that shares the
/// same analytic structure with the Fermi-Dirac distribution with
/// chemical potential mu and inverse temperature beta.
///
/// NOTE: All units (temperature, gap, deltaE) should be the same.
///	Without specification they should all be Hartree (au).
///
///
/// Example:
///   Pseudocode (MATLAB notation) using pole expansion to reconstruct the
///   electron density.  NOTE: mu is included in zshift.

///   Rho = zeros(N, 1);
///   for l = 1 : Npoles
///     Rho = Rho + diag(imag( zweight(l) * inv(H - zshift(l)*eye(N))));
///   end
///
/// Reference:
///
///   L. Lin, J. Lu, L. Ying and W. E, Pole-based approximation of the
///   Fermi-Dirac function, Chin. Ann. Math.  30B, 729, 2009
///
/// @param[in] func  input function to be expanded by pole expansion
/// @param[out] zshift Complex shift of poles.
/// @param[out] zweight Weight of poles.
/// @param[in] Npole the number of poles to be used.
/// @param[in] temp
/// @param[in] gap Energy gap defined to be min(abs(EV-mu)).  EV is the
/// eigenvalue set of Hamiltonian
/// @param[in] deltaE Spectrum width defined to be max(EV)-min(EV). EV
/// is the eigenvalue set of Hamiltonian.
/// @param[in] mu Chemical potential.
int GetPoleFunc(Complex (*func)(Complex, double, double),
    Complex* zshift, Complex* zweight, int* Npole,
    double* temp, double* gap, double* deltaE, double* mu){
  double beta;
  double M, mshift, m2, M2, kr, L, K, Kp, coef;
  Complex t, sn, cn, dn, z, dzdt, zsq, funczsq;
  double pi = atan(1.0)*4.0;
  int i, j, Npolehalf, flag=0;
  beta = 1.0 / ((*temp));

  if( (*Npole) % 2 != 0 ){
    fprintf(stderr, "Npole has to be an even number!\n");
    return 1;
  }

  Npolehalf   = (*Npole) / 2;
  M           = *deltaE;
  mshift      = pow((pi/beta), 2.0);
  m2          = mshift + (*gap)*(*gap);
  M2          = M*M;
  kr          = (sqrt(M2/m2)-1.0)/(sqrt(M2/m2)+1.0);
  L           = -log(kr)/pi;
  ellipkkp(&K, &Kp, &L);


  for( j = 0; j < Npolehalf; j++){
    t   = (-K + (0.5 + j) / Npolehalf * 2.0 * K) + Z_I * 0.5 * Kp;
    ellipjc(&sn, &cn, &dn, &t, &L, &flag);

    z    = sqrt(m2*M2) * (1.0/kr + sn) / (1.0/kr-sn) - mshift;

    /* General formula for zshift and zweight
     *
     zshift(j)  = zsqrt(j) + mu;
     zweight(j) = 2*K*sqrt(m2*M2)/(kr*pi*Npolehalf) /
     zsqrt(j) * dzdt(j) * func(zsqrt(j));
     */

    dzdt = cn * dn / ((1.0/kr-sn) * (1.0/kr-sn));

    coef = 2.0 * K * sqrt(m2*M2) / (kr*pi*Npolehalf);


    /* The first Npolehalf poles */
    zsq     = sqrt(z);

    funczsq = (*func)(zsq, beta, *mu);

    zshift[j]  = (*mu) + zsq;
    zweight[j] = funczsq * dzdt / zsq * coef;

    /* The second Npolehalf poles */
    zsq  = -sqrt(z);

    funczsq = (*func)(zsq, beta, *mu);

    zshift[j+Npolehalf]  = (*mu) + zsq;
    zweight[j+Npolehalf] = funczsq * dzdt / zsq * coef;

  }

  return 0;
}


/// @brief Compute the update of the pole expansion.
///
/// Generate the poles and weights for any function f that shares the
/// same analytic structure with the Fermi-Dirac distribution with
/// chemical potential mu and inverse temperature beta.
///
/// The shift is given at chemical potential mu, while the weight is
/// given at mu+dmu.
///
/// @note
///
/// - This is used to evaluate the number of electrons at mu+dmu without
/// recomputing the pole expansion.  dmu should be <b>on the order of
/// \f$k_BT\f$</b> in order to be accurate.
///
/// - All units (temperature, gap, deltaE) should be the same.
///	Without specification they should all be Hartree (au).
///
/// @param[in] func  input function to be expanded by pole expansion
/// @param[out] zshift Complex shift of poles, evaluated at mu.
/// @param[out] zweight Weight of poles, evaluated at mu+dmu.
/// @param[in] Npole the number of poles to be used.
/// @param[in] temp
/// @param[in] gap Energy gap defined to be min(abs(EV-mu)).  EV is the
/// eigenvalue set of Hamiltonian
/// @param[in] deltaE Spectrum width defined to be max(EV)-min(EV). EV
/// is the eigenvalue set of Hamiltonian.
/// @param[in] mu Chemical potential.
/// @param[in] dmu Update of chemical potential.
int GetPoleUpdateFunc(Complex (*func)(Complex, double, double),
    Complex* zshift, Complex* zweight, int* Npole,
    double* temp, double* gap, double* deltaE, double* mu,
    double* dmu){
  double beta;
  double M, mshift, m2, M2, kr, L, K, Kp, coef;
  Complex t, sn, cn, dn, z, dzdt, zsq, funczsq, zsqweight;
  double pi = atan(1.0)*4.0;
  int i, j, Npolehalf, flag=0;
  beta = 1.0 / ((*temp));

  if( (*Npole) % 2 != 0 ){
    fprintf(stderr, "Npole has to be an even number!\n");
    return 1;
  }

  Npolehalf   = (*Npole) / 2;
  M           = *deltaE;
  mshift      = pow((pi/beta), 2.0);
  m2          = mshift + (*gap)*(*gap);
  M2          = M*M;
  kr          = (sqrt(M2/m2)-1.0)/(sqrt(M2/m2)+1.0);
  L           = -log(kr)/pi;
  ellipkkp(&K, &Kp, &L);


  for( j = 0; j < Npolehalf; j++){
    t   = (-K + (0.5 + j) / Npolehalf * 2.0 * K) + Z_I * 0.5 * Kp;
    ellipjc(&sn, &cn, &dn, &t, &L, &flag);

    z    = sqrt(m2*M2) * (1.0/kr + sn) / (1.0/kr-sn) - mshift;

    /* General formula for zshift and zweight
     *
     zshift(j)  = zsqrt(j) + mu;
     zweight(j) = 2*K*sqrt(m2*M2)/(kr*pi*Npolehalf) /
     zsqrt(j) * dzdt(j) * func(zsqrt(j));
     */

    dzdt = cn * dn / ((1.0/kr-sn) * (1.0/kr-sn));

    coef = 2.0 * K * sqrt(m2*M2) / (kr*pi*Npolehalf);


    /* The first Npolehalf poles */
    /* dmu reflects the update of weight */
    zsq     = sqrt(z);
    zsqweight = sqrt(z)- (*dmu);

    funczsq = (*func)(zsqweight, beta, *mu);

    zshift[j]  = (*mu) + zsq;
    zweight[j] = funczsq * dzdt / zsq * coef;

    /* The second Npolehalf poles */
    /* dmu reflects the update of weight */
    zsq  = -sqrt(z);
    zsqweight  = -sqrt(z) - (*dmu);

    funczsq = (*func)(zsqweight, beta, *mu);

    zshift[j+Npolehalf]  = (*mu) + zsq;
    zweight[j+Npolehalf] = funczsq * dzdt / zsq * coef;

  }

  return 0;
}



// Wrapper function

int GetPoleDensity(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu){
  return GetPoleFunc(&fd, zshift,  zweight,
      &Npole, &temp, &gap, &deltaE, &mu);
}

int GetPoleDensityDrvMu(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu){
  return GetPoleFunc(&fdDrvMu, zshift,  zweight,
      &Npole, &temp, &gap, &deltaE, &mu);
}

int GetPoleDensityDrvT(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu){
  return GetPoleFunc(&fdDrvT, zshift,  zweight,
      &Npole, &temp, &gap, &deltaE, &mu);
}


int GetPoleHelmholtz(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu){
  return GetPoleFunc(&hmz, zshift,  zweight,
      &Npole, &temp, &gap, &deltaE, &mu);
}

int GetPoleForce(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu){
  return GetPoleFunc(&egy, zshift,  zweight,
      &Npole, &temp, &gap, &deltaE, &mu);
}

// Wrapper for updating formula for the pole weights
int GetPoleDensityUpdate(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu, double dmu){
  return GetPoleUpdateFunc(&fd, zshift,  zweight,
      &Npole, &temp, &gap, &deltaE, &mu, &dmu);
}

int GetPoleHelmholtzUpdate(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu, double dmu){
  return GetPoleUpdateFunc(&hmz, zshift,  zweight,
      &Npole, &temp, &gap, &deltaE, &mu, &dmu);
}

int GetPoleForceUpdate(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu, double dmu){
  return GetPoleUpdateFunc(&egy, zshift,  zweight,
      &Npole, &temp, &gap, &deltaE, &mu, &dmu);
}


} //  namespace PEXSI

