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
/// @file pole.hpp
/// @brief Pole expansion subroutines.
/// @date 2012-10-12
#ifndef _PEXSI_POLE_HPP_
#define _PEXSI_POLE_HPP_

#include "pexsi/environment.hpp"

namespace PEXSI{

/// @brief Pole expansion for the Fermi-Dirac operator.
///
/// This is the most commonly used subroutine for the pole expansion,
/// and can be used to compute the shifts and weights for calculating
/// the density matrix, the total energy, and the Hellman-Feynman
/// force.
///
/// This routine obtains the expansion
///
/// \f[
///    f_{\beta} (z) = \frac{2}{1+e^{\beta z}} \approx
///    \mathrm{Im} \sum_{l=1}^{P} \frac{\omega^{\rho}_l}{z-z_l}
/// \f]
///
/// @note
/// The unit of `temp`,`gap`,`deltaE`,`mu` must be the same.
///
///
/// @param[out]  zshift Dimension: Npole. The shifts \f$\{z_l\}\f$.
/// @param[out]  zweight Dimension: Npole. The weights \f$\{\omega^{\rho}_l\}\f$.
/// @param[in]   Npole   Number of poles. **Must be an even number**.
/// @param[in]   temp    Temperature. Temperature equals to
/// \f$1/\beta\f$.
/// @param[in]   gap     The spectral gap, defined as
/// \f$\min_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   deltaE  The spectral range, defined as
/// \f$\max_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   mu      The chemical potential.
///
/// @return
/// - = 0: successful exit.
/// - > 0: unsuccessful.
///
int GetPoleDensity(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu);

/// @brief Pole expansion for the derivative of the Fermi-Dirac
/// operator with respect to the chemical potential mu.
///
/// This routine can be used to evaluate the derivative of the number
/// of electrons with respect to the chemical potential for the
/// Newton step for updating the chemical potential.
///
/// Note that \f$f_{\beta}\f$ does not explicitly contain \f$\mu\f$,
/// so this routine actually computes the expansion
///
/// \f[
///    -\frac{\partial f_{\beta}}{\partial z} (z) =
///    2\beta \frac{e^{\beta z}}{(1+e^{\beta z})^2}
///    \approx \mathrm{Im} \sum_{l=1}^{P}
///    \frac{\omega^{\mu}_l}{z-z_l}
/// \f]
///
/// @note
/// The unit of `temp`,`gap`,`deltaE`,`mu` must be the same.
///
///
/// @param[out]  zshift Dimension: Npole. The shifts \f$\{z_l\}\f$.
/// @param[out]  zweight Dimension: Npole. The weights \f$\{\omega^{\mu}_l\}\f$.
/// @param[in]   Npole   Number of poles. **Must be an even number**.
/// @param[in]   temp    Temperature. Temperature equals to
/// \f$1/\beta\f$.
/// @param[in]   gap     The spectral gap, defined as
/// \f$\min_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   deltaE  The spectral range, defined as
/// \f$\max_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   mu      The chemical potential.
///
/// @return
/// - = 0: successful exit.
/// - > 0: unsuccessful.
int GetPoleDensityDrvMu(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu);

/// @brief Pole expansion for the derivative of the Fermi-Dirac
/// operator with respect to the temperature T \f$(1/\beta)\f$.
///
/// This routine can be used to extrapolate the number of electrons
/// from a finite temperature calculation to a zero temperature
/// calculation, using the derivative information.  However, this
/// functionality is not used anymore in the current version of
/// %PEXSI.
///
///
/// \f[
///    \frac{\partial f_{\beta}}{\partial (1/\beta)} (z) =
///    2 \beta^2 z \frac{e^{\beta z}}{(1+e^{\beta z})^2}
///    \approx \mathrm{Im} \sum_{l=1}^{P}
///    \frac{\omega^{T}_l}{z-z_l}
/// \f]
///
/// @note
/// The unit of `temp`,`gap`,`deltaE`,`mu` must be the same.
///
///
/// @param[out]  zshift Dimension: Npole. The shifts \f$\{z_l\}\f$.
/// @param[out]  zweight Dimension: Npole. The weights \f$\{\omega^{T}_l\}\f$.
/// @param[in]   Npole   Number of poles. **Must be an even number**.
/// @param[in]   temp    Temperature. Temperature equals to
/// \f$1/\beta\f$.
/// @param[in]   gap     The spectral gap, defined as
/// \f$\min_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   deltaE  The spectral range, defined as
/// \f$\max_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   mu      The chemical potential.
///
/// @return
/// - = 0: successful exit.
/// - > 0: unsuccessful.
int GetPoleDensityDrvT(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu);

/// @brief Pole expansion for the Helmholtz free energy function.
///
/// This routine can be used to compute the (Helmholtz) free energy
/// when finite temperature effect exists. This is especially
/// important for metallic system and other small gapped systems.
/// This routine expands the free energy function
///
/// \f[
///    f^{\mathcal{F}}_{\beta}(z) = -\frac{2}{\beta} \log
///    (1 + e^{-\beta z}) \approx \mathrm{Im} \sum_{l=1}^{P}
///    \frac{\omega^{\mathcal{F}}_l}{z-z_l}
/// \f]
///
/// @note
/// The unit of `temp`,`gap`,`deltaE`,`mu` must be the same.
///
///
/// @param[out]  zshift Dimension: Npole. The shifts \f$\{z_l\}\f$.
/// @param[out]  zweight Dimension: Npole. The weights \f$\{\omega^{\mathcal{F}}_l\}\f$.
/// @param[in]   Npole   Number of poles. **Must be an even number**.
/// @param[in]   temp    Temperature. Temperature equals to
/// \f$1/\beta\f$.
/// @param[in]   gap     The spectral gap, defined as
/// \f$\min_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   deltaE  The spectral range, defined as
/// \f$\max_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   mu      The chemical potential.
///
/// @return
/// - = 0: successful exit.
/// - > 0: unsuccessful.
int GetPoleHelmholtz(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu);


/// @brief Pole expansion for the energy density function.
///
/// This routine can be used to compute the Pulay contribution of the
/// atomic force in electronic structure calculations.  This term is
/// especially important when basis set is not complete and changes
/// with atomic positions.
/// This routine expands the free energy function
///
/// \f[
///    f^{E}_{\beta}(z) = (z+\mu) f_{\beta}(z)
///    \approx \mathrm{Im} \sum_{l=1}^{P}
///    \frac{\omega^{E}_l}{z-z_l}
/// \f]
///
/// @note
/// The unit of `temp`,`gap`,`deltaE`,`mu` must be the same.
///
///
/// @param[out]  zshift Dimension: Npole. The shifts \f$\{z_l\}\f$.
/// @param[out]  zweight Dimension: Npole. The weights \f$\{\omega^{E}_l\}\f$.
/// @param[in]   Npole   Number of poles. **Must be an even number**.
/// @param[in]   temp    Temperature. Temperature equals to
/// \f$1/\beta\f$.
/// @param[in]   gap     The spectral gap, defined as
/// \f$\min_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   deltaE  The spectral range, defined as
/// \f$\max_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   mu      The chemical potential.
///
/// @return
/// - = 0: successful exit.
/// - > 0: unsuccessful.
int GetPoleForce(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu);


/// @brief Pole expansion for the Fermi-Dirac operator and update the
/// weight at mu+dmu.
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
///
/// @param[out]  zshift Dimension: Npole. The shifts \f$\{z_l\}\f$.
/// @param[out]  zweight Dimension: Npole. The weights \f$\{\omega^{\rho}_l\}\f$.
/// @param[in]   Npole   Number of poles. **Must be an even number**.
/// @param[in]   temp    Temperature. Temperature equals to
/// \f$1/\beta\f$.
/// @param[in]   gap     The spectral gap, defined as
/// \f$\min_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   deltaE  The spectral range, defined as
/// \f$\max_{\varepsilon} |\varepsilon-\mu|\f$, where
/// \f$\varepsilon\f$ is an eigenvalue.
/// @param[in]   mu      The chemical potential.
/// @param[in] dmu Update of chemical potential.
///
/// @return
/// - = 0: successful exit.
/// - > 0: unsuccessful.
///
int GetPoleDensityUpdate(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu, double dmu);

/// @brief Pole expansion for the Helmholtz free energy function.
///
/// Similar to @ref GetPoleHelmholtz but obtain the weights using the
/// update formula
int GetPoleHelmholtzUpdate(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu, double dmu);

/// @brief Pole expansion for the energy density function.
///
/// Similar to @ref GetPoleForce but obtain the weights using the
/// update formula
int GetPoleForceUpdate(Complex* zshift, Complex* zweight,
    int Npole, double temp, double gap, double deltaE,
    double mu, double dmu);

}

#endif //_PEXSI_POLE_HPP_
