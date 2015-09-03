TODO List   {#pageTODO}
=========
@todo
- Add support of 64-bit integer.
- Control the output file in a better way.
- When the inertia counting keeps failing, return a nonzero info to make
  the external routine stop.
- Add a routine to estimate the spectral radius.
- Remove the release of TODO.md and test/ folders since not useful for
  users.

- Resolve the conflict in data format between symmetric and unsymmetric
  matrices.
- Add an sample matrix for the unsymmetric case.
- (Mathias) Update the functionality documentation (templated class DONE, conversion to SuperLU DONE, distinction between the symmetric and asymmetric PMatrix class).
- (Lin) Add the FORTRAN interface for the support of unsymmetric matrices.
- (future) Add expert user interface for PEXSI for k-point sampling and for spin
  related calculations.
- x Release of some utility subroutines?
- x Update software dependence in the config files
- (Mathias) More "standard" config/make/make install process --> cmake option?
- (Lin) MultiSolve routine does not work for complex arithmetic.
- (Mathias) Equibration option should be turn on for symm / asymm matrix? asymm for sure. CAN'T BE USED
- ROWPERM option for asymm matrices -- example for a ill-conditioned matrix: test the efficiency v.s. accuracy
