 /*! \brief C interface for the implementation of the elpa_allocate method
 *
 *  \param  none
 *  \result elpa_t handle
 */
 /*! \brief C interface for the implementation of the elpa_deallocate method
 *
 *  \param  elpa_t  handle of ELPA object to be deallocated
 *  \param  int*    error code
 *  \result void
 */
 /*! \brief C interface for the implementation of the elpa_load_settings method
 *
 *  \param elpa_t handle
 *  \param  char* filename
 */
 void elpa_load_settings(elpa_t handle, const char *filename, int *error);
 /*! \brief C interface for the implementation of the elpa_print_settings method
 *
 *  \param elpa_t handle
 *  \param  char* filename
 */
 void elpa_print_settings(elpa_t handle, int *error);
 /*! \brief C interface for the implementation of the elpa_store_settings method
 *
 *  \param elpa_t handle
 *  \param  char* filename
 */
 void elpa_store_settings(elpa_t handle, const char *filename, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_deallocate method
 *
 *  \param  elpa_autotune_impl_t  handle of ELPA autotune object to be deallocated
 *  \result void
 */
 /*! \brief C interface for the implementation of the elpa_setup method
 *
 *  \param  elpa_t  handle of the ELPA object which describes the problem to
 *                  be set up
 *  \result int     error code, which can be queried with elpa_strerr
 */
 int elpa_setup(elpa_t handle);
 /*! \brief C interface for the implementation of the elpa_set_integer method
 *  This method is available to the user as C generic elpa_set method
 *
 *  \param  handle  handle of the ELPA object for which a key/value pair should be set
 *  \param  name    the name of the key
 *  \param  value   the value to be set for the key
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
 void elpa_set_integer(elpa_t handle, const char *name, int value, int *error);
 /*! \brief C interface for the implementation of the elpa_get_integer method
 *  This method is available to the user as C generic elpa_get method
 *
 *  \param  handle  handle of the ELPA object for which a key/value pair should be queried
 *  \param  name    the name of the key
 *  \param  value   the value to be obtain for the key
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
 void elpa_get_integer(elpa_t handle, const char *name, int *value, int *error);
 /*! \brief C interface for the implementation of the elpa_set_float method
 *  This method is available to the user as C generic elpa_set method
 *
 *  \param  handle  handle of the ELPA object for which a key/value pair should be set
 *  \param  name    the name of the key
 *  \param  value   the value to be set for the key
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
 void elpa_set_float(elpa_t handle, const char *name, float value, int *error);
 /*! \brief C interface for the implementation of the elpa_get_float method
 *  This method is available to the user as C generic elpa_get method
 *
 *  \param  handle  handle of the ELPA object for which a key/value pair should be queried
 *  \param  name    the name of the key
 *  \param  value   the value to be obtain for the key
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
 void elpa_get_float(elpa_t handle, const char *name, float *value, int *error);
 /*! \brief C interface for the implementation of the elpa_set_double method
 *  This method is available to the user as C generic elpa_set method
 *
 *  \param  handle  handle of the ELPA object for which a key/value pair should be set
 *  \param  name    the name of the key
 *  \param  value   the value to be set for the key
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
 void elpa_set_double(elpa_t handle, const char *name, double value, int *error);
 /*! \brief C interface for the implementation of the elpa_get_double method
 *  This method is available to the user as C generic elpa_get method
 *
 *  \param  handle  handle of the ELPA object for which a key/value pair should be queried
 *  \param  name    the name of the key
 *  \param  value   the value to be obtain for the key
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
 void elpa_get_double(elpa_t handle, const char *name, double *value, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_set_api_version method
 *
 *  \param  elpa_t           handle: of the ELPA object which should be tuned
 *  \param  int              api_version: the version used for autotuning
 */
 void elpa_autotune_set_api_version(elpa_t handle, int api_version, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_setup method
 *
 *  \param  elpa_t           handle: of the ELPA object which should be tuned
 *  \param  int              level:  "thoroughness" of autotuning
 *  \param  int              domain: real/complex autotuning
 *  \result elpa_autotune_t  handle:  on the autotune object
 */
 elpa_autotune_t elpa_autotune_setup(elpa_t handle, int level, int domain, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_step method
 *
 *  \param  elpa_t           handle: of the ELPA object which should be tuned
 *  \param  elpa_autotune_t  autotune_handle: the autotuning object
 *  \param  error            int *error code
 *  \result int              unfinished:  describes whether autotuning finished (0) or not (1)
 */
 int elpa_autotune_step(elpa_t handle, elpa_autotune_t autotune_handle, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_print_state method
 *
 *  \param  elpa_t           handle: of the ELPA object which should be tuned
 *  \param  elpa_autotune_t  autotune_handle: the autotuning object
 *  \param  error            int *
 *  \result none
 */
 void elpa_autotune_print_state(elpa_t handle, elpa_autotune_t autotune_handle, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_save_state method
 *
 *  \param  elpa_t           handle: of the ELPA object which should be tuned
 *  \param  elpa_autotune_t  autotune_handle: the autotuning object
 *  \param  error            int *
 *  \result none
 */
 void elpa_autotune_save_state(elpa_t handle, elpa_autotune_t autotune_handle, const char *filename, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_load_state method
 *
 *  \param  elpa_t           handle: of the ELPA object which should be tuned
 *  \param  elpa_autotune_t  autotune_handle: the autotuning object
 *  \param  error            int *
 *  \result none
 */
 void elpa_autotune_load_state(elpa_t handle, elpa_autotune_t autotune_handle, const char *filename, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_set_best method
 *
 *  \param  elpa_t           handle: of the ELPA object which should be tuned
 *  \param  elpa_autotune_t  autotune_handle: the autotuning object
 *  \param  error            int *
 *  \result none
 */
 void elpa_autotune_set_best(elpa_t handle, elpa_autotune_t autotune_handle, int *error);
 /*! \brief C interface for the implementation of the elpa_autotune_print_best method
 *
 *  \param  elpa_t           handle: of the ELPA object which should be tuned
 *  \param  elpa_autotune_t  autotune_handle: the autotuning object
 *  \param  error            int *
 *  \result none
 */
 void elpa_autotune_print_best(elpa_t handle, elpa_autotune_t autotune_handle, int *error);
 void elpa_eigenvectors_all_host_arrays_d(elpa_t handle, double *a, double *ev, double *q, int *error);
 void elpa_eigenvectors_all_host_arrays_f(elpa_t handle, float *a, float *ev, float *q, int *error);
 void elpa_eigenvectors_all_host_arrays_dc(elpa_t handle, double complex *a, double *ev, double complex *q, int *error);
 void elpa_eigenvectors_all_host_arrays_fc(elpa_t handle, float complex *a, float *ev, float complex *q, int *error);
 void elpa_eigenvectors_device_pointer_d(elpa_t handle, double *a, double *ev, double *q, int *error);
 void elpa_eigenvectors_device_pointer_f(elpa_t handle, float *a, float *ev, float *q, int *error);
 void elpa_eigenvectors_device_pointer_dc(elpa_t handle, double complex *a, double *ev, double complex *q, int *error);
 void elpa_eigenvectors_device_pointer_fc(elpa_t handle, float complex *a, float *ev, float complex *q, int *error);
 void elpa_skew_eigenvectors_all_host_arrays_d(elpa_t handle, double *a, double *ev, double *q, int *error);
 void elpa_skew_eigenvectors_all_host_arrays_f(elpa_t handle, float *a, float *ev, float *q, int *error);
 void elpa_skew_eigenvectors_device_pointer_d(elpa_t handle, double *a, double *ev, double *q, int *error);
 void elpa_skew_eigenvectors_device_pointer_f(elpa_t handle, float *a, float *ev, float *q, int *error);
 void elpa_eigenvalues_all_host_arrays_d(elpa_t handle, double *a, double *ev, int *error);
 void elpa_eigenvalues_all_host_arrays_f(elpa_t handle, float *a, float *ev, int *error);
 void elpa_eigenvalues_all_host_arrays_dc(elpa_t handle, double complex *a, double *ev, int *error);
 void elpa_eigenvalues_all_host_arrays_fc(elpa_t handle, float complex *a, float *ev, int *error);
 void elpa_eigenvalues_device_pointer_d(elpa_t handle, double *a, double *ev, int *error);
 void elpa_eigenvalues_device_pointer_f(elpa_t handle, float *a, float *ev, int *error);
 void elpa_eigenvalues_device_pointer_dc(elpa_t handle, double complex *a, double *ev, int *error);
 void elpa_eigenvalues_device_pointer_fc(elpa_t handle, float complex *a, float *ev, int *error);
 void elpa_skew_eigenvalues_all_host_arrays_d(elpa_t handle, double *a, double *ev, int *error);
 void elpa_skew_eigenvalues_all_host_arrays_f(elpa_t handle, float *a, float *ev, int *error);
 void elpa_skew_eigenvalues_device_pointer_d(elpa_t handle, double *a, double *ev, int *error);
 void elpa_skew_eigenvalues_device_pointer_f(elpa_t handle, float *a, float *ev, int *error);
 void elpa_generalized_eigenvectors_d(elpa_t handle, double *a, double *b, double *ev, double *q,
 int is_already_decomposed, int *error);
 void elpa_generalized_eigenvectors_f(elpa_t handle, float *a, float *b, float *ev, float *q,
 int is_already_decomposed, int *error);
 void elpa_generalized_eigenvectors_dc(elpa_t handle, double complex *a, double complex *b, double *ev, double complex *q,
 int is_already_decomposed, int *error);
 void elpa_generalized_eigenvectors_fc(elpa_t handle, float complex *a, float complex *b, float *ev, float complex *q,
 int is_already_decomposed, int *error);
 void elpa_generalized_eigenvalues_d(elpa_t handle, double *a, double *b, double *ev,
 int is_already_decomposed, int *error);
 void elpa_generalized_eigenvalues_f(elpa_t handle, float *a, float *b, float *ev,
 int is_already_decomposed, int *error);
 void elpa_generalized_eigenvalues_dc(elpa_t handle, double complex *a, double complex *b, double *ev,
 int is_already_decomposed, int *error);
 void elpa_generalized_eigenvalues_fc(elpa_t handle, float complex *a, float complex *b, float *ev,
 int is_already_decomposed, int *error);
 void elpa_hermitian_multiply_d(elpa_t handle, char uplo_a, char uplo_c, int ncb, double *a, double *b, int nrows_b, int ncols_b, double *c, int nrows_c, int ncols_c, int *error);
 void elpa_hermitian_multiply_df(elpa_t handle, char uplo_a, char uplo_c, int ncb, float *a, float *b, int nrows_b, int ncols_b, float *c, int nrows_c, int ncols_c, int *error);
 void elpa_hermitian_multiply_dc(elpa_t handle, char uplo_a, char uplo_c, int ncb, double complex *a, double complex *b, int nrows_b, int ncols_b, double complex *c, int nrows_c, int ncols_c, int *error);
 void elpa_hermitian_multiply_fc(elpa_t handle, char uplo_a, char uplo_c, int ncb, float complex *a, float complex *b, int nrows_b, int ncols_b, float complex *c, int nrows_c, int ncols_c, int *error);
 void elpa_cholesky_d(elpa_t handle, double *a, int *error);
 void elpa_cholesky_f(elpa_t handle, float *a, int *error);
 void elpa_cholesky_dc(elpa_t handle, double complex *a, int *error);
 void elpa_cholesky_fc(elpa_t handle, float complex *a, int *error);
 void elpa_invert_trm_d(elpa_t handle, double *a, int *error);
 void elpa_invert_trm_f(elpa_t handle, float *a, int *error);
 void elpa_invert_trm_dc(elpa_t handle, double complex *a, int *error);
 void elpa_invert_trm_fc(elpa_t handle, float complex *a, int *error);
 int elpa_init(int api_version);
 #if OPTIONAL_C_ERROR_ARGUMENT == 1
 #define elpa_allocate(...) CONC(elpa_allocate, NARGS(__VA_ARGS__))(__VA_ARGS__)
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT == 1
 elpa_t elpa_allocate2(int *error);
 elpa_t elpa_allocate1();
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT == 1
 #define NARGS(...) NARGS_(__VA_ARGS__, 5, 4, 3, 2, 1, 0)
 #define NARGS_(_5, _4, _3, _2, _1, N, ...) N
 #define CONC(A, B) CONC_(A, B)
 #define CONC_(A, B) A##B
 #define elpa_deallocate(...) CONC(elpa_deallocate, NARGS(__VA_ARGS__))(__VA_ARGS__)
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT == 1
 void elpa_deallocate2(elpa_t handle, int *error);
 void elpa_deallocate1(elpa_t handle);
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT == 1
 #define elpa_autotune_deallocate(...) CONC(elpa_autotune_deallocate, NARGS(__VA_ARGS__))(__VA_ARGS__)
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT == 1
 void elpa_autotune_deallocate2(elpa_autotune_t handle, int *error);
 void elpa_autotune_deallocate1(elpa_autotune_t handle);
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT == 1
 #define elpa_uninit(...) CONC(elpa_uninit, NARGS(__VA_ARGS__))(__VA_ARGS__)
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT == 1
 void elpa_uninit1(int *error);
 void elpa_uninit0();
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT != 1
 elpa_t elpa_allocate(int *error);
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT != 1
 void elpa_deallocate(elpa_t handle, int *error);
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT != 1
 void elpa_autotune_deallocate(elpa_autotune_t handle, int *error);
 #endif
 #if OPTIONAL_C_ERROR_ARGUMENT != 1
 void elpa_uninit(int *error);
 #endif
