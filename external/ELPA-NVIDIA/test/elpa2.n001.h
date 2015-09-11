extern "C" __device__ void elpa2_reset_shared_block_(signed char* _ps_block,signed char* _pb_size);
extern "C" __device__ void elpa2_reset_shared_block_pair_(signed char* _ps_block_1,signed char* _ps_block_2,signed char* _pb_size);
extern "C" __device__ void elpa2_warp_reduce_(signed char* _ps_block);
extern "C" __global__ void elpa2_compute_hh_dotp_kernel_(signed char* _phh,signed char* _pv_dot,int n__V_nb,int n__V_n);
extern "C" __global__ void elpa2_extract_hh_tau_kernel_(signed char* _phh,signed char* _phh_tau,int n__V_nb,int n__V_n,unsigned int n__V_is_zero);
extern "C" __global__ void elpa2_compute_hh_trafo_single_kernel_(signed char* _pq,signed char* _phh,signed char* _phh_tau,int n__V_nb,int n__V_ldq,int n__V_ncols);
extern "C" __global__ void elpa2_compute_hh_trafo_kernel_(signed char* _pq,signed char* _phh,signed char* _phh_dot,signed char* _phh_tau,int n__V_nb,int n__V_ldq,int n__V_off,int n__V_ncols);
extern "C" __global__ void elpa2_my_pack_kernel_(int n__V_n_offset,int n__V_max_idx,int n__V_stripe_width,int n__V_a_dim2,int n__V_stripe_count,signed char* _psrc,signed char* _pdst,signed char* _pdst_sd);
extern "C" __global__ void elpa2_my_unpack_kernel_(int n__V_n_offset,int n__V_max_idx,int n__V_stripe_width,int n__V_a_dim2,int n__V_stripe_count,signed char* _psrc,signed char* _pdst,signed char* _psrc_sd);
