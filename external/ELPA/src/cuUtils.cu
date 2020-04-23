//    This file is part of ELPA.
//
//    The ELPA library was originally created by the ELPA consortium,
//    consisting of the following organizations:
//
//    - Max Planck Computing and Data Facility (MPCDF), formerly known as
//      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
//    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
//      Informatik,
//    - Technische Universität München, Lehrstuhl für Informatik mit
//      Schwerpunkt Wissenschaftliches Rechnen ,
//    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
//    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
//      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
//      and
//    - IBM Deutschland GmbH
//
//    This particular source code file contains additions, changes and
//    enhancements authored by Intel Corporation which is not part of
//    the ELPA consortium.
//
//    More information can be found here:
//    http://elpa.mpcdf.mpg.de/
//
//    ELPA is free software: you can redistribute it and/or modify
//    it under the terms of the version 3 of the license of the
//    GNU Lesser General Public License as published by the Free
//    Software Foundation.
//
//    ELPA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
//
//    ELPA reflects a substantial effort on the part of the original
//    ELPA consortium, and we ask you to respect the spirit of the
//    license that we chose: i.e., please contribute any changes you
//    may have back to the original ELPA library distribution, and keep
//    any derivatives of ELPA under the same license that we chose for
//    the original distribution, the GNU Lesser General Public License.
//
//
// --------------------------------------------------------------------------------------------------
//
// This file was originally written by NVIDIA
// and re-written by A. Marek, MPCDF

#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include <cuComplex.h>

#define MAX_BLOCK_SIZE 1024

__global__ void my_pack_c_kernel_real_double(const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int l_nev, double *src, double *dst, int i_off)
{
    int b_id = blockIdx.y;
    int t_id = threadIdx.x + i_off * blockDim.x;
    int dst_ind = b_id * stripe_width + t_id;

    if (dst_ind < max_idx)
    {

        *(dst + dst_ind + (l_nev * blockIdx.x)) = *(src + t_id + (stripe_width * (n_offset + blockIdx.x)) + (b_id * stripe_width * a_dim2));

    }
}

__global__ void my_unpack_c_kernel_real_double(const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int l_nev, double *src, double *dst, int i_off)
{
    int b_id = blockIdx.y;
    int t_id = threadIdx.x + i_off * blockDim.x;
    int src_ind = b_id * stripe_width + t_id;

    if (src_ind < max_idx)
    {

        *(dst + (t_id + ((n_offset + blockIdx.x) * stripe_width) + (b_id * stripe_width * a_dim2))) = *(src + src_ind + (blockIdx.x) * l_nev);

    }
}

__global__ void extract_hh_tau_c_kernel_real_double(double *hh, double *hh_tau, const int nbw, const int n, int val)
{
    int h_idx = (blockIdx.x) * blockDim.x + threadIdx.x;

    if (h_idx < n)
    {

        *(hh_tau + h_idx) = *(hh + (h_idx * nbw));

        if (val == 0)
        {
            *(hh + (h_idx * nbw)) = 1.0;
        }
        else
        {
            *(hh + (h_idx * nbw)) = 0.0;
        }
     }
}

extern "C" void launch_my_pack_c_kernel_real_double(const int row_count, const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int stripe_count, const int l_nev, double *a_dev, double *row_group_dev)
{
    cudaError_t err;
    dim3 grid_size = dim3(row_count, stripe_count, 1);
    int blocksize = stripe_width > MAX_BLOCK_SIZE ? MAX_BLOCK_SIZE : stripe_width;

    for (int i_off = 0; i_off < stripe_width / blocksize; i_off++)
    {

        my_pack_c_kernel_real_double<<<grid_size, blocksize>>>(n_offset, max_idx, stripe_width, a_dim2, l_nev, a_dev, row_group_dev, i_off);
    }

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n my pack_kernel failed %s \n", cudaGetErrorString(err));
    }
}

extern "C" void launch_extract_hh_tau_c_kernel_real_double(double *bcast_buffer_dev, double *hh_tau_dev, const int nbw, const int n, const int is_zero)
{
    cudaError_t err;
    int grid_size = 1 + (n - 1) / MAX_BLOCK_SIZE;

    extract_hh_tau_c_kernel_real_double<<<grid_size, MAX_BLOCK_SIZE>>>(bcast_buffer_dev, hh_tau_dev, nbw, n, is_zero);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n extract _kernel failed %s \n", cudaGetErrorString(err));
    }
}

extern "C" void launch_my_unpack_c_kernel_real_double(const int row_count, const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int stripe_count, const int l_nev, double *row_group_dev, double *a_dev)
{
    cudaError_t err;
    dim3 grid_size = dim3(row_count, stripe_count, 1);
    int blocksize = stripe_width > MAX_BLOCK_SIZE ? MAX_BLOCK_SIZE : stripe_width;

    for (int i_off = 0; i_off < stripe_width / blocksize; i_off++)
    {

        my_unpack_c_kernel_real_double<<<grid_size, blocksize>>>(n_offset, max_idx, stripe_width, a_dim2, l_nev, row_group_dev, a_dev, i_off);
    }

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n my_unpack_c_kernel failed %s \n", cudaGetErrorString(err));
    }
}

extern "C" int cuda_MemcpyDeviceToDevice(int val)
{
    val = cudaMemcpyDeviceToDevice;
    return val;
}

__global__ void my_pack_c_kernel_real_single(const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int l_nev, float *src, float *dst, int i_off)
{
    int b_id = blockIdx.y;
    int t_id = threadIdx.x + i_off * blockDim.x;
    int dst_ind = b_id * stripe_width + t_id;

    if (dst_ind < max_idx)
    {

        *(dst + dst_ind + (l_nev * blockIdx.x)) = *(src + t_id + (stripe_width * (n_offset + blockIdx.x)) + (b_id * stripe_width * a_dim2));

    }
}

__global__ void my_unpack_c_kernel_real_single(const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int l_nev, float *src, float *dst, int i_off)
{
    int b_id = blockIdx.y;
    int t_id = threadIdx.x + i_off * blockDim.x;
    int src_ind = b_id * stripe_width + t_id;

    if (src_ind < max_idx)
    {

        *(dst + (t_id + ((n_offset + blockIdx.x) * stripe_width) + (b_id * stripe_width * a_dim2))) = *(src + src_ind + (blockIdx.x) * l_nev);

    }
}

__global__ void extract_hh_tau_c_kernel_real_single(float *hh, float *hh_tau, const int nbw, const int n, int val)
{
    int h_idx = (blockIdx.x) * blockDim.x + threadIdx.x;

    if (h_idx < n)
    {

        *(hh_tau + h_idx) = *(hh + (h_idx * nbw));

        if (val == 0)
        {
            *(hh + (h_idx * nbw)) = 1.0;
        }
        else
        {
            *(hh + (h_idx * nbw)) = 0.0;
        }
     }
}

extern "C" void launch_my_pack_c_kernel_real_single(const int row_count, const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int stripe_count, const int l_nev, float *a_dev, float *row_group_dev)
{
    cudaError_t err;
    dim3 grid_size = dim3(row_count, stripe_count, 1);
    int blocksize = stripe_width > MAX_BLOCK_SIZE ? MAX_BLOCK_SIZE : stripe_width;

    for (int i_off = 0; i_off < stripe_width / blocksize; i_off++)
    {

        my_pack_c_kernel_real_single<<<grid_size, blocksize>>>(n_offset, max_idx, stripe_width, a_dim2, l_nev, a_dev, row_group_dev, i_off);
    }

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n my pack_kernel failed %s \n", cudaGetErrorString(err));
    }
}

extern "C" void launch_extract_hh_tau_c_kernel_real_single(float *bcast_buffer_dev, float *hh_tau_dev, const int nbw, const int n, const int is_zero)
{
    cudaError_t err;
    int grid_size = 1 + (n - 1) / MAX_BLOCK_SIZE;

    extract_hh_tau_c_kernel_real_single<<<grid_size, MAX_BLOCK_SIZE>>>(bcast_buffer_dev, hh_tau_dev, nbw, n, is_zero);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n extract _kernel failed %s \n", cudaGetErrorString(err));
    }
}

extern "C" void launch_my_unpack_c_kernel_real_single(const int row_count, const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int stripe_count, const int l_nev, float *row_group_dev, float *a_dev)
{
    cudaError_t err;
    dim3 grid_size = dim3(row_count, stripe_count, 1);
    int blocksize = stripe_width > MAX_BLOCK_SIZE ? MAX_BLOCK_SIZE : stripe_width;

    for (int i_off = 0; i_off < stripe_width / blocksize; i_off++)
    {

        my_unpack_c_kernel_real_single<<<grid_size, blocksize>>>(n_offset, max_idx, stripe_width, a_dim2, l_nev, row_group_dev, a_dev, i_off);
    }

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n my_unpack_c_kernel failed %s \n", cudaGetErrorString(err));
    }
}
__global__ void my_pack_c_kernel_complex_double(const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int l_nev, cuDoubleComplex *src, cuDoubleComplex *dst, int i_off)

{
    int b_id = blockIdx.y;
    int t_id = threadIdx.x + i_off * blockDim.x;
    int dst_ind = b_id * stripe_width + t_id;

    if (dst_ind < max_idx)
    {

        dst[dst_ind + (l_nev * blockIdx.x)].x = src[t_id + (stripe_width * (n_offset + blockIdx.x)) + (b_id * stripe_width * a_dim2)].x;
        dst[dst_ind + (l_nev * blockIdx.x)].y = src[t_id + (stripe_width * (n_offset + blockIdx.x)) + (b_id * stripe_width * a_dim2)].y;

    }
}
__global__ void my_unpack_c_kernel_complex_double(const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int l_nev, cuDoubleComplex *src, cuDoubleComplex *dst, int i_off)

{
    int b_id = blockIdx.y;
    int t_id = threadIdx.x + i_off * blockDim.x;
    int src_ind = b_id * stripe_width + t_id;

    if (src_ind < max_idx)
    {

        dst[t_id + ((n_offset + blockIdx.x) * stripe_width) + (b_id * stripe_width * a_dim2)].x = src[src_ind + (blockIdx.x) * l_nev].x;
        dst[t_id + ((n_offset + blockIdx.x) * stripe_width) + (b_id * stripe_width * a_dim2)].y = src[src_ind + (blockIdx.x) * l_nev].y;

    }
}
__global__ void extract_hh_tau_c_kernel_complex_double(cuDoubleComplex *hh, cuDoubleComplex *hh_tau, const int nbw, const int n, int val)

{
    int h_idx = (blockIdx.x) * blockDim.x + threadIdx.x;

    if (h_idx < n)
    {

        hh_tau[h_idx] = hh[h_idx * nbw];
        if (val == 0)
        {
            hh[(h_idx * nbw)].x = 1.0;
            hh[h_idx * nbw].y= 0.0;
        }
        else
        {
            hh[(h_idx * nbw)].x = 0.0;
            hh[h_idx * nbw].y =0.0;
        }

     }
}
extern "C" void launch_my_pack_c_kernel_complex_double(const int row_count, const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int stripe_count, const int l_nev, cuDoubleComplex *a_dev, cuDoubleComplex *row_group_dev)

{
    cudaError_t err;
    dim3 grid_size = dim3(row_count, stripe_count, 1);
    int blocksize = stripe_width > MAX_BLOCK_SIZE ? MAX_BLOCK_SIZE : stripe_width;

    for (int i_off = 0; i_off < stripe_width / blocksize; i_off++)
    {
        my_pack_c_kernel_complex_double<<<grid_size, blocksize>>>(n_offset, max_idx, stripe_width, a_dim2, l_nev, a_dev, row_group_dev, i_off);

    }

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n my pack_kernel failed %s \n", cudaGetErrorString(err));
    }
}
extern "C" void launch_extract_hh_tau_c_kernel_complex_double(cuDoubleComplex *bcast_buffer_dev, cuDoubleComplex *hh_tau_dev, const int nbw, const int n, const int is_zero)

{
    cudaError_t err;
    int grid_size = 1 + (n - 1) / MAX_BLOCK_SIZE;
    extract_hh_tau_c_kernel_complex_double<<<grid_size, MAX_BLOCK_SIZE>>>(bcast_buffer_dev, hh_tau_dev, nbw, n, is_zero);

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n extract _kernel failed %s \n", cudaGetErrorString(err));
    }
}
extern "C" void launch_my_unpack_c_kernel_complex_double(const int row_count, const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int stripe_count, const int l_nev, cuDoubleComplex *row_group_dev, cuDoubleComplex *a_dev)

{
    cudaError_t err;
    dim3 grid_size = dim3(row_count, stripe_count, 1);
    int blocksize = stripe_width > MAX_BLOCK_SIZE ? MAX_BLOCK_SIZE : stripe_width;

    for (int i_off = 0; i_off < stripe_width / blocksize; i_off++)
    {
        my_unpack_c_kernel_complex_double<<<grid_size, blocksize>>>(n_offset, max_idx, stripe_width, a_dim2, l_nev, row_group_dev, a_dev, i_off);

    }

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n my_unpack_c_kernel failed %s \n", cudaGetErrorString(err));
    }
}
__global__ void my_pack_c_kernel_complex_single(const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int l_nev, cuFloatComplex *src, cuFloatComplex *dst, int i_off)

{
    int b_id = blockIdx.y;
    int t_id = threadIdx.x + i_off * blockDim.x;
    int dst_ind = b_id * stripe_width + t_id;

    if (dst_ind < max_idx)
    {

        dst[dst_ind + (l_nev * blockIdx.x)].x = src[t_id + (stripe_width * (n_offset + blockIdx.x)) + (b_id * stripe_width * a_dim2)].x;
        dst[dst_ind + (l_nev * blockIdx.x)].y = src[t_id + (stripe_width * (n_offset + blockIdx.x)) + (b_id * stripe_width * a_dim2)].y;

    }
}
__global__ void my_unpack_c_kernel_complex_single(const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int l_nev, cuFloatComplex *src, cuFloatComplex *dst, int i_off)

{
    int b_id = blockIdx.y;
    int t_id = threadIdx.x + i_off * blockDim.x;
    int src_ind = b_id * stripe_width + t_id;

    if (src_ind < max_idx)
    {

        dst[t_id + ((n_offset + blockIdx.x) * stripe_width) + (b_id * stripe_width * a_dim2)].x = src[src_ind + (blockIdx.x) * l_nev].x;
        dst[t_id + ((n_offset + blockIdx.x) * stripe_width) + (b_id * stripe_width * a_dim2)].y = src[src_ind + (blockIdx.x) * l_nev].y;

    }
}
__global__ void extract_hh_tau_c_kernel_complex_single(cuFloatComplex *hh, cuFloatComplex *hh_tau, const int nbw, const int n, int val)

{
    int h_idx = (blockIdx.x) * blockDim.x + threadIdx.x;

    if (h_idx < n)
    {

        hh_tau[h_idx] = hh[h_idx * nbw];
        if (val == 0)
        {
            hh[(h_idx * nbw)].x = 1.0;
            hh[h_idx * nbw].y= 0.0;
        }
        else
        {
            hh[(h_idx * nbw)].x = 0.0;
            hh[h_idx * nbw].y =0.0;
        }

     }
}
extern "C" void launch_my_pack_c_kernel_complex_single(const int row_count, const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int stripe_count, const int l_nev, cuFloatComplex *a_dev, cuFloatComplex *row_group_dev)

{
    cudaError_t err;
    dim3 grid_size = dim3(row_count, stripe_count, 1);
    int blocksize = stripe_width > MAX_BLOCK_SIZE ? MAX_BLOCK_SIZE : stripe_width;

    for (int i_off = 0; i_off < stripe_width / blocksize; i_off++)
    {
        my_pack_c_kernel_complex_single<<<grid_size, blocksize>>>(n_offset, max_idx, stripe_width, a_dim2, l_nev, a_dev, row_group_dev, i_off);

    }

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n my pack_kernel failed %s \n", cudaGetErrorString(err));
    }
}
extern "C" void launch_extract_hh_tau_c_kernel_complex_single(cuFloatComplex *bcast_buffer_dev, cuFloatComplex *hh_tau_dev, const int nbw, const int n, const int is_zero)

{
    cudaError_t err;
    int grid_size = 1 + (n - 1) / MAX_BLOCK_SIZE;
    extract_hh_tau_c_kernel_complex_single<<<grid_size, MAX_BLOCK_SIZE>>>(bcast_buffer_dev, hh_tau_dev, nbw, n, is_zero);

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n extract _kernel failed %s \n", cudaGetErrorString(err));
    }
}
extern "C" void launch_my_unpack_c_kernel_complex_single(const int row_count, const int n_offset, const int max_idx, const int stripe_width, const int a_dim2, const int stripe_count, const int l_nev, cuFloatComplex *row_group_dev, cuFloatComplex *a_dev)

{
    cudaError_t err;
    dim3 grid_size = dim3(row_count, stripe_count, 1);
    int blocksize = stripe_width > MAX_BLOCK_SIZE ? MAX_BLOCK_SIZE : stripe_width;

    for (int i_off = 0; i_off < stripe_width / blocksize; i_off++)
    {
        my_unpack_c_kernel_complex_single<<<grid_size, blocksize>>>(n_offset, max_idx, stripe_width, a_dim2, l_nev, row_group_dev, a_dev, i_off);

    }

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("\n my_unpack_c_kernel failed %s \n", cudaGetErrorString(err));
    }
}
