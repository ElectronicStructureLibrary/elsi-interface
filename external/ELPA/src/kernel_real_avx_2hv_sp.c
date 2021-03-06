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
// This file contains the compute intensive kernels for the Householder transformations.
// It should be compiled with the highest possible optimization level.
//
// On Intel Nehalem or Intel Westmere or AMD Magny Cours use -O3 -msse3
// On Intel Sandy Bridge use -O3 -mavx
//
// Copyright of the original code rests with the authors inside the ELPA
// consortium. The copyright of any additional modifications shall rest
// with their original authors, but shall adhere to the licensing terms
// distributed along with the original code in the file "COPYING".
//
// Author: Alexander Heinecke (alexander.heinecke@mytum.de)
// Adapted for building a shared-library by Andreas Marek, MPCDF (andreas.marek@mpcdf.mpg.de)
// --------------------------------------------------------------------------------------------------

#include <x86intrin.h>
#include <stdio.h>
#include <stdlib.h>

#define __forceinline __attribute__((always_inline)) static

#define offset 8
#define __AVX_DATATYPE __m256
#define _AVX_LOAD _mm256_load_ps
#define _AVX_STORE _mm256_store_ps
#define _AVX_ADD _mm256_add_ps
#define _AVX_MUL _mm256_mul_ps
#define _AVX_BROADCAST _mm256_broadcast_ss
#define _AVX_XOR _mm256_xor_ps

#define _AVX_FMA _mm256_FMA_ps

//Forward declaration
__forceinline void hh_trafo_kernel_8_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s);
__forceinline void hh_trafo_kernel_16_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s);
__forceinline void hh_trafo_kernel_24_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s);
__forceinline void hh_trafo_kernel_32_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s);
__forceinline void hh_trafo_kernel_40_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s);
__forceinline void hh_trafo_kernel_48_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s);

void double_hh_trafo_real_avx_avx2_2hv_single(float* q, float* hh, int* pnb, int* pnq, int* pldq, int* pldh);

void double_hh_trafo_real_avx_avx2_2hv_single(float* q, float* hh, int* pnb, int* pnq, int* pldq, int* pldh)
{
        int i;
        int nb = *pnb;
        int nq = *pldq;
        int ldq = *pldq;
        int ldh = *pldh;
        int worked_on;

        worked_on = 0;

        // calculating scalar product to compute
        // 2 householder vectors simultaneously
        float s = hh[(ldh)+1]*1.0f;
        #pragma ivdep
        for (i = 2; i < nb; i++)
        {
                s += hh[i-1] * hh[(i+ldh)];
        }

        // Production level kernel calls with padding
        for (i = 0; i < nq-40; i+=48)
        {
                hh_trafo_kernel_48_AVX_2hv_single(&q[i], hh, nb, ldq, ldh, s);
                worked_on += 48;
        }

        if (nq == i)
        {
                return;
        }

        if (nq-i == 40)
        {
                hh_trafo_kernel_40_AVX_2hv_single(&q[i], hh, nb, ldq, ldh, s);
                worked_on += 40;
        }

        if (nq-i == 32)
        {
                hh_trafo_kernel_32_AVX_2hv_single(&q[i], hh, nb, ldq, ldh, s);
                worked_on += 32;
        }

        if (nq-i == 24)
        {
                hh_trafo_kernel_24_AVX_2hv_single(&q[i], hh, nb, ldq, ldh, s);
                worked_on += 24;
        }

        if (nq-i == 16)
        {
                hh_trafo_kernel_16_AVX_2hv_single(&q[i], hh, nb, ldq, ldh, s);
                worked_on += 16;
        }

        if (nq-i == 8)
        {
                hh_trafo_kernel_8_AVX_2hv_single(&q[i], hh, nb, ldq, ldh, s);
                worked_on += 8;
        }
}

/**
 * Unrolled kernel that computes
 * 48 rows of Q simultaneously, a
 * matrix Vector product with two householder
 * vectors + a rank 2 update is performed
 */
 __forceinline void hh_trafo_kernel_48_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s)
{
        /////////////////////////////////////////////////////
        // Matrix Vector Multiplication, Q [24 x nb+1] * hh
        // hh contains two householder vectors, with offset 1
        /////////////////////////////////////////////////////
        int i;
        // Needed bit mask for floating point sign flip
        __AVX_DATATYPE sign = (__AVX_DATATYPE)_mm256_set1_epi32(0x80000000);

        __AVX_DATATYPE x1 = _AVX_LOAD(&q[ldq]);
        __AVX_DATATYPE x2 = _AVX_LOAD(&q[ldq+offset]);
        __AVX_DATATYPE x3 = _AVX_LOAD(&q[ldq+2*offset]);
        __AVX_DATATYPE x4 = _AVX_LOAD(&q[ldq+3*offset]);
        __AVX_DATATYPE x5 = _AVX_LOAD(&q[ldq+4*offset]);
        __AVX_DATATYPE x6 = _AVX_LOAD(&q[ldq+5*offset]);

        __AVX_DATATYPE h1 = _AVX_BROADCAST(&hh[ldh+1]);
        __AVX_DATATYPE h2;

        __AVX_DATATYPE q1 = _AVX_LOAD(q);
        __AVX_DATATYPE y1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        __AVX_DATATYPE q2 = _AVX_LOAD(&q[offset]);
        __AVX_DATATYPE y2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        __AVX_DATATYPE q3 = _AVX_LOAD(&q[2*offset]);
        __AVX_DATATYPE y3 = _AVX_ADD(q3, _AVX_MUL(x3, h1));
        __AVX_DATATYPE q4 = _AVX_LOAD(&q[3*offset]);
        __AVX_DATATYPE y4 = _AVX_ADD(q4, _AVX_MUL(x4, h1));
        __AVX_DATATYPE q5 = _AVX_LOAD(&q[4*offset]);
        __AVX_DATATYPE y5 = _AVX_ADD(q5, _AVX_MUL(x5, h1));
        __AVX_DATATYPE q6 = _AVX_LOAD(&q[5*offset]);
        __AVX_DATATYPE y6 = _AVX_ADD(q6, _AVX_MUL(x6, h1));

        for(i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
                y1 = _AVX_ADD(y1, _AVX_MUL(q1,h2));
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
                y2 = _AVX_ADD(y2, _AVX_MUL(q2,h2));
                q3 = _AVX_LOAD(&q[(i*ldq)+2*offset]);
                x3 = _AVX_ADD(x3, _AVX_MUL(q3,h1));
                y3 = _AVX_ADD(y3, _AVX_MUL(q3,h2));
                q4 = _AVX_LOAD(&q[(i*ldq)+3*offset]);
                x4 = _AVX_ADD(x4, _AVX_MUL(q4,h1));
                y4 = _AVX_ADD(y4, _AVX_MUL(q4,h2));
                q5 = _AVX_LOAD(&q[(i*ldq)+4*offset]);
                x5 = _AVX_ADD(x5, _AVX_MUL(q5,h1));
                y5 = _AVX_ADD(y5, _AVX_MUL(q5,h2));
                q6 = _AVX_LOAD(&q[(i*ldq)+5*offset]);
                x6 = _AVX_ADD(x6, _AVX_MUL(q6,h1));
                y6 = _AVX_ADD(y6, _AVX_MUL(q6,h2));

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
        q3 = _AVX_LOAD(&q[(nb*ldq)+2*offset]);
        x3 = _AVX_ADD(x3, _AVX_MUL(q3,h1));
        q4 = _AVX_LOAD(&q[(nb*ldq)+3*offset]);
        x4 = _AVX_ADD(x4, _AVX_MUL(q4,h1));
        q5 = _AVX_LOAD(&q[(nb*ldq)+4*offset]);
        x5 = _AVX_ADD(x5, _AVX_MUL(q5,h1));
        q6 = _AVX_LOAD(&q[(nb*ldq)+5*offset]);
        x6 = _AVX_ADD(x6, _AVX_MUL(q6,h1));

        /////////////////////////////////////////////////////
        // Rank-2 update of Q [24 x nb+1]
        /////////////////////////////////////////////////////

        __AVX_DATATYPE tau1 = _AVX_BROADCAST(hh);
        __AVX_DATATYPE tau2 = _AVX_BROADCAST(&hh[ldh]);
        __AVX_DATATYPE vs = _AVX_BROADCAST(&s);

        h1 = _AVX_XOR(tau1, sign);
        x1 = _AVX_MUL(x1, h1);
        x2 = _AVX_MUL(x2, h1);
        x3 = _AVX_MUL(x3, h1);
        x4 = _AVX_MUL(x4, h1);
        x5 = _AVX_MUL(x5, h1);
        x6 = _AVX_MUL(x6, h1);

        h1 = _AVX_XOR(tau2, sign);
        h2 = _AVX_MUL(h1, vs);
        y1 = _AVX_ADD(_AVX_MUL(y1,h1), _AVX_MUL(x1,h2));
        y2 = _AVX_ADD(_AVX_MUL(y2,h1), _AVX_MUL(x2,h2));
        y3 = _AVX_ADD(_AVX_MUL(y3,h1), _AVX_MUL(x3,h2));
        y4 = _AVX_ADD(_AVX_MUL(y4,h1), _AVX_MUL(x4,h2));
        y5 = _AVX_ADD(_AVX_MUL(y5,h1), _AVX_MUL(x5,h2));
        y6 = _AVX_ADD(_AVX_MUL(y6,h1), _AVX_MUL(x6,h2));

        q1 = _AVX_LOAD(q);
        q1 = _AVX_ADD(q1, y1);
        _AVX_STORE(q,q1);
        q2 = _AVX_LOAD(&q[offset]);
        q2 = _AVX_ADD(q2, y2);
        _AVX_STORE(&q[offset],q2);
        q3 = _AVX_LOAD(&q[2*offset]);
        q3 = _AVX_ADD(q3, y3);
        _AVX_STORE(&q[2*offset],q3);
        q4 = _AVX_LOAD(&q[3*offset]);
        q4 = _AVX_ADD(q4, y4);
        _AVX_STORE(&q[3*offset],q4);
        q5 = _AVX_LOAD(&q[4*offset]);
        q5 = _AVX_ADD(q5, y5);
        _AVX_STORE(&q[4*offset],q5);
        q6 = _AVX_LOAD(&q[5*offset]);
        q6 = _AVX_ADD(q6, y6);
        _AVX_STORE(&q[5*offset],q6);

        h2 = _AVX_BROADCAST(&hh[ldh+1]);
        q1 = _AVX_LOAD(&q[ldq]);
        q1 = _AVX_ADD(q1, _AVX_ADD(x1, _AVX_MUL(y1, h2)));
        _AVX_STORE(&q[ldq],q1);
        q2 = _AVX_LOAD(&q[ldq+offset]);
        q2 = _AVX_ADD(q2, _AVX_ADD(x2, _AVX_MUL(y2, h2)));
        _AVX_STORE(&q[ldq+offset],q2);
        q3 = _AVX_LOAD(&q[ldq+2*offset]);
        q3 = _AVX_ADD(q3, _AVX_ADD(x3, _AVX_MUL(y3, h2)));
        _AVX_STORE(&q[ldq+2*offset],q3);
        q4 = _AVX_LOAD(&q[ldq+3*offset]);
        q4 = _AVX_ADD(q4, _AVX_ADD(x4, _AVX_MUL(y4, h2)));
        _AVX_STORE(&q[ldq+3*offset],q4);
        q5 = _AVX_LOAD(&q[ldq+4*offset]);
        q5 = _AVX_ADD(q5, _AVX_ADD(x5, _AVX_MUL(y5, h2)));
        _AVX_STORE(&q[ldq+4*offset],q5);
        q6 = _AVX_LOAD(&q[ldq+5*offset]);
        q6 = _AVX_ADD(q6, _AVX_ADD(x6, _AVX_MUL(y6, h2)));
        _AVX_STORE(&q[ldq+5*offset],q6);

        for (i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                q1 = _AVX_ADD(q1, _AVX_ADD(_AVX_MUL(x1,h1), _AVX_MUL(y1, h2)));
                _AVX_STORE(&q[i*ldq],q1);
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                q2 = _AVX_ADD(q2, _AVX_ADD(_AVX_MUL(x2,h1), _AVX_MUL(y2, h2)));
                _AVX_STORE(&q[(i*ldq)+offset],q2);
                q3 = _AVX_LOAD(&q[(i*ldq)+2*offset]);
                q3 = _AVX_ADD(q3, _AVX_ADD(_AVX_MUL(x3,h1), _AVX_MUL(y3, h2)));
                _AVX_STORE(&q[(i*ldq)+2*offset],q3);
                q4 = _AVX_LOAD(&q[(i*ldq)+3*offset]);
                q4 = _AVX_ADD(q4, _AVX_ADD(_AVX_MUL(x4,h1), _AVX_MUL(y4, h2)));
                _AVX_STORE(&q[(i*ldq)+3*offset],q4);
                q5 = _AVX_LOAD(&q[(i*ldq)+4*offset]);
                q5 = _AVX_ADD(q5, _AVX_ADD(_AVX_MUL(x5,h1), _AVX_MUL(y5, h2)));
                _AVX_STORE(&q[(i*ldq)+4*offset],q5);
                q6 = _AVX_LOAD(&q[(i*ldq)+5*offset]);
                q6 = _AVX_ADD(q6, _AVX_ADD(_AVX_MUL(x6,h1), _AVX_MUL(y6, h2)));
                _AVX_STORE(&q[(i*ldq)+5*offset],q6);

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        q1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        _AVX_STORE(&q[nb*ldq],q1);
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        q2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        _AVX_STORE(&q[(nb*ldq)+offset],q2);
        q3 = _AVX_LOAD(&q[(nb*ldq)+2*offset]);
        q3 = _AVX_ADD(q3, _AVX_MUL(x3, h1));
        _AVX_STORE(&q[(nb*ldq)+2*offset],q3);
        q4 = _AVX_LOAD(&q[(nb*ldq)+3*offset]);
        q4 = _AVX_ADD(q4, _AVX_MUL(x4, h1));
        _AVX_STORE(&q[(nb*ldq)+3*offset],q4);
        q5 = _AVX_LOAD(&q[(nb*ldq)+4*offset]);
        q5 = _AVX_ADD(q5, _AVX_MUL(x5, h1));
        _AVX_STORE(&q[(nb*ldq)+4*offset],q5);
        q6 = _AVX_LOAD(&q[(nb*ldq)+5*offset]);
        q6 = _AVX_ADD(q6, _AVX_MUL(x6, h1));
        _AVX_STORE(&q[(nb*ldq)+5*offset],q6);

}

/**
 * Unrolled kernel that computes
 * 40 rows of Q simultaneously, a
 * matrix Vector product with two householder
 * vectors + a rank 2 update is performed
 */
 __forceinline void hh_trafo_kernel_40_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s)
{
        /////////////////////////////////////////////////////
        // Matrix Vector Multiplication, Q [24 x nb+1] * hh
        // hh contains two householder vectors, with offset 1
        /////////////////////////////////////////////////////
        int i;
        // Needed bit mask for floating point sign flip
        __AVX_DATATYPE sign = (__AVX_DATATYPE)_mm256_set1_epi32(0x80000000);

        __AVX_DATATYPE x1 = _AVX_LOAD(&q[ldq]);
        __AVX_DATATYPE x2 = _AVX_LOAD(&q[ldq+offset]);
        __AVX_DATATYPE x3 = _AVX_LOAD(&q[ldq+2*offset]);
        __AVX_DATATYPE x4 = _AVX_LOAD(&q[ldq+3*offset]);
        __AVX_DATATYPE x5 = _AVX_LOAD(&q[ldq+4*offset]);

        __AVX_DATATYPE h1 = _AVX_BROADCAST(&hh[ldh+1]);
        __AVX_DATATYPE h2;

        __AVX_DATATYPE q1 = _AVX_LOAD(q);
        __AVX_DATATYPE y1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        __AVX_DATATYPE q2 = _AVX_LOAD(&q[offset]);
        __AVX_DATATYPE y2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        __AVX_DATATYPE q3 = _AVX_LOAD(&q[2*offset]);
        __AVX_DATATYPE y3 = _AVX_ADD(q3, _AVX_MUL(x3, h1));
        __AVX_DATATYPE q4 = _AVX_LOAD(&q[3*offset]);
        __AVX_DATATYPE y4 = _AVX_ADD(q4, _AVX_MUL(x4, h1));
        __AVX_DATATYPE q5 = _AVX_LOAD(&q[4*offset]);
        __AVX_DATATYPE y5 = _AVX_ADD(q5, _AVX_MUL(x5, h1));

        for(i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
                y1 = _AVX_ADD(y1, _AVX_MUL(q1,h2));
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
                y2 = _AVX_ADD(y2, _AVX_MUL(q2,h2));
                q3 = _AVX_LOAD(&q[(i*ldq)+2*offset]);
                x3 = _AVX_ADD(x3, _AVX_MUL(q3,h1));
                y3 = _AVX_ADD(y3, _AVX_MUL(q3,h2));
                q4 = _AVX_LOAD(&q[(i*ldq)+3*offset]);
                x4 = _AVX_ADD(x4, _AVX_MUL(q4,h1));
                y4 = _AVX_ADD(y4, _AVX_MUL(q4,h2));
                q5 = _AVX_LOAD(&q[(i*ldq)+4*offset]);
                x5 = _AVX_ADD(x5, _AVX_MUL(q5,h1));
                y5 = _AVX_ADD(y5, _AVX_MUL(q5,h2));

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
        q3 = _AVX_LOAD(&q[(nb*ldq)+2*offset]);
        x3 = _AVX_ADD(x3, _AVX_MUL(q3,h1));
        q4 = _AVX_LOAD(&q[(nb*ldq)+3*offset]);
        x4 = _AVX_ADD(x4, _AVX_MUL(q4,h1));
        q5 = _AVX_LOAD(&q[(nb*ldq)+4*offset]);
        x5 = _AVX_ADD(x5, _AVX_MUL(q5,h1));

        /////////////////////////////////////////////////////
        // Rank-2 update of Q [24 x nb+1]
        /////////////////////////////////////////////////////

        __AVX_DATATYPE tau1 = _AVX_BROADCAST(hh);
        __AVX_DATATYPE tau2 = _AVX_BROADCAST(&hh[ldh]);
        __AVX_DATATYPE vs = _AVX_BROADCAST(&s);

        h1 = _AVX_XOR(tau1, sign);
        x1 = _AVX_MUL(x1, h1);
        x2 = _AVX_MUL(x2, h1);
        x3 = _AVX_MUL(x3, h1);
        x4 = _AVX_MUL(x4, h1);
        x5 = _AVX_MUL(x5, h1);

        h1 = _AVX_XOR(tau2, sign);
        h2 = _AVX_MUL(h1, vs);
        y1 = _AVX_ADD(_AVX_MUL(y1,h1), _AVX_MUL(x1,h2));
        y2 = _AVX_ADD(_AVX_MUL(y2,h1), _AVX_MUL(x2,h2));
        y3 = _AVX_ADD(_AVX_MUL(y3,h1), _AVX_MUL(x3,h2));
        y4 = _AVX_ADD(_AVX_MUL(y4,h1), _AVX_MUL(x4,h2));
        y5 = _AVX_ADD(_AVX_MUL(y5,h1), _AVX_MUL(x5,h2));

        q1 = _AVX_LOAD(q);
        q1 = _AVX_ADD(q1, y1);
        _AVX_STORE(q,q1);
        q2 = _AVX_LOAD(&q[offset]);
        q2 = _AVX_ADD(q2, y2);
        _AVX_STORE(&q[offset],q2);
        q3 = _AVX_LOAD(&q[2*offset]);
        q3 = _AVX_ADD(q3, y3);
        _AVX_STORE(&q[2*offset],q3);
        q4 = _AVX_LOAD(&q[3*offset]);
        q4 = _AVX_ADD(q4, y4);
        _AVX_STORE(&q[3*offset],q4);
        q5 = _AVX_LOAD(&q[4*offset]);
        q5 = _AVX_ADD(q5, y5);
        _AVX_STORE(&q[4*offset],q5);

        h2 = _AVX_BROADCAST(&hh[ldh+1]);
        q1 = _AVX_LOAD(&q[ldq]);
        q1 = _AVX_ADD(q1, _AVX_ADD(x1, _AVX_MUL(y1, h2)));
        _AVX_STORE(&q[ldq],q1);
        q2 = _AVX_LOAD(&q[ldq+offset]);
        q2 = _AVX_ADD(q2, _AVX_ADD(x2, _AVX_MUL(y2, h2)));
        _AVX_STORE(&q[ldq+offset],q2);
        q3 = _AVX_LOAD(&q[ldq+2*offset]);
        q3 = _AVX_ADD(q3, _AVX_ADD(x3, _AVX_MUL(y3, h2)));
        _AVX_STORE(&q[ldq+2*offset],q3);
        q4 = _AVX_LOAD(&q[ldq+3*offset]);
        q4 = _AVX_ADD(q4, _AVX_ADD(x4, _AVX_MUL(y4, h2)));
        _AVX_STORE(&q[ldq+3*offset],q4);
        q5 = _AVX_LOAD(&q[ldq+4*offset]);
        q5 = _AVX_ADD(q5, _AVX_ADD(x5, _AVX_MUL(y5, h2)));
        _AVX_STORE(&q[ldq+4*offset],q5);

        for (i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                q1 = _AVX_ADD(q1, _AVX_ADD(_AVX_MUL(x1,h1), _AVX_MUL(y1, h2)));
                _AVX_STORE(&q[i*ldq],q1);
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                q2 = _AVX_ADD(q2, _AVX_ADD(_AVX_MUL(x2,h1), _AVX_MUL(y2, h2)));
                _AVX_STORE(&q[(i*ldq)+offset],q2);
                q3 = _AVX_LOAD(&q[(i*ldq)+2*offset]);
                q3 = _AVX_ADD(q3, _AVX_ADD(_AVX_MUL(x3,h1), _AVX_MUL(y3, h2)));
                _AVX_STORE(&q[(i*ldq)+2*offset],q3);
                q4 = _AVX_LOAD(&q[(i*ldq)+3*offset]);
                q4 = _AVX_ADD(q4, _AVX_ADD(_AVX_MUL(x4,h1), _AVX_MUL(y4, h2)));
                _AVX_STORE(&q[(i*ldq)+3*offset],q4);
                q5 = _AVX_LOAD(&q[(i*ldq)+4*offset]);
                q5 = _AVX_ADD(q5, _AVX_ADD(_AVX_MUL(x5,h1), _AVX_MUL(y5, h2)));
                _AVX_STORE(&q[(i*ldq)+4*offset],q5);

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        q1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        _AVX_STORE(&q[nb*ldq],q1);
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        q2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        _AVX_STORE(&q[(nb*ldq)+offset],q2);
        q3 = _AVX_LOAD(&q[(nb*ldq)+2*offset]);
        q3 = _AVX_ADD(q3, _AVX_MUL(x3, h1));
        _AVX_STORE(&q[(nb*ldq)+2*offset],q3);
        q4 = _AVX_LOAD(&q[(nb*ldq)+3*offset]);
        q4 = _AVX_ADD(q4, _AVX_MUL(x4, h1));
        _AVX_STORE(&q[(nb*ldq)+3*offset],q4);
        q5 = _AVX_LOAD(&q[(nb*ldq)+4*offset]);
        q5 = _AVX_ADD(q5, _AVX_MUL(x5, h1));
        _AVX_STORE(&q[(nb*ldq)+4*offset],q5);

}

/**
 * Unrolled kernel that computes
 * 32 rows of Q simultaneously, a
 * matrix Vector product with two householder
 * vectors + a rank 2 update is performed
 */
 __forceinline void hh_trafo_kernel_32_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s)
{
        /////////////////////////////////////////////////////
        // Matrix Vector Multiplication, Q [16 x nb+1] * hh
        // hh contains two householder vectors, with offset 1
        /////////////////////////////////////////////////////
        int i;
        // Needed bit mask for floating point sign flip
        __AVX_DATATYPE sign = (__AVX_DATATYPE)_mm256_set1_epi32(0x80000000);

        __AVX_DATATYPE x1 = _AVX_LOAD(&q[ldq]);
        __AVX_DATATYPE x2 = _AVX_LOAD(&q[ldq+offset]);
        __AVX_DATATYPE x3 = _AVX_LOAD(&q[ldq+2*offset]);
        __AVX_DATATYPE x4 = _AVX_LOAD(&q[ldq+3*offset]);

        __AVX_DATATYPE h1 = _AVX_BROADCAST(&hh[ldh+1]);
        __AVX_DATATYPE h2;

        __AVX_DATATYPE q1 = _AVX_LOAD(q);
        __AVX_DATATYPE y1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        __AVX_DATATYPE q2 = _AVX_LOAD(&q[offset]);
        __AVX_DATATYPE y2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        __AVX_DATATYPE q3 = _AVX_LOAD(&q[2*offset]);
        __AVX_DATATYPE y3 = _AVX_ADD(q3, _AVX_MUL(x3, h1));
        __AVX_DATATYPE q4 = _AVX_LOAD(&q[3*offset]);
        __AVX_DATATYPE y4 = _AVX_ADD(q4, _AVX_MUL(x4, h1));

        for(i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
                y1 = _AVX_ADD(y1, _AVX_MUL(q1,h2));
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
                y2 = _AVX_ADD(y2, _AVX_MUL(q2,h2));
                q3 = _AVX_LOAD(&q[(i*ldq)+2*offset]);
                x3 = _AVX_ADD(x3, _AVX_MUL(q3,h1));
                y3 = _AVX_ADD(y3, _AVX_MUL(q3,h2));
                q4 = _AVX_LOAD(&q[(i*ldq)+3*offset]);
                x4 = _AVX_ADD(x4, _AVX_MUL(q4,h1));
                y4 = _AVX_ADD(y4, _AVX_MUL(q4,h2));

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
        q3 = _AVX_LOAD(&q[(nb*ldq)+2*offset]);
        x3 = _AVX_ADD(x3, _AVX_MUL(q3,h1));
        q4 = _AVX_LOAD(&q[(nb*ldq)+3*offset]);
        x4 = _AVX_ADD(x4, _AVX_MUL(q4,h1));

        /////////////////////////////////////////////////////
        // Rank-2 update of Q [16 x nb+1]
        /////////////////////////////////////////////////////

        __AVX_DATATYPE tau1 = _AVX_BROADCAST(hh);
        __AVX_DATATYPE tau2 = _AVX_BROADCAST(&hh[ldh]);
        __AVX_DATATYPE vs = _AVX_BROADCAST(&s);

        h1 = _AVX_XOR(tau1, sign);
        x1 = _AVX_MUL(x1, h1);
        x2 = _AVX_MUL(x2, h1);
        x3 = _AVX_MUL(x3, h1);
        x4 = _AVX_MUL(x4, h1);
        h1 = _AVX_XOR(tau2, sign);
        h2 = _AVX_MUL(h1, vs);
        y1 = _AVX_ADD(_AVX_MUL(y1,h1), _AVX_MUL(x1,h2));
        y2 = _AVX_ADD(_AVX_MUL(y2,h1), _AVX_MUL(x2,h2));
        y3 = _AVX_ADD(_AVX_MUL(y3,h1), _AVX_MUL(x3,h2));
        y4 = _AVX_ADD(_AVX_MUL(y4,h1), _AVX_MUL(x4,h2));

        q1 = _AVX_LOAD(q);
        q1 = _AVX_ADD(q1, y1);
        _AVX_STORE(q,q1);
        q2 = _AVX_LOAD(&q[offset]);
        q2 = _AVX_ADD(q2, y2);
        _AVX_STORE(&q[offset],q2);
        q3 = _AVX_LOAD(&q[2*offset]);
        q3 = _AVX_ADD(q3, y3);
        _AVX_STORE(&q[2*offset],q3);
        q4 = _AVX_LOAD(&q[3*offset]);
        q4 = _AVX_ADD(q4, y4);
        _AVX_STORE(&q[3*offset],q4);

        h2 = _AVX_BROADCAST(&hh[ldh+1]);
        q1 = _AVX_LOAD(&q[ldq]);
        q1 = _AVX_ADD(q1, _AVX_ADD(x1, _AVX_MUL(y1, h2)));
        _AVX_STORE(&q[ldq],q1);
        q2 = _AVX_LOAD(&q[ldq+offset]);
        q2 = _AVX_ADD(q2, _AVX_ADD(x2, _AVX_MUL(y2, h2)));
        _AVX_STORE(&q[ldq+offset],q2);
        q3 = _AVX_LOAD(&q[ldq+2*offset]);
        q3 = _AVX_ADD(q3, _AVX_ADD(x3, _AVX_MUL(y3, h2)));
        _AVX_STORE(&q[ldq+2*offset],q3);
        q4 = _AVX_LOAD(&q[ldq+3*offset]);
        q4 = _AVX_ADD(q4, _AVX_ADD(x4, _AVX_MUL(y4, h2)));
        _AVX_STORE(&q[ldq+3*offset],q4);

        for (i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                q1 = _AVX_ADD(q1, _AVX_ADD(_AVX_MUL(x1,h1), _AVX_MUL(y1, h2)));
                _AVX_STORE(&q[i*ldq],q1);
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                q2 = _AVX_ADD(q2, _AVX_ADD(_AVX_MUL(x2,h1), _AVX_MUL(y2, h2)));
                _AVX_STORE(&q[(i*ldq)+offset],q2);
                q3 = _AVX_LOAD(&q[(i*ldq)+2*offset]);
                q3 = _AVX_ADD(q3, _AVX_ADD(_AVX_MUL(x3,h1), _AVX_MUL(y3, h2)));
                _AVX_STORE(&q[(i*ldq)+2*offset],q3);
                q4 = _AVX_LOAD(&q[(i*ldq)+3*offset]);
                q4 = _AVX_ADD(q4, _AVX_ADD(_AVX_MUL(x4,h1), _AVX_MUL(y4, h2)));
                _AVX_STORE(&q[(i*ldq)+3*offset],q4);

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        q1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        _AVX_STORE(&q[nb*ldq],q1);
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        q2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        _AVX_STORE(&q[(nb*ldq)+offset],q2);
        q3 = _AVX_LOAD(&q[(nb*ldq)+2*offset]);
        q3 = _AVX_ADD(q3, _AVX_MUL(x3, h1));
        _AVX_STORE(&q[(nb*ldq)+2*offset],q3);
        q4 = _AVX_LOAD(&q[(nb*ldq)+3*offset]);
        q4 = _AVX_ADD(q4, _AVX_MUL(x4, h1));
        _AVX_STORE(&q[(nb*ldq)+3*offset],q4);

}

/**
 * Unrolled kernel that computes
 * 24 rows of Q simultaneously, a
 * matrix Vector product with two householder
 * vectors + a rank 2 update is performed
 */
 __forceinline void hh_trafo_kernel_24_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s)
{
        /////////////////////////////////////////////////////
        // Matrix Vector Multiplication, Q [16 x nb+1] * hh
        // hh contains two householder vectors, with offset 1
        /////////////////////////////////////////////////////
        int i;
        // Needed bit mask for floating point sign flip
        __AVX_DATATYPE sign = (__AVX_DATATYPE)_mm256_set1_epi32(0x80000000);

        __AVX_DATATYPE x1 = _AVX_LOAD(&q[ldq]);
        __AVX_DATATYPE x2 = _AVX_LOAD(&q[ldq+offset]);
        __AVX_DATATYPE x3 = _AVX_LOAD(&q[ldq+2*offset]);

        __AVX_DATATYPE h1 = _AVX_BROADCAST(&hh[ldh+1]);
        __AVX_DATATYPE h2;

        __AVX_DATATYPE q1 = _AVX_LOAD(q);
        __AVX_DATATYPE y1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        __AVX_DATATYPE q2 = _AVX_LOAD(&q[offset]);
        __AVX_DATATYPE y2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        __AVX_DATATYPE q3 = _AVX_LOAD(&q[2*offset]);
        __AVX_DATATYPE y3 = _AVX_ADD(q3, _AVX_MUL(x3, h1));

        for(i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
                y1 = _AVX_ADD(y1, _AVX_MUL(q1,h2));
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
                y2 = _AVX_ADD(y2, _AVX_MUL(q2,h2));
                q3 = _AVX_LOAD(&q[(i*ldq)+2*offset]);
                x3 = _AVX_ADD(x3, _AVX_MUL(q3,h1));
                y3 = _AVX_ADD(y3, _AVX_MUL(q3,h2));

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
        q3 = _AVX_LOAD(&q[(nb*ldq)+2*offset]);
        x3 = _AVX_ADD(x3, _AVX_MUL(q3,h1));

        /////////////////////////////////////////////////////
        // Rank-2 update of Q [16 x nb+1]
        /////////////////////////////////////////////////////

        __AVX_DATATYPE tau1 = _AVX_BROADCAST(hh);
        __AVX_DATATYPE tau2 = _AVX_BROADCAST(&hh[ldh]);
        __AVX_DATATYPE vs = _AVX_BROADCAST(&s);

        h1 = _AVX_XOR(tau1, sign);
        x1 = _AVX_MUL(x1, h1);
        x2 = _AVX_MUL(x2, h1);
        x3 = _AVX_MUL(x3, h1);
        h1 = _AVX_XOR(tau2, sign);
        h2 = _AVX_MUL(h1, vs);
        y1 = _AVX_ADD(_AVX_MUL(y1,h1), _AVX_MUL(x1,h2));
        y2 = _AVX_ADD(_AVX_MUL(y2,h1), _AVX_MUL(x2,h2));
        y3 = _AVX_ADD(_AVX_MUL(y3,h1), _AVX_MUL(x3,h2));

        q1 = _AVX_LOAD(q);
        q1 = _AVX_ADD(q1, y1);
        _AVX_STORE(q,q1);
        q2 = _AVX_LOAD(&q[offset]);
        q2 = _AVX_ADD(q2, y2);
        _AVX_STORE(&q[offset],q2);
        q3 = _AVX_LOAD(&q[2*offset]);
        q3 = _AVX_ADD(q3, y3);
        _AVX_STORE(&q[2*offset],q3);

        h2 = _AVX_BROADCAST(&hh[ldh+1]);
        q1 = _AVX_LOAD(&q[ldq]);
        q1 = _AVX_ADD(q1, _AVX_ADD(x1, _AVX_MUL(y1, h2)));
        _AVX_STORE(&q[ldq],q1);
        q2 = _AVX_LOAD(&q[ldq+offset]);
        q2 = _AVX_ADD(q2, _AVX_ADD(x2, _AVX_MUL(y2, h2)));
        _AVX_STORE(&q[ldq+offset],q2);
        q3 = _AVX_LOAD(&q[ldq+2*offset]);
        q3 = _AVX_ADD(q3, _AVX_ADD(x3, _AVX_MUL(y3, h2)));
        _AVX_STORE(&q[ldq+2*offset],q3);

        for (i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                q1 = _AVX_ADD(q1, _AVX_ADD(_AVX_MUL(x1,h1), _AVX_MUL(y1, h2)));
                _AVX_STORE(&q[i*ldq],q1);
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                q2 = _AVX_ADD(q2, _AVX_ADD(_AVX_MUL(x2,h1), _AVX_MUL(y2, h2)));
                _AVX_STORE(&q[(i*ldq)+offset],q2);
                q3 = _AVX_LOAD(&q[(i*ldq)+2*offset]);
                q3 = _AVX_ADD(q3, _AVX_ADD(_AVX_MUL(x3,h1), _AVX_MUL(y3, h2)));
                _AVX_STORE(&q[(i*ldq)+2*offset],q3);

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        q1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        _AVX_STORE(&q[nb*ldq],q1);
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        q2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        _AVX_STORE(&q[(nb*ldq)+offset],q2);
        q3 = _AVX_LOAD(&q[(nb*ldq)+2*offset]);
        q3 = _AVX_ADD(q3, _AVX_MUL(x3, h1));
        _AVX_STORE(&q[(nb*ldq)+2*offset],q3);

}

/**
 * Unrolled kernel that computes
 * matrix Vector product with two householder
 * vectors + a rank 2 update is performed
 */
 __forceinline void hh_trafo_kernel_16_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s)
{
        /////////////////////////////////////////////////////
        // Matrix Vector Multiplication, Q [8 x nb+1] * hh
        // hh contains two householder vectors, with offset 1
        /////////////////////////////////////////////////////
        int i;
        // Needed bit mask for floating point sign flip
        __AVX_DATATYPE sign = (__AVX_DATATYPE)_mm256_set1_epi32(0x80000000);

        __AVX_DATATYPE x1 = _AVX_LOAD(&q[ldq]);
        __AVX_DATATYPE x2 = _AVX_LOAD(&q[ldq+offset]);
        __AVX_DATATYPE h1 = _AVX_BROADCAST(&hh[ldh+1]);
        __AVX_DATATYPE h2;

        __AVX_DATATYPE q1 = _AVX_LOAD(q);
        __AVX_DATATYPE y1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        __AVX_DATATYPE q2 = _AVX_LOAD(&q[offset]);
        __AVX_DATATYPE y2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));

        for(i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
                y1 = _AVX_ADD(y1, _AVX_MUL(q1,h2));
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));
                y2 = _AVX_ADD(y2, _AVX_MUL(q2,h2));

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        x2 = _AVX_ADD(x2, _AVX_MUL(q2,h1));

        /////////////////////////////////////////////////////
        // Rank-2 update of Q [8 x nb+1]
        /////////////////////////////////////////////////////

        __AVX_DATATYPE tau1 = _AVX_BROADCAST(hh);
        __AVX_DATATYPE tau2 = _AVX_BROADCAST(&hh[ldh]);
        __AVX_DATATYPE vs = _AVX_BROADCAST(&s);

        h1 = _AVX_XOR(tau1, sign);
        x1 = _AVX_MUL(x1, h1);
        x2 = _AVX_MUL(x2, h1);
        h1 = _AVX_XOR(tau2, sign);
        h2 = _AVX_MUL(h1, vs);
        y1 = _AVX_ADD(_AVX_MUL(y1,h1), _AVX_MUL(x1,h2));
        y2 = _AVX_ADD(_AVX_MUL(y2,h1), _AVX_MUL(x2,h2));

        q1 = _AVX_LOAD(q);
        q1 = _AVX_ADD(q1, y1);
        _AVX_STORE(q,q1);
        q2 = _AVX_LOAD(&q[offset]);
        q2 = _AVX_ADD(q2, y2);
        _AVX_STORE(&q[offset],q2);

        h2 = _AVX_BROADCAST(&hh[ldh+1]);
        q1 = _AVX_LOAD(&q[ldq]);
        q1 = _AVX_ADD(q1, _AVX_ADD(x1, _AVX_MUL(y1, h2)));
        _AVX_STORE(&q[ldq],q1);
        q2 = _AVX_LOAD(&q[ldq+offset]);
        q2 = _AVX_ADD(q2, _AVX_ADD(x2, _AVX_MUL(y2, h2)));
        _AVX_STORE(&q[ldq+offset],q2);

        for (i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);
                q1 = _AVX_ADD(q1, _AVX_ADD(_AVX_MUL(x1,h1), _AVX_MUL(y1, h2)));
                _AVX_STORE(&q[i*ldq],q1);
                q2 = _AVX_LOAD(&q[(i*ldq)+offset]);
                q2 = _AVX_ADD(q2, _AVX_ADD(_AVX_MUL(x2,h1), _AVX_MUL(y2, h2)));
                _AVX_STORE(&q[(i*ldq)+offset],q2);

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);
        q1 = _AVX_LOAD(&q[nb*ldq]);
        q1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        _AVX_STORE(&q[nb*ldq],q1);
        q2 = _AVX_LOAD(&q[(nb*ldq)+offset]);
        q2 = _AVX_ADD(q2, _AVX_MUL(x2, h1));
        _AVX_STORE(&q[(nb*ldq)+offset],q2);

}

/**
 * Unrolled kernel that computes
 * 8 rows of Q simultaneously, a
* matrix Vector product with two householder
 * vectors + a rank 2 update is performed
 */
 __forceinline void hh_trafo_kernel_8_AVX_2hv_single(float* q, float* hh, int nb, int ldq, int ldh, float s)
{
        /////////////////////////////////////////////////////
        // Matrix Vector Multiplication, Q [4 x nb+1] * hh
        // hh contains two householder vectors, with offset 1
        /////////////////////////////////////////////////////
        int i;
        // Needed bit mask for floating point sign flip
        __m256 sign = (__m256)_mm256_set1_epi32(0x80000000);

        __AVX_DATATYPE x1 = _AVX_LOAD(&q[ldq]);
        __AVX_DATATYPE h1 = _AVX_BROADCAST(&hh[ldh+1]);
        __AVX_DATATYPE h2;

        __AVX_DATATYPE q1 = _AVX_LOAD(q);
        __AVX_DATATYPE y1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));

        for(i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);

                q1 = _AVX_LOAD(&q[i*ldq]);

                x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));
                y1 = _AVX_ADD(y1, _AVX_MUL(q1,h2));
        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);

        q1 = _AVX_LOAD(&q[nb*ldq]);
        x1 = _AVX_ADD(x1, _AVX_MUL(q1,h1));

        /////////////////////////////////////////////////////
        // Rank-2 update of Q [4 x nb+1]
        /////////////////////////////////////////////////////

        __AVX_DATATYPE tau1 = _AVX_BROADCAST(hh);
        __AVX_DATATYPE tau2 = _AVX_BROADCAST(&hh[ldh]);
        __AVX_DATATYPE vs = _AVX_BROADCAST(&s);

        h1 = _AVX_XOR(tau1, sign);
        x1 = _AVX_MUL(x1, h1);
        h1 = _AVX_XOR(tau2, sign);
        h2 = _AVX_MUL(h1, vs);
        y1 = _AVX_ADD(_AVX_MUL(y1,h1), _AVX_MUL(x1,h2));

        q1 = _AVX_LOAD(q);
        q1 = _AVX_ADD(q1, y1);
        _AVX_STORE(q,q1);
        h2 = _AVX_BROADCAST(&hh[ldh+1]);

        q1 = _AVX_LOAD(&q[ldq]);
        q1 = _AVX_ADD(q1, _AVX_ADD(x1, _AVX_MUL(y1, h2)));
        _AVX_STORE(&q[ldq],q1);

        for (i = 2; i < nb; i++)
        {
                h1 = _AVX_BROADCAST(&hh[i-1]);
                h2 = _AVX_BROADCAST(&hh[ldh+i]);
                q1 = _AVX_LOAD(&q[i*ldq]);

                q1 = _AVX_ADD(q1, _AVX_ADD(_AVX_MUL(x1,h1), _AVX_MUL(y1, h2)));
                _AVX_STORE(&q[i*ldq],q1);

        }

        h1 = _AVX_BROADCAST(&hh[nb-1]);

        q1 = _AVX_LOAD(&q[nb*ldq]);
        q1 = _AVX_ADD(q1, _AVX_MUL(x1, h1));
        _AVX_STORE(&q[nb*ldq],q1);
}
