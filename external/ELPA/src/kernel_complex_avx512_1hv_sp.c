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
// Author: Andreas Marek (andreas.marek@mpcdf.mpg.de)
// --------------------------------------------------------------------------------------------------

#include <complex.h>
#include <x86intrin.h>
#include <stdio.h>
#include <stdlib.h>

#define __forceinline __attribute__((always_inline))

#define __AVX512_DATATYPE __m512
#define _AVX512_LOAD _mm512_load_ps
#define _AVX512_STORE _mm512_store_ps
#define _AVX512_SET1 _mm512_set1_ps
#define _AVX512_MUL _mm512_mul_ps
#define _AVX512_ADD _mm512_add_ps
#define _AVX512_SHUFFLE _mm512_shuffle_ps
#define _AVX512_XOR _mm512_xor_ps
#define _AVX512_XOR_EPI _mm512_xor_epi32
#define _SHUFFLE 0xb1

#define _mm512_FMADDSUB_ps(a,b,c) _mm512_fmaddsub_ps(a,b,c)
#define _mm512_FMSUBADD_ps(a,b,c) _mm512_fmsubadd_ps(a,b,c)

#define _AVX512_FMADDSUB _mm512_FMADDSUB_ps
#define _AVX512_FMSUBADD _mm512_FMSUBADD_ps

//Forward declaration
static  __forceinline void hh_trafo_complex_kernel_48_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq);
static  __forceinline void hh_trafo_complex_kernel_40_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq);
static  __forceinline void hh_trafo_complex_kernel_32_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq);
static  __forceinline void hh_trafo_complex_kernel_24_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq);
static  __forceinline void hh_trafo_complex_kernel_16_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq);
static  __forceinline void hh_trafo_complex_kernel_8_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq);

void single_hh_trafo_complex_avx512_1hv_single(float complex* q, float complex* hh, int* pnb, int* pnq, int* pldq)
{
        int i;
        int nb = *pnb;
        int nq = *pldq;
        int ldq = *pldq;
        int worked_on;
        //int ldh = *pldh;

        worked_on = 0;

        for (i = 0; i < nq-40; i+=48)
        {
                hh_trafo_complex_kernel_48_AVX512_1hv_single(&q[i], hh, nb, ldq);
                worked_on += 48;
        }
        if (nq == i)
        {
                return;
        }

        if (nq-i == 40)
        {
                hh_trafo_complex_kernel_40_AVX512_1hv_single(&q[i], hh, nb, ldq);
                worked_on += 40;
        }

        if (nq-i == 32)
        {
                hh_trafo_complex_kernel_32_AVX512_1hv_single(&q[i], hh, nb, ldq);
                worked_on += 32;
        }

        if (nq-i == 24)
        {
                hh_trafo_complex_kernel_24_AVX512_1hv_single(&q[i], hh, nb, ldq);
                worked_on += 24;
        }

        if (nq-i == 16)
        {
                hh_trafo_complex_kernel_16_AVX512_1hv_single(&q[i], hh, nb, ldq);
                worked_on += 16;
        }

        if (nq-i == 8)
        {
                hh_trafo_complex_kernel_8_AVX512_1hv_single(&q[i], hh, nb, ldq);
                worked_on += 8;
        }
}

static __forceinline void hh_trafo_complex_kernel_48_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq)
{

        float* q_dbl = (float*)q;
        float* hh_dbl = (float*)hh;
        __AVX512_DATATYPE x1, x2, x3, x4, x5, x6;
        __AVX512_DATATYPE q1, q2, q3, q4, q5, q6;
        __AVX512_DATATYPE h1_real, h1_imag;
        __AVX512_DATATYPE tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
        int i=0;

        __AVX512_DATATYPE sign = (__AVX512_DATATYPE)_mm512_set1_epi32(0x80000000);

#define offset 16

        x1 = _AVX512_LOAD(&q_dbl[0]);    // complex 1, 2, 3, 4
        x2 = _AVX512_LOAD(&q_dbl[offset]);    // complex 5, 6, 7, 8
        x3 = _AVX512_LOAD(&q_dbl[2*offset]);   // complex 9, 10, 11, 12
        x4 = _AVX512_LOAD(&q_dbl[3*offset]);   // complex 13, 14, 15, 16
        x5 = _AVX512_LOAD(&q_dbl[4*offset]);   // complex 17, 18, 19, 20
        x6 = _AVX512_LOAD(&q_dbl[5*offset]);   // complex 21, 22, 23, 24

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);
                q3 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+2*offset]);
                q4 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+3*offset]);
                q5 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+4*offset]);
                q6 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+5*offset]);

                tmp1 = _AVX512_MUL(h1_imag, q1);

                x1 = _AVX512_ADD(x1, _AVX512_FMSUBADD(h1_real, q1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, q2);

                x2 = _AVX512_ADD(x2, _AVX512_FMSUBADD(h1_real, q2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                tmp3 = _AVX512_MUL(h1_imag, q3);

                x3 = _AVX512_ADD(x3, _AVX512_FMSUBADD(h1_real, q3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE)));

                tmp4 = _AVX512_MUL(h1_imag, q4);

                x4 = _AVX512_ADD(x4, _AVX512_FMSUBADD(h1_real, q4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE)));

                tmp5 = _AVX512_MUL(h1_imag, q5);

                x5 = _AVX512_ADD(x5, _AVX512_FMSUBADD(h1_real, q5, _AVX512_SHUFFLE(tmp5, tmp5, _SHUFFLE)));

                tmp6 = _AVX512_MUL(h1_imag, q6);

                x6 = _AVX512_ADD(x6, _AVX512_FMSUBADD(h1_real, q6, _AVX512_SHUFFLE(tmp6, tmp6, _SHUFFLE)));
        }

        h1_real = _AVX512_SET1(hh_dbl[0]);
        h1_imag = _AVX512_SET1(hh_dbl[1]);

        h1_real = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_real, (__m512i) sign);
        h1_imag = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_imag, (__m512i) sign);

        tmp1 = _AVX512_MUL(h1_imag, x1);

        x1 = _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE));

        tmp2 = _AVX512_MUL(h1_imag, x2);

        x2 = _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE));

        tmp3 = _AVX512_MUL(h1_imag, x3);

        x3 = _AVX512_FMADDSUB(h1_real, x3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE));

        tmp4 = _AVX512_MUL(h1_imag, x4);

        x4 = _AVX512_FMADDSUB(h1_real, x4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE));

        tmp5 = _AVX512_MUL(h1_imag, x5);

        x5 = _AVX512_FMADDSUB(h1_real, x5, _AVX512_SHUFFLE(tmp5, tmp5, _SHUFFLE));

        tmp6 = _AVX512_MUL(h1_imag, x6);

        x6 = _AVX512_FMADDSUB(h1_real, x6, _AVX512_SHUFFLE(tmp6, tmp6, _SHUFFLE));

        q1 = _AVX512_LOAD(&q_dbl[0]);
        q2 = _AVX512_LOAD(&q_dbl[offset]);
        q3 = _AVX512_LOAD(&q_dbl[2*offset]);
        q4 = _AVX512_LOAD(&q_dbl[3*offset]);
        q5 = _AVX512_LOAD(&q_dbl[4*offset]);
        q6 = _AVX512_LOAD(&q_dbl[5*offset]);

        q1 = _AVX512_ADD(q1, x1);
        q2 = _AVX512_ADD(q2, x2);
        q3 = _AVX512_ADD(q3, x3);
        q4 = _AVX512_ADD(q4, x4);
        q5 = _AVX512_ADD(q5, x5);
        q6 = _AVX512_ADD(q6, x6);

        _AVX512_STORE(&q_dbl[0], q1);
        _AVX512_STORE(&q_dbl[offset], q2);
        _AVX512_STORE(&q_dbl[2*offset], q3);
        _AVX512_STORE(&q_dbl[3*offset], q4);
        _AVX512_STORE(&q_dbl[4*offset], q5);
        _AVX512_STORE(&q_dbl[5*offset], q6);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);
                q3 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+2*offset]);
                q4 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+3*offset]);
                q5 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+4*offset]);
                q6 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+5*offset]);

                tmp1 = _AVX512_MUL(h1_imag, x1);

                q1 = _AVX512_ADD(q1, _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, x2);

                q2 = _AVX512_ADD(q2, _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                tmp3 = _AVX512_MUL(h1_imag, x3);

                q3 = _AVX512_ADD(q3, _AVX512_FMADDSUB(h1_real, x3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE)));

                tmp4 = _AVX512_MUL(h1_imag, x4);

                q4 = _AVX512_ADD(q4, _AVX512_FMADDSUB(h1_real, x4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE)));

                tmp5 = _AVX512_MUL(h1_imag, x5);

                q5 = _AVX512_ADD(q5, _AVX512_FMADDSUB(h1_real, x5, _AVX512_SHUFFLE(tmp5, tmp5, _SHUFFLE)));

                tmp6 = _AVX512_MUL(h1_imag, x6);

                q6 = _AVX512_ADD(q6, _AVX512_FMADDSUB(h1_real, x6, _AVX512_SHUFFLE(tmp6, tmp6, _SHUFFLE)));

                _AVX512_STORE(&q_dbl[(2*i*ldq)+0], q1);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+offset], q2);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+2*offset], q3);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+3*offset], q4);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+4*offset], q5);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+5*offset], q6);
        }
}

static __forceinline void hh_trafo_complex_kernel_40_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq)
{

        float* q_dbl = (float*)q;
        float* hh_dbl = (float*)hh;
        __AVX512_DATATYPE x1, x2, x3, x4, x5, x6;
        __AVX512_DATATYPE q1, q2, q3, q4, q5, q6;
        __AVX512_DATATYPE h1_real, h1_imag;
        __AVX512_DATATYPE tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
        int i=0;

        __AVX512_DATATYPE sign = (__AVX512_DATATYPE)_mm512_set1_epi32(0x80000000);

#define offset 16

        x1 = _AVX512_LOAD(&q_dbl[0]);    // complex 1, 2, 3, 4
        x2 = _AVX512_LOAD(&q_dbl[offset]);    // complex 5, 6, 7, 8
        x3 = _AVX512_LOAD(&q_dbl[2*offset]);   // complex 9, 10, 11, 12
        x4 = _AVX512_LOAD(&q_dbl[3*offset]);   // complex 13, 14, 15, 16
        x5 = _AVX512_LOAD(&q_dbl[4*offset]);   // complex 17, 18, 19, 20

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);
                q3 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+2*offset]);
                q4 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+3*offset]);
                q5 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+4*offset]);

                tmp1 = _AVX512_MUL(h1_imag, q1);

                x1 = _AVX512_ADD(x1, _AVX512_FMSUBADD(h1_real, q1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, q2);

                x2 = _AVX512_ADD(x2, _AVX512_FMSUBADD(h1_real, q2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                tmp3 = _AVX512_MUL(h1_imag, q3);

                x3 = _AVX512_ADD(x3, _AVX512_FMSUBADD(h1_real, q3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE)));

                tmp4 = _AVX512_MUL(h1_imag, q4);

                x4 = _AVX512_ADD(x4, _AVX512_FMSUBADD(h1_real, q4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE)));

                tmp5 = _AVX512_MUL(h1_imag, q5);

                x5 = _AVX512_ADD(x5, _AVX512_FMSUBADD(h1_real, q5, _AVX512_SHUFFLE(tmp5, tmp5, _SHUFFLE)));

        }

        h1_real = _AVX512_SET1(hh_dbl[0]);
        h1_imag = _AVX512_SET1(hh_dbl[1]);

        h1_real = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_real, (__m512i) sign);
        h1_imag = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_imag, (__m512i) sign);

        tmp1 = _AVX512_MUL(h1_imag, x1);

        x1 = _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE));

        tmp2 = _AVX512_MUL(h1_imag, x2);

        x2 = _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE));

        tmp3 = _AVX512_MUL(h1_imag, x3);

        x3 = _AVX512_FMADDSUB(h1_real, x3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE));

        tmp4 = _AVX512_MUL(h1_imag, x4);

        x4 = _AVX512_FMADDSUB(h1_real, x4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE));

        tmp5 = _AVX512_MUL(h1_imag, x5);

        x5 = _AVX512_FMADDSUB(h1_real, x5, _AVX512_SHUFFLE(tmp5, tmp5, _SHUFFLE));

        q1 = _AVX512_LOAD(&q_dbl[0]);
        q2 = _AVX512_LOAD(&q_dbl[offset]);
        q3 = _AVX512_LOAD(&q_dbl[2*offset]);
        q4 = _AVX512_LOAD(&q_dbl[3*offset]);
        q5 = _AVX512_LOAD(&q_dbl[4*offset]);

        q1 = _AVX512_ADD(q1, x1);
        q2 = _AVX512_ADD(q2, x2);
        q3 = _AVX512_ADD(q3, x3);
        q4 = _AVX512_ADD(q4, x4);
        q5 = _AVX512_ADD(q5, x5);

        _AVX512_STORE(&q_dbl[0], q1);
        _AVX512_STORE(&q_dbl[offset], q2);
        _AVX512_STORE(&q_dbl[2*offset], q3);
        _AVX512_STORE(&q_dbl[3*offset], q4);
        _AVX512_STORE(&q_dbl[4*offset], q5);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);
                q3 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+2*offset]);
                q4 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+3*offset]);
                q5 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+4*offset]);

                tmp1 = _AVX512_MUL(h1_imag, x1);

                q1 = _AVX512_ADD(q1, _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, x2);

                q2 = _AVX512_ADD(q2, _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                tmp3 = _AVX512_MUL(h1_imag, x3);

                q3 = _AVX512_ADD(q3, _AVX512_FMADDSUB(h1_real, x3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE)));

                tmp4 = _AVX512_MUL(h1_imag, x4);

                q4 = _AVX512_ADD(q4, _AVX512_FMADDSUB(h1_real, x4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE)));

                tmp5 = _AVX512_MUL(h1_imag, x5);

                q5 = _AVX512_ADD(q5, _AVX512_FMADDSUB(h1_real, x5, _AVX512_SHUFFLE(tmp5, tmp5, _SHUFFLE)));

                _AVX512_STORE(&q_dbl[(2*i*ldq)+0], q1);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+offset], q2);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+2*offset], q3);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+3*offset], q4);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+4*offset], q5);
        }
}

static __forceinline void hh_trafo_complex_kernel_32_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq)
{

        float* q_dbl = (float*)q;
        float* hh_dbl = (float*)hh;

        __AVX512_DATATYPE x1, x2, x3, x4;
        __AVX512_DATATYPE q1, q2, q3, q4;
        __AVX512_DATATYPE h1_real, h1_imag;
        __AVX512_DATATYPE tmp1, tmp2, tmp3, tmp4;
        int i=0;

        __AVX512_DATATYPE sign = (__AVX512_DATATYPE)_mm512_set1_epi32(0x80000000);

        x1 = _AVX512_LOAD(&q_dbl[0]);   // complex 1 2 3 4
        x2 = _AVX512_LOAD(&q_dbl[offset]);
        x3 = _AVX512_LOAD(&q_dbl[2*offset]);
        x4 = _AVX512_LOAD(&q_dbl[3*offset]);  // comlex 13 14 15 16

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);
                q3 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+2*offset]);
                q4 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+3*offset]);

                tmp1 = _AVX512_MUL(h1_imag, q1);

                x1 = _AVX512_ADD(x1, _AVX512_FMSUBADD(h1_real, q1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, q2);

                x2 = _AVX512_ADD(x2, _AVX512_FMSUBADD(h1_real, q2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                tmp3 = _AVX512_MUL(h1_imag, q3);

                x3 = _AVX512_ADD(x3, _AVX512_FMSUBADD(h1_real, q3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE)));

                tmp4 = _AVX512_MUL(h1_imag, q4);

                x4 = _AVX512_ADD(x4, _AVX512_FMSUBADD(h1_real, q4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE)));
        }

	h1_real = _AVX512_SET1(hh_dbl[0]);
	h1_imag = _AVX512_SET1(hh_dbl[1]);

        h1_real = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_real, (__m512i) sign);
        h1_imag = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_imag, (__m512i) sign);

        tmp1 = _AVX512_MUL(h1_imag, x1);

        x1 = _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE));

        tmp2 = _AVX512_MUL(h1_imag, x2);

        x2 = _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE));

        tmp3 = _AVX512_MUL(h1_imag, x3);

        x3 = _AVX512_FMADDSUB(h1_real, x3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE));

        tmp4 = _AVX512_MUL(h1_imag, x4);

        x4 = _AVX512_FMADDSUB(h1_real, x4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE));

        q1 = _AVX512_LOAD(&q_dbl[0]);
        q2 = _AVX512_LOAD(&q_dbl[offset]);
        q3 = _AVX512_LOAD(&q_dbl[2*offset]);
        q4 = _AVX512_LOAD(&q_dbl[3*offset]);

        q1 = _AVX512_ADD(q1, x1);
        q2 = _AVX512_ADD(q2, x2);
        q3 = _AVX512_ADD(q3, x3);
        q4 = _AVX512_ADD(q4, x4);

        _AVX512_STORE(&q_dbl[0], q1);
        _AVX512_STORE(&q_dbl[offset], q2);
        _AVX512_STORE(&q_dbl[2*offset], q3);
        _AVX512_STORE(&q_dbl[3*offset], q4);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);
                q3 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+2*offset]);
                q4 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+3*offset]);

                tmp1 = _AVX512_MUL(h1_imag, x1);

                q1 = _AVX512_ADD(q1, _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, x2);

                q2 = _AVX512_ADD(q2, _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                tmp3 = _AVX512_MUL(h1_imag, x3);

                q3 = _AVX512_ADD(q3, _AVX512_FMADDSUB(h1_real, x3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE)));

                tmp4 = _AVX512_MUL(h1_imag, x4);

                q4 = _AVX512_ADD(q4, _AVX512_FMADDSUB(h1_real, x4, _AVX512_SHUFFLE(tmp4, tmp4, _SHUFFLE)));

                _AVX512_STORE(&q_dbl[(2*i*ldq)+0], q1);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+offset], q2);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+2*offset], q3);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+3*offset], q4);
        }
}

static __forceinline void hh_trafo_complex_kernel_24_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq)
{

        float* q_dbl = (float*)q;
        float* hh_dbl = (float*)hh;

        __AVX512_DATATYPE x1, x2, x3, x4;
        __AVX512_DATATYPE q1, q2, q3, q4;
        __AVX512_DATATYPE h1_real, h1_imag;
        __AVX512_DATATYPE tmp1, tmp2, tmp3, tmp4;
        int i=0;

        __AVX512_DATATYPE sign = (__AVX512_DATATYPE)_mm512_set1_epi32(0x80000000);

        x1 = _AVX512_LOAD(&q_dbl[0]);   // complex 1 2 3 4
        x2 = _AVX512_LOAD(&q_dbl[offset]);
        x3 = _AVX512_LOAD(&q_dbl[2*offset]);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);
                q3 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+2*offset]);

                tmp1 = _AVX512_MUL(h1_imag, q1);

                x1 = _AVX512_ADD(x1, _AVX512_FMSUBADD(h1_real, q1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, q2);

                x2 = _AVX512_ADD(x2, _AVX512_FMSUBADD(h1_real, q2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                tmp3 = _AVX512_MUL(h1_imag, q3);

                x3 = _AVX512_ADD(x3, _AVX512_FMSUBADD(h1_real, q3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE)));

        }

        h1_real = _AVX512_SET1(hh_dbl[0]);
        h1_imag = _AVX512_SET1(hh_dbl[1]);

        h1_real = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_real, (__m512i) sign);
        h1_imag = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_imag, (__m512i) sign);

        tmp1 = _AVX512_MUL(h1_imag, x1);

        x1 = _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE));

        tmp2 = _AVX512_MUL(h1_imag, x2);

        x2 = _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE));

        tmp3 = _AVX512_MUL(h1_imag, x3);

        x3 = _AVX512_FMADDSUB(h1_real, x3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE));

        q1 = _AVX512_LOAD(&q_dbl[0]);
        q2 = _AVX512_LOAD(&q_dbl[offset]);
        q3 = _AVX512_LOAD(&q_dbl[2*offset]);

        q1 = _AVX512_ADD(q1, x1);
        q2 = _AVX512_ADD(q2, x2);
        q3 = _AVX512_ADD(q3, x3);

        _AVX512_STORE(&q_dbl[0], q1);
        _AVX512_STORE(&q_dbl[offset], q2);
        _AVX512_STORE(&q_dbl[2*offset], q3);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);
                q3 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+2*offset]);

                tmp1 = _AVX512_MUL(h1_imag, x1);

                q1 = _AVX512_ADD(q1, _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, x2);

                q2 = _AVX512_ADD(q2, _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                tmp3 = _AVX512_MUL(h1_imag, x3);

                q3 = _AVX512_ADD(q3, _AVX512_FMADDSUB(h1_real, x3, _AVX512_SHUFFLE(tmp3, tmp3, _SHUFFLE)));

                _AVX512_STORE(&q_dbl[(2*i*ldq)+0], q1);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+offset], q2);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+2*offset], q3);
        }
}

static __forceinline void hh_trafo_complex_kernel_16_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq)
{

        float* q_dbl = (float*)q;
        float* hh_dbl = (float*)hh;
        __AVX512_DATATYPE x1, x2;
        __AVX512_DATATYPE q1, q2;
        __AVX512_DATATYPE h1_real, h1_imag;
        __AVX512_DATATYPE tmp1, tmp2;
        int i=0;

        __AVX512_DATATYPE sign = (__AVX512_DATATYPE)_mm512_set1_epi32(0x80000000);

        x1 = _AVX512_LOAD(&q_dbl[0]);
        x2 = _AVX512_LOAD(&q_dbl[offset]);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);

                tmp1 = _AVX512_MUL(h1_imag, q1);
                x1 = _AVX512_ADD(x1, _AVX512_FMSUBADD(h1_real, q1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));
                tmp2 = _AVX512_MUL(h1_imag, q2);
                x2 = _AVX512_ADD(x2, _AVX512_FMSUBADD(h1_real, q2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));
        }

        h1_real = _AVX512_SET1(hh_dbl[0]);
        h1_imag = _AVX512_SET1(hh_dbl[1]);

        h1_real = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_real, (__m512i) sign);
        h1_imag = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_imag, (__m512i) sign);

        tmp1 = _AVX512_MUL(h1_imag, x1);
        x1 = _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE));

        tmp2 = _AVX512_MUL(h1_imag, x2);

        x2 = _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE));

        q1 = _AVX512_LOAD(&q_dbl[0]);
        q2 = _AVX512_LOAD(&q_dbl[offset]);

        q1 = _AVX512_ADD(q1, x1);
        q2 = _AVX512_ADD(q2, x2);

        _AVX512_STORE(&q_dbl[0], q1);
        _AVX512_STORE(&q_dbl[offset], q2);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);
                q2 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+offset]);

                tmp1 = _AVX512_MUL(h1_imag, x1);
                q1 = _AVX512_ADD(q1, _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                tmp2 = _AVX512_MUL(h1_imag, x2);

                q2 = _AVX512_ADD(q2, _AVX512_FMADDSUB(h1_real, x2, _AVX512_SHUFFLE(tmp2, tmp2, _SHUFFLE)));

                _AVX512_STORE(&q_dbl[(2*i*ldq)+0], q1);
                _AVX512_STORE(&q_dbl[(2*i*ldq)+offset], q2);
        }
}

static __forceinline void hh_trafo_complex_kernel_8_AVX512_1hv_single(float complex* q, float complex* hh, int nb, int ldq)
{

        float* q_dbl = (float*)q;
        float* hh_dbl = (float*)hh;
        __AVX512_DATATYPE x1, x2;
        __AVX512_DATATYPE q1, q2;
        __AVX512_DATATYPE h1_real, h1_imag;
        __AVX512_DATATYPE tmp1, tmp2;
        int i=0;

        __AVX512_DATATYPE sign = (__AVX512_DATATYPE)_mm512_set1_epi32(0x80000000);

        x1 = _AVX512_LOAD(&q_dbl[0]);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);

                tmp1 = _AVX512_MUL(h1_imag, q1);
                x1 = _AVX512_ADD(x1, _AVX512_FMSUBADD(h1_real, q1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));
        }

	h1_real = _AVX512_SET1(hh_dbl[0]);
	h1_imag = _AVX512_SET1(hh_dbl[1]);

        h1_real = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_real, (__m512i) sign);
        h1_imag = (__AVX512_DATATYPE) _AVX512_XOR_EPI((__m512i) h1_imag, (__m512i) sign);

        tmp1 = _AVX512_MUL(h1_imag, x1);
        x1 = _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE));

        q1 = _AVX512_LOAD(&q_dbl[0]);

        q1 = _AVX512_ADD(q1, x1);

        _AVX512_STORE(&q_dbl[0], q1);

        for (i = 1; i < nb; i++)
        {
                h1_real = _AVX512_SET1(hh_dbl[i*2]);
                h1_imag = _AVX512_SET1(hh_dbl[(i*2)+1]);

                q1 = _AVX512_LOAD(&q_dbl[(2*i*ldq)+0]);

                tmp1 = _AVX512_MUL(h1_imag, x1);
                q1 = _AVX512_ADD(q1, _AVX512_FMADDSUB(h1_real, x1, _AVX512_SHUFFLE(tmp1, tmp1, _SHUFFLE)));

                _AVX512_STORE(&q_dbl[(2*i*ldq)+0], q1);
        }
}
