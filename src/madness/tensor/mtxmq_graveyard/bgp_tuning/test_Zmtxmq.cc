/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <complex>

typedef std::complex<double> double_complex;

#ifdef __bgp__
extern "C" {
#include "dcmf.h"
}
#endif

#include "rdtsc.h"

#ifdef HPM
#include <mpi.h>
extern "C" void HPM_Init(void);           // initialize the UPC unit
extern "C" void HPM_Start(char *label);   // start counting in a block marked by the label
extern "C" void HPM_Stop(char *label);    // stop counting in a block marked by the label
extern "C" void HPM_Print(void);          // print counters for all blocks
extern "C" void HPM_Print_Flops(void);
#endif

#define ZGEMM_FNAME zgemm_

extern "C" void ZGEMM_FNAME(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                            const double_complex *alpha, const double_complex *a, const int *lda,
                            const double_complex *b, const int *ldb, const double_complex *beta,
                            double_complex *c, const int *ldc, int la, int lb);

void mTxmq(long dimi, long dimj, long dimk, double_complex* c, const double_complex* a, const double_complex* b);
void mTxm_tune(long dimi, long dimj, long dimk, double_complex* c, const double_complex* a, const double_complex* b);

void mTxm_zgemm(long ni, long nj, long nk, double_complex* c, const double_complex* a, const double_complex*b ) {
  int fni=ni;
  int fnj=nj;
  int fnk=nk;
  double_complex one=1.0;
  ZGEMM_FNAME("n","t",&fnj,&fni,&fnk,&one,b,&fnj,a,&fni,&one,c,&fnj,1,1);
}

double_complex ran()
{
  static unsigned long seed = 76521;
  seed = seed *1812433253 + 12345;
  return ((double_complex) (seed & 0x7fffffff)) * 4.6566128752458e-10;
}

void ran_fill(int n, double_complex *a) {
    while (n--) *a++ = ran();
}

void mTxm(long dimi, long dimj, long dimk,
          double_complex* c, const double_complex* a, const double_complex* b) {
    int i, j, k;
    for (k=0; k<dimk; ++k) {
        for (j=0; j<dimj; ++j) {
            for (i=0; i<dimi; ++i) {
                c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
            }
        }
    }
}

void timer(const char* s, long ni, long nj, long nk, double_complex *a, double_complex *b, double_complex *c) {
  double fastest=0.0, fastest_zgemm=0.0, fastest_tune=0.0;

  double nflop = 2.0*ni*nj*nk;
  long loop;

#ifdef HPM
  HPM_Start("mTxmq");
#endif
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxmq(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    if (rate > fastest) fastest = rate;
  }
#ifdef HPM
  HPM_Stop("mTxmq");
#endif

#ifdef HPM
  HPM_Start("mTxmq_zgemm");
#endif
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxm_zgemm(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    if (rate > fastest_zgemm) fastest_zgemm = rate;
  }
#ifdef HPM
  HPM_Stop("mTxmq_zgemm");
#endif

#ifdef HPM
  HPM_Start("mTxmq_tune");
#endif
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxm_tune(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    if (rate > fastest_tune) fastest_tune = rate;
  }
#ifdef HPM
  HPM_Stop("mTxmq_tune");
#endif

  printf("%20s %3ld %3ld %3ld %8.2f %8.2f %8.2f\n",s, ni,nj,nk, fastest, fastest_zgemm, fastest_tune);
}

void trantimer(const char* s, long ni, long nj, long nk, double_complex *a, double_complex *b, double_complex *c) {
  double fastest=0.0, fastest_zgemm=0.0, fastest_tune=0.0;

  double nflop = 3.0*2.0*ni*nj*nk;
  long loop;

#ifdef HPM
  HPM_Start("mTxmq");
#endif
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxmq(ni,nj,nk,c,a,b);
    mTxmq(ni,nj,nk,a,c,b);
    mTxmq(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    if (rate > fastest) fastest = rate;
  }
#ifdef HPM
  HPM_Stop("mTxmq");
#endif

#ifdef HPM
  HPM_Start("mTxmq_zgemm");
#endif
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxm_zgemm(ni,nj,nk,c,a,b);
    mTxm_zgemm(ni,nj,nk,a,c,b);
    mTxm_zgemm(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    if (rate > fastest_zgemm) fastest_zgemm = rate;
  }
#ifdef HPM
  HPM_Stop("mTxmq_zgemm");
#endif

#ifdef HPM
  HPM_Start("mTxmq_tune");
#endif
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxm_tune(ni,nj,nk,c,a,b);
    mTxm_tune(ni,nj,nk,a,c,b);
    mTxm_tune(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    if (rate > fastest_tune) fastest_tune = rate;
  }
#ifdef HPM
  HPM_Stop("mTxmq_tune");
#endif

  printf("%20s %3ld %3ld %3ld %8.2f %8.2f %8.2f\n",s, ni,nj,nk, fastest, fastest_zgemm, fastest_tune);
}

int main(int argc, char **argv) {
    const long nimax=30*30;
    const long njmax=100;
    const long nkmax=100;
    long ni, nj, nk, i, m;

    double_complex *a, *b, *c, *d;

#ifdef HPM
    MPI_Init(&argc, &argv);
    HPM_Init();
#endif

    posix_memalign((void **) &a, 16, nkmax*nimax*sizeof(double_complex));
    posix_memalign((void **) &b, 16, nkmax*njmax*sizeof(double_complex));
    posix_memalign((void **) &c, 16, nimax*njmax*sizeof(double_complex));
    posix_memalign((void **) &d, 16, nimax*njmax*sizeof(double_complex));

    ran_fill(nkmax*nimax, a);
    ran_fill(nkmax*njmax, b);

/*     ni = nj = nk = 2; */
/*     for (i=0; i<ni*nj; ++i) d[i] = c[i] = 0.0; */
/*     mTxm (ni,nj,nk,c,a,b); */
/*     mTxmq(ni,nj,nk,d,a,b); */
/*     for (i=0; i<ni; ++i) { */
/*       long j; */
/*       for (j=0; j<nj; ++j) { */
/* 	printf("%2ld %2ld %.6f %.6f\n", i, j, c[i*nj+j], d[i*nj+j]); */
/*       } */
/*     } */
/*     return 0; */

    printf("Starting to test Zmtxmq ... \n");
    for (ni=2; ni<60; ni+=2) {
        for (nj=2; nj<100; nj+=6) {
            for (nk=2; nk<100; nk+=6) {
                for (i=0; i<ni*nj; ++i) d[i] = c[i] = 0.0;
                mTxm (ni,nj,nk,c,a,b);
                mTxmq(ni,nj,nk,d,a,b);
                for (i=0; i<ni*nj; ++i) {
                    double err = std::abs(d[i]-c[i]);
                    /* This test is sensitive to the compilation options.
                       Be sure to have the reference code above compiled
                       -msse2 -fpmath=sse if using GCC.  Otherwise, to
                       pass the test you may need to change the threshold
                       to circa 1e-13.
                    */
                    if (err > 1e-15) {
                        printf("test_Zmtxmq: error %ld %ld %ld %e\n",ni,nj,nk,err);
                        exit(1);
                    }
                }
            }
        }
    }
    printf("... OK!\n");

    for (ni=2; ni<60; ni+=2) timer("(m*m)T*(m*m)", ni,ni,ni,a,b,c);
    for (m=2; m<=30; m+=2) timer("(m*m,m)T*(m*m)", m*m,m,m,a,b,c);
    for (m=2; m<=30; m+=2) trantimer("tran(m,m,m)", m*m,m,m,a,b,c);
    for (m=2; m<=20; m+=2) timer("(20*20,20)T*(20,m)", 20*20,m,20,a,b,c);

#ifdef HPM
    HPM_Print();
    HPM_Print_Flops();
    MPI_Finalize();
#endif

    return 0;
}
