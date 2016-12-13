
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel
void matvec_single(__global const double *A, __global const double *x,
                   __global double *b, int n) {

  int gid = get_global_id(0);
  const int stride = get_global_size(0);

  double sum;
  while(gid<n) {

    sum = 0.0;
    for(int i=0; i<n; ++i)
      sum += A[gid*n+i]*x[i];
    b[gid] = sum;

    gid += stride;

  }

}

// A is nxn, M is nxk, B is nxk
__kernel
void matmat_single(__global const double *A, __global const double *M,
                   __global double *B, int n, int k) {

  int gid = get_global_id(0);
  const int stride = get_global_size(0);

  double sum;
  while(gid<n) {

    for(int j=0; j<k; ++j) {

      sum = 0.0;
      for(int i=0; i<n; ++i)
        sum += A[gid*n+i]*M[k*i+j];

      B[gid*k+j] = sum;

    }

    gid += stride;

  }

}

