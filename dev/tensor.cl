
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

__kernel
void matvec_tile(__global const double *A, __global const double *x,
                 __global double *b, int n, int tile_size) {

  // Assume n%tile_size==0

  // gid is i_outer
  int gid = get_global_id(0)*tile_size;
  const int stride = get_global_size(0)*tile_size;
  const int T = tile_size;

  double sum;
  while(gid<n) {

    for(int j_outer=0; j_outer<n/T; ++j_outer) {
      for(int i_inner=0; i_inner<T; ++i_inner) {

          sum = 0.0;
          for(int j_inner=0; j_inner<T; ++j_inner)
            sum += A[(gid+i_inner)*n+j_outer*T+j_inner]*x[j_outer*T+j_inner];

          b[gid+i_inner] += sum;

      }
    }

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

// A is nxn, M is nxk, B is nxk
__kernel
void matmat_tile(__global const double *A, __global const double *M,
                 __global double *B, int n, int k, int tile_size) {

  // Assume n%tile_size==0

  // gid is i_outer
  int gid = get_global_id(0)*tile_size;
  const int stride = get_global_size(0)*tile_size;
  const int T = tile_size;

  double sum;
  while(gid<n) {

    for(int j_outer=0; j_outer<n/T; ++j_outer) {
      for(int k_outer=0; k_outer<k/T; ++k_outer) {
        for(int i_inner=0; i_inner<T; ++i_inner) {
          for(int k_inner=0; k_inner<T; ++k_inner) {

            sum = 0.0;
            for(int j_inner=0; j_inner<T; ++j_inner)
              sum += A[(gid+i_inner)*n+j_outer*T+j_inner]*M[(j_outer*T+j_inner)*k+k_outer*T+k_inner];

            B[(gid+i_inner)*k+k_outer*T+k_inner] += sum;

          }
        }

      }
    }

    gid += stride;

  }

}

// A is nxn, M is nxk, B is nxk
__kernel
void tensor_IIA_single(__global const double *A, __global const double *M,
                       __global double *B, int n) {

  int gid = get_global_id(0);
  const int stride = get_global_size(0);

  double sum;
  while(gid<n) {

    for(int k=0; k<n; ++k) {
      for(int j=0; j<n; ++j) {

        sum = 0.0;
        for(int i=0; i<n; ++i)
          sum += A[n*gid+i]*M[k*n*n+n*j+i];

        B[k*n*n+n*j+gid] = sum;

      }

    }

    gid += stride;

  }

}

// A is nxn, M is nxk, B is nxk
__kernel
void tensor_IAI_single(__global const double *A, __global const double *M,
                       __global double *B, int n) {

  int gid = get_global_id(0);
  const int stride = get_global_size(0);

  double sum;
  while(gid<n) {

    for(int k=0; k<n; ++k) {
      for(int j=0; j<n; ++j) {

        sum = 0.0;
        for(int i=0; i<n; ++i)
          sum += A[n*gid+i]*M[k*n*n+n*i+j];

        B[k*n*n+n*gid+j] = sum;

      }

    }

    gid += stride;

  }

}

// A is nxn, M is nxk, B is nxk
__kernel
void tensor_AII_single(__global const double *A, __global const double *M,
                       __global double *B, int n) {

  int gid = get_global_id(0);
  const int stride = get_global_size(0);

  double sum;
  while(gid<n) {

    for(int k=0; k<n; ++k) {
      for(int j=0; j<n; ++j) {

        sum = 0.0;
        for(int i=0; i<n; ++i)
          sum += A[n*gid+i]*M[i*n*n+n*k+j];

        B[gid*n*n+n*k+j] = sum;

      }

    }

    gid += stride;

  }

}

__kernel
void tensor_tile(__global const double *A, __global const double *M,
                 __global double *B, int n, int tile_size, int imat) {

  // Assume n%tile_size==0

  // gid is i_outer*T
  int gid = get_global_id(0)*tile_size;
  const int stride = get_global_size(0)*tile_size;
  const int T = tile_size;

  int ni, nj, nk, nm;
  if(imat==0) {
    ni = 1;
    nj = 1;
    nk = n*n;
    nm = n;
  }
  else if(imat==1) {
    ni = n;
    nj = n;
    nk = 1;
    nm = n*n;
  }
  else if(imat==2) {
    ni = n*n;
    nj = n*n;
    nk = 1;
    nm = n;
  }

  double sum;
  while(gid<n) {

    for(int m=0; m<n; ++m) {

      for(int j_outer=0; j_outer<n/T; ++j_outer) {
        for(int k_outer=0; k_outer<n/T; ++k_outer) {
          for(int i_inner=0; i_inner<T; ++i_inner) {
            for(int k_inner=0; k_inner<T; ++k_inner) {

              sum = 0.0;
              for(int j_inner=0; j_inner<T; ++j_inner)
                sum += A[(gid+i_inner)*n+j_outer*T+j_inner]*M[m*nm+
                                                              (j_outer*T+j_inner)*nj+
                                                              (k_outer*T+k_inner)*nk];

              B[m*nm+(gid+i_inner)*ni+(k_outer*T+k_inner)*nk] += sum;

            }
          }

        }
      }

    }

    gid += stride;

  }

}

// __kernel
// void tensor_IIA_tile(__global const double *A, __global const double *M,
//                      __global double *B, int n, int tile_size) {

//   // Assume n%tile_size==0

//   // gid is i_outer*T
//   int gid = get_global_id(0)*tile_size;
//   const int stride = get_global_size(0)*tile_size;
//   const int T = tile_size;


//   double sum;
//   while(gid<n) {

//     for(int m=0; m<n; ++m) {

//       for(int j_outer=0; j_outer<n/T; ++j_outer) {
//         for(int k_outer=0; k_outer<n/T; ++k_outer) {
//           for(int i_inner=0; i_inner<T; ++i_inner) {
//             for(int k_inner=0; k_inner<T; ++k_inner) {

//               sum = 0.0;
//               for(int j_inner=0; j_inner<T; ++j_inner)
//                 sum += A[(gid+i_inner)*n+j_outer*T+j_inner]*M[m*n+
//                                                               (j_outer*T+j_inner)+
//                                                               (k_outer*T+k_inner)*n*n];

//               B[m*n+(gid+i_inner)+(k_outer*T+k_inner)*n*n] += sum;

//             }
//           }

//         }
//       }

//     }

//     gid += stride;

//   }

// }

// __kernel
// void tensor_IAI_tile(__global const double *A, __global const double *M,
//                      __global double *B, int n, int tile_size) {

//   // Assume n%tile_size==0

//   // gid is i_outer
//   int gid = get_global_id(0)*tile_size;
//   const int stride = get_global_size(0)*tile_size;
//   const int T = tile_size;

//   double sum;
//   while(gid<n) {

//     for(int m=0; m<n; ++m) {

//       for(int j_outer=0; j_outer<n/T; ++j_outer) {
//         for(int k_outer=0; k_outer<n/T; ++k_outer) {
//           for(int i_inner=0; i_inner<T; ++i_inner) {
//             for(int k_inner=0; k_inner<T; ++k_inner) {

//               sum = 0.0;
//               for(int j_inner=0; j_inner<T; ++j_inner)
//                 sum += A[(gid+i_inner)*n+j_outer*T+j_inner]*M[(j_outer*T+j_inner)*n+
//                                                               k_outer*T+k_inner+
//                                                               m*n*n];

//               B[(gid+i_inner)*n+k_outer*T+k_inner+m*n*n] += sum;

//             }
//           }

//         }
//       }

//     }

//     gid += stride;

//   }

// }

// __kernel
// void tensor_AII_tile(__global const double *A, __global const double *M,
//                      __global double *B, int n, int tile_size) {

//   // Assume n%tile_size==0

//   // gid is i_outer
//   int gid = get_global_id(0)*tile_size;
//   const int stride = get_global_size(0)*tile_size;
//   const int T = tile_size;

//   double sum;
//   while(gid<n) {

//     for(int m=0; m<n; ++m) {

//       for(int j_outer=0; j_outer<n/T; ++j_outer) {
//         for(int k_outer=0; k_outer<n/T; ++k_outer) {
//           for(int i_inner=0; i_inner<T; ++i_inner) {
//             for(int k_inner=0; k_inner<T; ++k_inner) {

//               sum = 0.0;
//               for(int j_inner=0; j_inner<T; ++j_inner)
//                 sum += A[(gid+i_inner)*n+j_outer*T+j_inner]*M[m*n+
//                                                               (j_outer*T+j_inner)*n*n+
//                                                               (k_outer*T+k_inner)];

//               B[m*n+(gid+i_inner)*n*n+(k_outer*T+k_inner)] += sum;

//             }
//           }

//         }
//       }

//     }

//     gid += stride;

//   }

// }
