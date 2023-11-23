#include "../src/nd_array.h"

#include <time.h>
#include <unistd.h>

#define N  10
#define M  12
#define P  14
#define DOF  3

int main(void)
{
    /* test alloc.c*/ //init, uninit, malloc, calloc, free, init_tranpose, init_slice, init_strip_zims

    nd_arr_z in, out; 

    nd_init_z(&in, 4, nd_idx{N,M,P,DOF});
    nd_init_z(&out, 4, nd_idx{N,M,P,DOF});

    nd_malloc_z(&in);
    nd_malloc_z(&out);

    unsigned int i,j,k;

    for(i=0;i<N;i++){
        for(j=0;j<M;j++){
            for(k=0;k<P;k++){
                //printf("ijk\n");
                *nd_ele_z(&in,nd_idx{i,j,k,0}) = 30.+12.*sin(2*3.1415926535*i/((double)N))*sin(2*3.1415926535*j/((double)M))*sin(2*3.1415926535*k/((double)P))*I;
                *nd_ele_z(&in,nd_idx{i,j,k,1}) = 42.0f;
                *nd_ele_z(&in,nd_idx{i,j,k,2}) = 1.0f;
            }
        }
    }

    /*Perform the fft*/
    nd_fft_z(&in,&out, 3, nd_idx{0,1,2}, -1, true);

    printf("Test : %10e %10e \n", creal(*nd_ele_z(&in,nd_idx{7,4,11,0}) ), cimag(*nd_ele_z(&in,nd_idx{2,9,5,2}) )*10000 );
    printf("Test : %10e %10e \n", creal(*nd_ele_z(&out,nd_idx{7,4,11,0}) )*10000, cimag(*nd_ele_z(&out,nd_idx{2,9,5,2}) )*10000 );

    nd_set_all_z(&in,0.0);

    nd_fft_z(&out,&in, 3, nd_idx{0,1,2}, +1, false);

    printf("Test : %10e %10e \n", creal(*nd_ele_z(&in,nd_idx{7,4,11,0}) ), cimag(*nd_ele_z(&in,nd_idx{2,9,5,2}) )*10000 );
    printf("Test : %10e %10e \n", creal(*nd_ele_z(&out,nd_idx{7,4,11,0}) )*10000, cimag(*nd_ele_z(&out,nd_idx{2,9,5,2}) )*10000 );

    nd_destroy_z(&in);
    nd_destroy_z(&out);

    return 0;
}