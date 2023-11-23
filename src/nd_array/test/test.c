#include "../src/nd_array.h"

#include <time.h>
#include <unistd.h>

enum { NS_PER_SECOND = 1000000000 };

void sub_timespec(struct timespec t1, struct timespec t2, struct timespec *td)
{
    td->tv_nsec = t2.tv_nsec - t1.tv_nsec;
    td->tv_sec  = t2.tv_sec - t1.tv_sec;
    if (td->tv_sec > 0 && td->tv_nsec < 0)
    {
        td->tv_nsec += NS_PER_SECOND;
        td->tv_sec--;
    }
    else if (td->tv_sec < 0 && td->tv_nsec > 0)
    {
        td->tv_nsec -= NS_PER_SECOND;
        td->tv_sec++;
    }
}



int main(void)
{
    /* test alloc.c*/ //init, uninit, malloc, calloc, free, init_tranpose, init_slice, init_strip_dims

    nd_arr_c read_arr, write_arr, read_sub_arr; 

    struct timespec start, finish, delta; // timing vars

    clock_gettime(CLOCK_REALTIME, &start);

    /*  TEST-1  */
    printf("*************** Running Test-1 ................................ \n");
    nd_init_c(&read_arr, 0, NULL);
    nd_init_c(&read_sub_arr, 0, NULL);

    nd_read_c("nc.temp", "exc_elph", &read_arr);
    nd_read_sub_c("nc.temp", "exc_elph", &read_sub_arr, ND_ALL, nd_idx{1,3,1}, nd_idx{173,356,7}, nd_idx{324,894,9} ); //nd_idx{0,ND_END,1}, nd_idx{0,ND_END,1},nd_idx{0,ND_END,1}  );
    //:,1:3,173:356:7, 324:894:9
    nd_init_c(&write_arr, *(read_arr.rank), read_arr.dims);
    
    nd_malloc_c(&write_arr);
    
    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

    clock_gettime(CLOCK_REALTIME, &start);


    nd_copy_c(&read_arr,&write_arr);

    nd_write_c("nc.temp2", "exc_elph", &write_arr, (char * [4]) {"nq", "modes", "Sf", "Si"},(size_t[]){ 1,1,16,16,2});  
    nd_write_c("nc.temp_sub", "exc_elph", &read_sub_arr, (char * [4]) {"nq", "modes", "Sf", "Si"},NULL);   
    //:, 1:2:5, 13:23:3, 17: 64:1
    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

    clock_gettime(CLOCK_REALTIME, &start);
    printf("*************** Running Test-2 ................................ \n");

    /*  TEST-2  */

    printf("Size of array read = %lld \n" ,nd_size_c(&write_arr));

    float complex element = *nd_ele_c(&write_arr, nd_idx{0,2,999,999});
    printf("%10e + %10ej \n", (double) crealf(element),(double) cimagf(element));

    element = *nd_ele_c(&write_arr, nd_idx{0,1,348,749});
    printf("%10e + %10ej \n", (double) crealf(element),(double) cimagf(element));

     /*  TEST-3  */
    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

    clock_gettime(CLOCK_REALTIME, &start);

    printf("*************** Running Test-3 ................................ \n");
    nd_arr_c reshaped_arr, transpose_arr, sliced_arr, stripped_arr;

    nd_init_c(&reshaped_arr, 5, nd_idx{6,500,2,125,4} );
    /** Reshape **/
    nd_reshape_c(&write_arr,&reshaped_arr);

    element = *nd_ele_c(&reshaped_arr, nd_idx{3,127,1,73,2});
    printf("Reshaped element : %10e + %10ej \n", (double) crealf(element),(double) cimagf(element));

    /** Transpose **/
    nd_init_tranpose_c(&reshaped_arr, nd_idx{4,2,0,1,3}, &transpose_arr); // 4,2,6,500,125
    nd_malloc_c(&transpose_arr);
    nd_tranpose_c(&reshaped_arr, nd_idx{4,2,0,1,3}, &transpose_arr);

    element = *nd_ele_c(&transpose_arr, nd_idx{2,1,3,273,39});
    printf("Transposed element : %10e + %10ej \n", (double) crealf(element),(double) cimagf(element));

    /** Slice **/
    nd_init_slice_c(nd_idx{2,0,3,179,53}, nd_idx{4,1,5,500,117}, nd_idx{1,1,1,3,5}, &transpose_arr, &sliced_arr );
    nd_malloc_c(&sliced_arr );
    
    nd_slice_c(nd_idx{2,0,3,179,53}, nd_idx{4,1,5,500,117}, nd_idx{1,1,1,3,5}, &transpose_arr, &sliced_arr);
    
    element = *nd_ele_c(&sliced_arr, nd_idx{1,0,37,7});
    printf("Sliced element : %10e + %10ej \n", (double) crealf(element),(double) cimagf(element));

    /** Strip **/
    nd_init_strip_dims_c(&sliced_arr, (ND_int)2, &stripped_arr /*nd_idx{1,0}*/);
    nd_strip_dims_c(&sliced_arr, (ND_int)2, nd_idx{1,0}, &stripped_arr);

    element = *nd_ele_c(&stripped_arr, nd_idx{13,5});
    printf("Stripped element : %10e + %10ej \n", (double) crealf(element),(double) cimagf(element));
    //
    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

    clock_gettime(CLOCK_REALTIME, &start);
 
     /*  TEST-4  */
    //clock_gettime(CLOCK_REALTIME, &start);
    printf("*************** Running Test-4 ................................ \n");
    nd_arr_c mat_prod, einsum_test, sum_test, einsum_test2;

    nd_init_c(&mat_prod,(ND_int)2, write_arr.dims+2 );
    nd_malloc_c(&mat_prod);

    nd_set_all_c(&mat_prod, 0.0f + 0.0f*I);
    nd_matmul_c('T', 'C', &write_arr, &write_arr, &mat_prod, 1.0f + 0.0f*I, 0.0f + 0.0f*I, nd_idx{0,1},  nd_idx{0,2}, NULL /*nd_idx{},  nd_idx{}*/);

    element = *nd_ele_c(&mat_prod, nd_idx{239,739});
    printf("Matmul element : %10e + %10ej \n", (double) crealf(element),(double) cimagf(element));

    nd_set_all_c(&mat_prod, 0.0f + 0.0f*I);

    nd_matmulX_c ('T', 'C', nd_ele_c(&write_arr, nd_idx{0,1,0,0}), nd_ele_c(&write_arr, nd_idx{0,2,0,0}), mat_prod.data,
                1.0f , 0.0f, write_arr.dims[3], write_arr.dims[3], mat_prod.dims[1],mat_prod.dims[0], mat_prod.dims[1], write_arr.dims[2]);

    float complex element2 = *nd_ele_c(&mat_prod, nd_idx{239,739});

    if (element == element2) printf("MatmulX test passed \n");
    else printf("MatmulX test failed \n");
    
    nd_init_c(&einsum_test,*(write_arr.rank)-1,write_arr.dims);
    nd_malloc_c(&einsum_test);

    nd_set_all_c(&einsum_test, 0.0f + 0.0f*I);
    //nd_einsum_c("ijkk,ijlk->ijl",&write_arr,&write_arr,&einsum_test,1.0f + 0.0f*I, 0.0f + 0.0f*I);

    //element = *nd_ele_c(&einsum_test, nd_idx{0,2,473});
    //printf("Einsum element : %10e + %10ej \n", (double) crealf(element),(double) cimagf(element));

    nd_init_c(&einsum_test2,0,NULL);
    nd_init_c(&sum_test,0,NULL);
    
    nd_malloc_c(&sum_test);
    nd_malloc_c(&einsum_test2);

    nd_set_all_c(&einsum_test2, 0.0f + 0.0f*I);
    nd_set_all_c(&sum_test, 0.0f + 0.0f*I);

    //nd_einsum_c("ijkx,ijky->",&write_arr,&write_arr,&einsum_test2,1.0f + 0.0f*I, 0.0f + 0.0f*I);
    //nd_sum_c("ijl","",&einsum_test,&sum_test,1.0f + 0.0f*I, 0.0f + 0.0f*I);

    //element = *nd_ele_c(&einsum_test2, nd_idx{0});
    //printf("Einsum element 2 : %10e + %10ej \n", (double) crealf(element),(double) cimagf(element));
    //
    //element = *nd_ele_c(&sum_test, nd_idx{0});
    //printf("Sum element : %10e + %10ej \n", (double) crealf(element),(double) cimagf(element));

    nd_free_c(&sum_test);
    nd_free_c(&einsum_test2);
    nd_free_c(&mat_prod);
    nd_free_c(&einsum_test);
    nd_free_c(&read_arr);
    nd_free_c(&write_arr);
    nd_free_c(&transpose_arr);
    nd_free_c(&sliced_arr);
    nd_free_c(&read_sub_arr);

    nd_uninit_c(&read_sub_arr);
    nd_uninit_c(&sum_test);
    nd_uninit_c(&einsum_test2);
    nd_uninit_c(&mat_prod);
    nd_uninit_c(&einsum_test);
    nd_uninit_c(&reshaped_arr);
    nd_uninit_c(&transpose_arr);
    nd_uninit_c(&sliced_arr);
    nd_uninit_c(&stripped_arr);
    nd_uninit_c(&read_arr);
    nd_uninit_c(&write_arr);

    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);
    
    printf(" ************** All Test Sucessfull :) *************** \n");

}









// int main(void){

//     nd_arr_c temp_array, einsum_test;
//     nd_arr_c mat_prod, sliced_arr;
//     //nd_arr_c temp_array_T, mat_prod;
//     ND_int pppp = 1000;

//     nd_malloc_c(&einsum_test,(ND_int) 3, nd_idx{1,3,pppp});
//     //nd_malloc_c(&einsum_test,(ND_int) 0, NULL);

//     nd_malloc_c(&mat_prod, (ND_int) 2, nd_idx{pppp,pppp});

//     nd_read_c("/Users/murali/phd/one_phonon_raman/si/bse/si_data/nscf_wo/raman/ndb.BS_elph", "exc_elph", &temp_array);

//     nd_write_c("nc.temp", "elph_ex", &temp_array, (char * [4]) {"nq", "modes", "Sf", "Si"});

//     // nd_malloc_c(&sliced_arr, (ND_int) 3, nd_idx{2,150,42} );

//     // nd_slice_c(nd_idx{0,1,400,200}, nd_idx{1,3,700,700}, nd_idx{1,1,2,12}, &temp_array, &sliced_arr );

//     // double complex pp1 = (double complex) nd_ele_c(&sliced_arr, nd_idx {1,127,37})[0] ; 

//     nd_malloc_c(&sliced_arr, (ND_int) 3, nd_idx{2,75,5} );

//     nd_slice_c(nd_idx{0,1,400,200}, nd_idx{1,3,700,205}, nd_idx{1,1,4,1}, &temp_array, &sliced_arr );

//     double complex pp1 = (double complex) nd_ele_c(&sliced_arr, nd_idx {1,57,3})[0] ; 

//     printf("%10e + %10e I\n", creal(pp1), cimag(pp1) );

//     struct timespec start, finish, delta; // timing vars
    
//     clock_gettime(CLOCK_REALTIME, &start); // timing start

//     nd_set_all_c(&mat_prod, 0.0f + 0.0f*I);
//     nd_matmul_c('T', 'C', &temp_array, &temp_array, &mat_prod, 1.0f + 0.0f*I, 0.0f + 0.0f*I, nd_idx{0,1},  nd_idx{0,2}, NULL /*nd_idx{},  nd_idx{}*/);

//     nd_einsum_c("ijkk,ijlk->ijl",&temp_array,&temp_array,&einsum_test,1.0f + 0.0f*I, 0.0f + 0.0f*I);

//     double complex pp = (double complex) nd_ele_c(&mat_prod, nd_idx {50,500})[0] ; 
//     printf("%10e + %10e I\n", creal(pp), cimag(pp) );

//     clock_gettime(CLOCK_REALTIME, &finish); //timing end:
//     sub_timespec(start, finish, &delta); // timing end:
//     printf("Time : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

//     for (ND_int i=0; i<1; ++i){
//         pp = (double complex) nd_ele_c(&einsum_test, nd_idx {0,0,i})[0] ; 
//         printf("%10e + %10e I\n", creal(pp), cimag(pp) );
//     }

//     nd_free_c(&sliced_arr);
//     nd_free_c(&mat_prod);
//     nd_free_c(&einsum_test);
//     // nd_free_c(&temp_array_T);
//     nd_free_c(&temp_array);


//     return 0;

// }

