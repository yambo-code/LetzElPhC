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

    nd_arr_s read_arr, write_arr, read_sub_arr; 

    struct timespec start, finish, delta; // timing vars

    clock_gettime(CLOCK_REALTIME, &start);

    /*  TEST-1  */
    printf("*************** Running Test-1 ................................ \n");
    nd_init_s(&read_arr, 0, NULL);
    nd_init_s(&read_sub_arr, 0, NULL);

    nd_read_s("nc.temp", "exc_elph", &read_arr);
    nd_read_sub_s("nc.temp", "exc_elph", &read_sub_arr, ND_ALL, nd_idx{1,3,1}, nd_idx{173,356,7}, nd_idx{324,894,9}, nd_idx{0,ND_END,1} ); //nd_idx{0,ND_END,1}, nd_idx{0,ND_END,1},nd_idx{0,ND_END,1}  );
    //:,1:3,173:356:7, 324:894:9
    nd_init_s(&write_arr, *(read_arr.rank), read_arr.dims);
    
    nd_malloc_s(&write_arr);
    
    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

    clock_gettime(CLOCK_REALTIME, &start);


    nd_copy_s(&read_arr,&write_arr);

    nd_write_s("nc.temp_2s", "exc_elph", &write_arr, (char * []) {"nq", "modes", "Sf", "Si", "re_im"},(size_t[]){ 1,1,512,512,2});  
    nd_write_s("nc.temp_2ssub", "exc_elph", &read_sub_arr, (char * []) {"nq", "modes", "Sf", "Si", "re_im"},NULL);   
    //:, 1:2:5, 13:23:3, 17: 64:1
    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

    clock_gettime(CLOCK_REALTIME, &start);
    printf("*************** Running Test-2 ................................ \n");

    /*  TEST-2  */

    printf("Size of array read = %lld \n" ,nd_size_s(&write_arr));

    float complex element = *nd_ele_s(&write_arr, nd_idx{0,2,999,999,0});
    printf("%10e \n", (double) element);

    element = *nd_ele_s(&write_arr, nd_idx{0,1,348,749,1});
    printf("%10e \n", (double) element);

     /*  TEST-3  */
    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

    clock_gettime(CLOCK_REALTIME, &start);

    printf("*************** Running Test-3 ................................ \n");
    nd_arr_s reshaped_arr, transpose_arr, sliced_arr, stripped_arr;

    nd_init_s(&reshaped_arr, 5, nd_idx{6,500,4,125,4} );
    /** Reshape **/
    nd_reshape_s(&write_arr,&reshaped_arr);

    element = *nd_ele_s(&reshaped_arr, nd_idx{3,127,1,73,2});
    printf("Reshaped element : %10e \n", (double) element);

    /** Transpose **/
    nd_init_tranpose_s(&reshaped_arr, nd_idx{4,2,0,1,3}, &transpose_arr); // 4,2,6,500,125
    nd_malloc_s(&transpose_arr);
    nd_tranpose_s(&reshaped_arr, nd_idx{4,2,0,1,3}, &transpose_arr);

    element = *nd_ele_s(&transpose_arr, nd_idx{2,1,3,273,39});
    printf("Transposed element : %10e \n", (double) element);

    /** Slice **/
    nd_init_slice_s(nd_idx{2,0,3,179,53}, nd_idx{4,1,5,500,117}, nd_idx{1,1,1,3,5}, &transpose_arr, &sliced_arr );
    nd_malloc_s(&sliced_arr );
    
    nd_slice_s(nd_idx{2,0,3,179,53}, nd_idx{4,1,5,500,117}, nd_idx{1,1,1,3,5}, &transpose_arr, &sliced_arr);
    
    element = *nd_ele_s(&sliced_arr, nd_idx{1,0,37,7});
    printf("Sliced element : %10e \n", (double) element);

    /** Strip **/
    nd_init_strip_dims_s(&sliced_arr, (ND_int)2, &stripped_arr /*nd_idx{1,0}*/);
    nd_strip_dims_s(&sliced_arr, (ND_int)2, nd_idx{1,0}, &stripped_arr);

    element = *nd_ele_s(&stripped_arr, nd_idx{13,5});
    printf("Stripped element : %10e \n", (double) element);
    //
    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);

    clock_gettime(CLOCK_REALTIME, &start);
 
     /*  TEST-4  */
    //clock_gettime(CLOCK_REALTIME, &start);
    printf("*************** Running Test-4 ................................ \n");
    nd_arr_s mat_prod, einsum_test, sum_test, einsum_test2;

    nd_init_s(&mat_prod,(ND_int)2, write_arr.dims+2 );
    nd_malloc_s(&mat_prod);

    nd_set_all_s(&mat_prod, 0.0f );
    nd_matmul_s('N', 'T', &write_arr, &write_arr, &mat_prod, 1.0f , 0.0f , nd_idx{0,1,300},  nd_idx{0,2,579}, NULL /*nd_idx{},  nd_idx{}*/);

    element = *nd_ele_s(&mat_prod, nd_idx{239,739});
    printf("Matmul element : %10e \n", (double) element);

    float element2 = 0.0f;

    nd_set_all_s(&mat_prod, 0.0f );
    nd_matmulX_s ('N', 'T', nd_ele_s(&write_arr, nd_idx{0,1,300,0,0}), nd_ele_s(&write_arr, nd_idx{0,2,579,0,0}), mat_prod.data,
                1.0f , 0.0f, write_arr.dims[4], write_arr.dims[4], mat_prod.dims[1],mat_prod.dims[0], mat_prod.dims[1], write_arr.dims[4]);

    element2 = *nd_ele_s(&mat_prod, nd_idx{239,739});
    if (element == element2) printf("MatmulX test passed \n");

    nd_init_s(&einsum_test,*(write_arr.rank)-2,write_arr.dims);
    nd_malloc_s(&einsum_test);

    nd_set_all_s(&einsum_test, 0.0f );
    //nd_einsum_s("ijkkp,ijlkp->ijl",&write_arr,&write_arr,&einsum_test,1.0f , 0.0f );

    //element = *nd_ele_s(&einsum_test, nd_idx{0,2,473});
    //printf("Einsum element : %10e \n", (double) element);

    nd_init_s(&einsum_test2,0,NULL);
    nd_init_s(&sum_test,0,NULL);
    
    nd_malloc_s(&sum_test);
    nd_malloc_s(&einsum_test2);

    nd_set_all_s(&einsum_test2, 0.0f );
    nd_set_all_s(&sum_test, 0.0f );

    //nd_einsum_s("ijkkp,ijlkp->",&write_arr,&write_arr,&einsum_test2,1.0f , 0.0f );
    //nd_sum_s("ijl","",&einsum_test,&sum_test,1.0f, 0.0f);

    //element = *nd_ele_s(&einsum_test2, nd_idx{0});
    //printf("Einsum element 2 : %10e \n", (double) element);
    //
    //element = *nd_ele_s(&sum_test, nd_idx{0});
    //printf("Sum element : %10e \n", (double) element);

    nd_free_s(&sum_test);
    nd_free_s(&einsum_test2);
    nd_free_s(&mat_prod);
    nd_free_s(&einsum_test);
    nd_free_s(&read_arr);
    nd_free_s(&write_arr);
    nd_free_s(&transpose_arr);
    nd_free_s(&sliced_arr);
    nd_free_s(&read_sub_arr);

    nd_uninit_s(&read_sub_arr);
    nd_uninit_s(&sum_test);
    nd_uninit_s(&einsum_test2);
    nd_uninit_s(&mat_prod);
    nd_uninit_s(&einsum_test);
    nd_uninit_s(&reshaped_arr);
    nd_uninit_s(&transpose_arr);
    nd_uninit_s(&sliced_arr);
    nd_uninit_s(&stripped_arr);
    nd_uninit_s(&read_arr);
    nd_uninit_s(&write_arr);

    clock_gettime(CLOCK_REALTIME, &finish); //timing end:
    sub_timespec(start, finish, &delta); // timing end:
    printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec, delta.tv_nsec);
    
    printf(" ************** All Test Sucessfull :) *************** \n");

}







