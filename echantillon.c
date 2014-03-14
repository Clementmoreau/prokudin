//
//  echantillon.c
//  
//
//  Created by Clément on 14/03/2014.
//
//

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

// Program to get a nxn image from a larger image
// n : size of the output image
// i,j : up - left corner of the output image

void echantillon_size_n(float *in, int w, int n, int i, int j, float *out)
{
    for(int k = 0 ; k < n ; k++)
        for(int l = 0 ; l < n ; l++)
        {
            int ii=i+k;
            int jj=j+l;
            out[k*n+l]=in[ii*w+jj];
        }
}

int main(int argc, char *argv[])
{
    // process input arguments
    if (argc != 6) {
        fprintf(stderr, "usage:\n \t imageIn n i j imageOut \n");
        //                         0     1   2 3 4 5
    }
    
    char *filename_ImgIn1 = argv[1];
    int n = atof(argv[2]);
    int i = atof(argv[3]);
    int j = atof(argv[4]);
    char *filename_ImgOut = argv[5];
    
    //read input images
    int w, h;
    float *im = iio_read_image_float(filename_ImgIn1, &w, &h);
    
    //allocate space for output image
    float *out = malloc(n*n*sizeof(float));
    
    echantillon_size_n(im, w, n, i, j, out);
    
    //save image
    iio_save_image_float(filename_ImgOut, out, n, n);
    
    //cleanup and exit
    free(out);
    
    return 0;

}

