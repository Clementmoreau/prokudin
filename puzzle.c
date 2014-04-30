//
//  puzzle.c
//  
//
//  Created by Clément on 30/04/2014.
//
//

#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iio.h"

static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return 0;
	else
		return x[j*w+i];
}

void puzzle(float **im, int w1, int h1, int n, float *out)
{
    int w = w1*n;
    int h = h1*n;
    
    for (int i = 0 ; i < n ; i++)
        for (int j = 0 ; j < n ; j++)
            for (int k = 0 ; k < w1 ; k++)
                for (int l = 0 ; l < h1 ; l++)
                {
                    out[(i*h1+l)*w+j*w1+k] = getpixel_0(im[j*n+i], w1, h1, k, l);
                }
    
    //other way to fill the image
    
    /*
     for (int i = 0 ; i < w ; i++)
     for (int j = 0 ; j < h ; j++)
     {
     int k=floor(j/w1);
     int l=floor(i/h1);
     out[j*w+i]=getpixel_0(im[k*n+l], w1, h1, i-k*w1, j-l*h1);
     }
     */
}

int main(int argc, char **argv)
{
    //process input arguments
    if (argc != 4)
    {
        fprintf(stderr, "usage:\n \t tile n image Out");
        return 0;
    }
    
    char *filename_tile = argv[1];
    int n = atoi(argv[2]);
    char *filename_imgOut = argv[3];
    
    //read input images
    int w1, h1;
    char filenames_tiles[n*n][256];
    
    for(int i = 0 ; i < n*n ; i++)
    {
        sprintf(filenames_tiles[i], "%s_%d.png", filename_tile, i+1);
    }
    
    float *tile[n*n];
    
    for (int i = 0 ; i < n*n ; i++)
    {
        tile[i] = iio_read_image_float(filenames_tiles[i], &w1, &h1);
    }
    
    //allocate space for rebuilt image
    
    int w2 = w1*n;
    int h2 = h1*n;
    
    float *out = malloc(w2*h2*sizeof(float));
    
    //rebuild image
    puzzle(tile, w1, h1, n, out);
    
    //save outputs
    iio_save_image_float(filename_imgOut, out, w2, h2);
    
    //cleanup and exit
    free(out);
    
    return 0;
    
}



    
















