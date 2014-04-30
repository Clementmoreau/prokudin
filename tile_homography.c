//
//  tile_homography.c
//  
//
//  Created by Clément on 17/04/2014.
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

// function to get n^2 subimages from the initial image

void cut_n_parts(float *im, int w, int h, int n, float **out)
{
    int w1 = floor(w/n);
    int h1 = floor(h/n);
    
    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++)
            for (int k = 0 ; k < w1 ; k++)
                for (int l = 0 ; l < h1 ; l++)
                {
                    out[j*n+i][l*w1+k] = getpixel_0(im, w, h, j*w1+k, i*h1+l);
                }
}

//function to get the initial image back from the n^2 tiles



int main(int argc, char **argv)
{
    //process input arguments
    if (argc != 4)
    {
        fprintf(stderr, "usage:\n \t imageIn n tile \n");
        //                         0     1   2    3
        return 0;
    }
    
    char *filename_ImgIn = argv[1];
    int n = atoi(argv[2]);
    char *fileout = argv[3];
    
    char filename_tile[n*n][256];
    
    for(int i = 0 ; i < n*n ; i++)
    {
        sprintf(filename_tile[i], "%s_%d.png", fileout, i+1);
    }
    
    
    //read input images
    int w, h;
    float *im = iio_read_image_float(filename_ImgIn, &w, &h);
    
    int w1 = floor(w/n);
    int h1 = floor(h/n);
    
    //allocate space for subimages
    float *tiles[n*n];
    for(int i = 0 ; i < n*n ; i++)
    {
        tiles[i] = malloc(w1*h1*sizeof(float));
    }
    
    //create tiles and rebuild image
    cut_n_parts(im, w, h, n, tiles);
    
    //save outputs
    
    for (int i = 0 ; i < n*n ; i++)
    {
        iio_save_image_float(filename_tile[i], tiles[i], w1, h1);
    }
    
    //cleanup and exit
    for(int i = 0 ; i < n*n ; i++)
    {
        free(tiles[i]);
    }
    
    return 0;

}









