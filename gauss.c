//
//  gauss.c
//  
//
//  Created by Clément on 28/02/2014.
//
//

// program to apply a gaussian filter 5x5
// method : use discrete values

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"


// auxiliary function to get the value of an image at any point (i,j)
// (points outisde the original domain get the value 0)
//
// x: image data
// w: width
// h: height
// i: horizontal position
// j: vertical position
//
// return value: color of the requested pixel
//
float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return 0;
	else
		return x[j*w+i];
}


// apply the gaussian filter
//
// in: input image
// w: width
// h: height
// out: output image, to be filled-in
//
void gfilter(float *out, float *in, int w, int h)
{
	for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            out[j*w+i] = floor((359*getpixel_0(in, w, h, i, j)+164*(getpixel_0(in, w, h, i-1, j)+getpixel_0(in, w, h, i, j-1)+getpixel_0(in, w, h, i+1, j)+getpixel_0(in, w, h, i, j+1))+75*(getpixel_0(in, w, h, i+1, j+1)+getpixel_0(in, w, h, i-1, j-1)+getpixel_0(in, w, h, i-1, j+1)+getpixel_0(in, w, h, i+1, j-1))+16*(getpixel_0(in, w, h, i+2, j)+getpixel_0(in, w, h, i, j+2)+getpixel_0(in, w, h, i-2, j)+getpixel_0(in, w, h, i, j-2))+7*(getpixel_0(in, w, h, i-2, j+1)+getpixel_0(in, w, h, i-2, j-1)+getpixel_0(in, w, h, i-1, j+2)+getpixel_0(in, w, h, i-1, j-2)+getpixel_0(in, w, h, i+1, j+2)+getpixel_0(in, w, h, i+1, j-2)+getpixel_0(in, w, h, i+2, j+1)+getpixel_0(in, w, h, i+2, j-1))+getpixel_0(in, w, h, i+2, j+2)+getpixel_0(in, w, h, i+2, j-2)+getpixel_0(in, w, h, i-2, j+2)+getpixel_0(in, w, h, i-2, j-2))/1444);
        }
}


// main function
int main(int argc, char **argv)
{
	// process input arguments
	if (argc != 3) {
		fprintf(stderr, "usage:\n\t%s image filtered_image \n", *argv);
		//                          0 1     2
		return 1;
	}
	char *filename_image = argv[1];
	char *filename_filtered_image = argv[2];
    
	// read input image
	int w, h;
	float *im = iio_read_image_float(filename_image, &w, &h);
    
	// allocate space for the output image
	float *out = malloc(w*h*sizeof(float));
    
	// run the algorithm
	gfilter(out, im, w, h);
    
	// save the output image
	iio_save_image_float(filename_filtered_image, out, w, h);
    
	// cleanup and exit
	free(im);
	free(out);
	return 0;
}

