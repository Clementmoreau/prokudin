//
//  gauss.c
//
//
//  Created by Clément on 28/02/2014.
//
//

// program to create the zoom-out pyramid

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


// zoom-out by a factor 2
//
// in: input image
// iw: input image width
// ih: input image height
// out: output image to be filled-in
// ow: output image width (supplied by the user)
// oh: output image height (supplied by the user)
//
static void zoom_out_by_factor_two(float *out, int ow, int oh,
                                   float *in, int iw, int ih)
{
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
        for (int i = 0; i < ow; i++)
        {
            float a[4];
            a[0] = getpixel_0(in, iw, ih, 2*i, 2*j);
            a[1] = getpixel_0(in, iw, ih, 2*i+1, 2*j);
            a[2] = getpixel_0(in, iw, ih, 2*i, 2*j+1);
            a[3] = getpixel_0(in, iw, ih, 2*i+1, 2*j+1);
            out[ow*j + i] = (a[0] + a[1] + a[2] + a[3])/4;
        }
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
	if (argc != 9) {
		fprintf(stderr, "usage:\n\t%s image out1 out2 out3 out4 out5 out6 out7 \n", *argv);
		//                          0 1     2    3    4    5    6    7    8
		return 1;
	}
	char *filename_image = argv[1];
	char *filename_zo1 = argv[2];
    char *filename_zo2 = argv[3];
    char *filename_zo3 = argv[4];
    char *filename_zo4 = argv[5];
    char *filename_zo5 = argv[6];
    char *filename_zo6 = argv[7];
    char *filename_zo7 = argv[8];
    
	// read input image
	int w, h;
	float *im = iio_read_image_float(filename_image, &w, &h);
    
    // compute size of the output images
    float w1[8];
    float h1[8];
    w1[0]=w;
    h1[0]=h;
    for (int j=1; j<8; j++)
    {
        w1[j]=ceil(w1[j-1]/2.0);
        h1[j]=ceil(h1[j-1]/2.0);
    }
    
	// allocate space for the output images
	float *out1 = malloc(w1[1]*h1[1]*sizeof(float));
    float *img = malloc(w1[0]*h1[0]*sizeof(float));
    
    float *out2 = malloc(w1[2]*h1[2]*sizeof(float));
    float *out1g = malloc(w1[1]*h1[1]*sizeof(float));

    float *out3 = malloc(w1[3]*h1[3]*sizeof(float));
    float *out2g = malloc(w1[2]*h1[2]*sizeof(float));

    float *out4 = malloc(w1[4]*h1[4]*sizeof(float));
    float *out3g = malloc(w1[3]*h1[3]*sizeof(float));

    float *out5 = malloc(w1[5]*h1[5]*sizeof(float));
    float *out4g = malloc(w1[4]*h1[4]*sizeof(float));

    float *out6 = malloc(w1[6]*h1[6]*sizeof(float));
    float *out5g = malloc(w1[5]*h1[5]*sizeof(float));

    float *out7 = malloc(w1[7]*h1[7]*sizeof(float));
    float *out6g = malloc(w1[6]*h1[6]*sizeof(float));

    
	// run the algorithm
	gfilter(img, im, w1[0], h1[0]);
    zoom_out_by_factor_two(out1, w1[1], h1[1], img, w1[0], h1[0]);
    
    gfilter(out1g, out1, w1[1], h1[1]);
    zoom_out_by_factor_two(out2, w1[2], h1[2], out1g, w1[1], h1[1]);
    
    gfilter(out2g, out2, w1[2], h1[2]);
    zoom_out_by_factor_two(out3, w1[3], h1[3], out2g, w1[2], h1[2]);
    
    gfilter(out3g, out3, w1[3], h1[3]);
    zoom_out_by_factor_two(out4, w1[4], h1[4], out3g, w1[3], h1[3]);
    
    gfilter(out4g, out4, w1[4], h1[4]);
    zoom_out_by_factor_two(out5, w1[5], h1[5], out4g, w1[4], h1[4]);
    
    gfilter(out5g, out5, w1[5], h1[5]);
    zoom_out_by_factor_two(out6, w1[6], h1[6], out5g, w1[5], h1[5]);
    
    gfilter(out6g, out6, w1[6], h1[6]);
    zoom_out_by_factor_two(out7, w1[7], h1[7], out6g, w1[6], h1[6]);
    
	// save the output images
	iio_save_image_float(filename_zo1, out1, w1[1], h1[1]);
    iio_save_image_float(filename_zo2, out2, w1[2], h1[2]);
    iio_save_image_float(filename_zo3, out3, w1[3], h1[3]);
    iio_save_image_float(filename_zo4, out4, w1[4], h1[4]);
    iio_save_image_float(filename_zo5, out5, w1[5], h1[5]);
    iio_save_image_float(filename_zo6, out6, w1[6], h1[6]);
    iio_save_image_float(filename_zo7, out7, w1[7], h1[7]);
    
	// cleanup and exit
	free(im);
	free(img);
    free(out1);
    free(out1g);
    free(out2);
    free(out2g);
    free(out3);
    free(out3g);
    free(out4);
    free(out4g);
    free(out5);
    free(out5g);
    free(out6);
    free(out6g);
    free(out7);
    
	return 0;
}

