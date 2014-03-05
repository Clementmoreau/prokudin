// fonction pour calculer les derivees directionnelles

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


// apply a translation to the given image
//
// in: input image
// w: width
// h: height
// dx: horizontal displacement
// dy: vertical displacement
// out: output image, to be filled-in
//
void apply_translation(float *out, int dx, int dy, float *in, int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int ii = i - dx;
		int jj = j - dy;
		out[j*w+i] = getpixel_0(in, w, h, ii, jj);
	}
}

// main function
int main(int argc, char **argv)
{
	// process input arguments
	if (argc != 6) {
		fprintf(stderr, "usage:\n\t%s image der1 der2 der3 der4\n", *argv);
		//                          0 1     2     3    4    5
		return 1;
	}
	char *filename_image = argv[1];
	char *filename_der1 = argv[2];
	char *filename_der2 = argv[3];
	char *filename_der3 = argv[4];
	char *filename_der4 = argv[5];

	// read input image
	int w, h;
	float *im = iio_read_image_float(filename_image, &w, &h);

	// allocate space for the output image
	float *out1 = malloc(w*h*sizeof(float));
	float *out2 = malloc(w*h*sizeof(float));
	float *out3 = malloc(w*h*sizeof(float));
	float *out4 = malloc(w*h*sizeof(float));
	float *dec1 = malloc(w*h*sizeof(float));
	float *dec2 = malloc(w*h*sizeof(float));
	float *dec3 = malloc(w*h*sizeof(float));
	float *dec4 = malloc(w*h*sizeof(float));

	// creer les images translatees
	apply_translation(dec1, 1, 0, im, w, h);
	apply_translation(dec2, 0, 1, im, w, h);
	apply_translation(dec3, -1, 1, im, w, h);
	apply_translation(dec4, -1, -1, im, w, h);

	// calculer les derivees
	for(int j = 0 ; j < h ; j++)
    for(int i = 0 ; i < w ; i++)
	{
	    out1[j*w+i] = (dec1[j*w+i] - im[j*w+i])*(dec1[j*w+i] - im[j*w+i]);
        out2[j*w+i] = (dec2[j*w+i] - im[j*w+i])*(dec2[j*w+i] - im[j*w+i]);
        out3[j*w+i] = (dec3[j*w+i] - im[j*w+i])*(dec3[j*w+i] - im[j*w+i])/2;
        out4[j*w+i] = (dec4[j*w+i] - im[j*w+i])*(dec4[j*w+i] - im[j*w+i])/2;
	}


	// save the output image
	iio_save_image_float(filename_der1, out1, w, h);
	iio_save_image_float(filename_der2, out2, w, h);
	iio_save_image_float(filename_der3, out3, w, h);
	iio_save_image_float(filename_der4, out4, w, h);

	// cleanup and exit
	free(im);
	free(out1);
	free(out2);
	free(out3);
	free(out4);
	free(dec1);
	free(dec2);
	free(dec3);
	free(dec4);
	return 0;
}
