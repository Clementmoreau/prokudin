//
//  homography_combination.c
//  
//
//  Created by Clément on 02/06/2014.
//
//

#include <stdio.h>
#include <math.h>
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

static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
                           + x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
                                + x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

typedef float (*getpixel_operator)(float*,int,int,int,int);

float bicubic_interpolation(float *img, int w, int h, float x, float y)
{
	x -= 1;
	y -= 1;
    
	getpixel_operator p = getpixel_0;
    
	int ix = floor(x);
	int iy = floor(y);
	float c[4][4];
	for (int j = 0; j < 4; j++)
    for (int i = 0; i < 4; i++)
    c[i][j] = p(img, w, h, ix + i, iy + j);
	return bicubic_interpolation_cell(c, x - ix, y - iy);
}

// compute the distance between two pixels in an image

float dist(int x1, int y1, int x2, int y2)
{
    float distx = abs(x2-x1)*abs(x2-x1);
    float disty = abs(y2-y1)*abs(y2-y1);
    return sqrt(distx + disty);
}

/* function to get n^2 subimages from the initial image

void cut_n_parts(float *im, int w, int h, int n, float *coord, float **out)
{
    int w1 = floor(2*w/(n+1));
    int h1 = floor(2*h/(n+1));
    
    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++)
            for (int k = 0 ; k < w1 ; k++)
                for (int l = 0 ; l < h1 ; l++)
                {
                    out[j*n+i][l*w1+k] = getpixel_0(im, w, h, j*w1/2+k, i*h1/2+l);
                    coord = {coord ; (j+1)*w1/2 ; (i+1)*h1/2};
                }
}
*/

// get the coordinates of the centers of the nxn subimages
void coordinates_centers(int w, int h, int n, float *coordx, float *coordy)
{
    int w1 = floor(w/(n+1));
    int h1 = floor(h/(n+1));
    
    for (int k = 0 ; k < n ; k++)
        for (int l = 0 ; l < n ; l++)
        {
            coordx[k*n+l] = (k+1)*w1;
            coordy[k*n+l] = (l+1)*h1;
        }
    
}

//get the inverse of a 3x3 homography matrix

void inv_homography(float h[9], float inv[9])
{
    float det = h[0]*h[4]*h[8]+h[1]*h[5]*h[6]+h[2]*h[3]*h[7]-h[2]*h[4]*h[6]-h[1]*h[3]*h[8]-h[0]*h[5]*h[7];
    inv[0] = (h[4]*h[8]-h[5]*h[7])/det;
    inv[1] = (h[2]*h[7]-h[1]*h[8])/det;
    inv[2] = (h[1]*h[5]-h[2]*h[4])/det;
    inv[3] = (h[5]*h[6]-h[3]*h[8])/det;
    inv[4] = (h[0]*h[8]-h[2]*h[6])/det;
    inv[5] = (h[2]*h[3]-h[0]*h[5])/det;
    inv[6] = (h[3]*h[7]-h[4]*h[6])/det;
    inv[7] = (h[1]*h[6]-h[0]*h[7])/det;
    inv[8] = (h[0]*h[4]-h[1]*h[3])/det;
}

//get the image by homography for a destination pixel x, y
float destination_pixel(float inv[9], float *im, int w, int h, int x, int y)
{
    float p[3] = {x, y, 1};
    float p1[3];
    p1[0] = p[0]*inv[0] + p[1]*inv[1] + p[2]*inv[2];
    p1[1] = p[0]*inv[3] + p[1]*inv[4] + p[2]*inv[5];
    p1[2] = p[0]*inv[6] + p[1]*inv[7] + p[2]*inv[8];
    
    return bicubic_interpolation(im, w, h, p1[0]/p1[2], p1[1]/p1[2]);
}

void destination_position(float inv[9], int x, int y, float out[2])
{
    float p[3] = {x, y, 1};
    float p1[3];
    p1[0] = p[0]*inv[0] + p[1]*inv[1] + p[2]*inv[2];
    p1[1] = p[0]*inv[3] + p[1]*inv[4] + p[2]*inv[5];
    p1[2] = p[0]*inv[6] + p[1]*inv[7] + p[2]*inv[8];
    
    out[0] = p1[0]/p1[2];
    out[1] = p1[1]/p1[2];
}


// apply the transformation

// nxn : number of tiles, thus number of homographies
// coord : coordinates of the centers of the tiles
// w, h : size of the image
// im : nxn warped images
// out : registered image

/*
void apply_transformation2(int n, float *coordx, float *coordy, float *im, float *hom, int w, int h, float *out)
{
    int delta = (w+h)/(n+2) ;
    int nn = n*n ;
    
    float inv[nn][9];
    
    for(int k = 0 ; k < nn ; k++)
    {
        inv_homography(hom+9*k, inv[k]);
    }
    
    for(int i = 0 ; i < w ; i++)
        for(int j = 0 ; j < h ; j++)
        {
            float S = 0;
            float S2 = 0;
            for(int k = 0 ; k < nn ; k++)
            {
                float p = destination_pixel(inv[k], im, w, h, i, j);
                S = S + p*exp(-dist(i, j, coordx[k], coordy[k])/delta);
                S2 = S2 + exp(-dist(i, j, coordx[k], coordy[k])/delta);
            }
            out[j*w+i] = S/S2 ;
        }
}
 */

void compute_transformation(int n, float *coordx, float *coordy, float *hom, float *outx, float *outy, int w, int h, float delta)
{
    int nn = n*n ;

    for(int i = 0 ; i < w ; i++)
        for(int j = 0 ; j < h ; j++)
        {
            float Sx = 0;
            float Sy = 0;
            float S2 = 0;
            for(int k = 0 ; k < nn ; k++)
            {
                float p[2];
                destination_position(hom+9*k, i, j, p);
                float W = exp(-dist(i, j, coordx[k], coordy[k])/delta);
                Sx = Sx + p[0]*W;
                Sy = Sy + p[1]*W;
                S2 = S2 + W;
            }
            outx[j*w+i] = Sx/S2 ;
            outy[j*w+i] = Sy/S2 ;
        }

}

void apply_transformation(float *im, float *Tx, float *Ty, int w, int h, float *out)
{
    for(int i = 0 ; i < w ; i++)
        for(int j = 0 ; j < h ; j++)
        {
            float dx = Tx[j*w+i];
            float dy = Ty[j*w+i];
            out[j*w+i] = bicubic_interpolation(im, w, h, dx, dy);
        }
    
}

int main(int argc, char **argv)
{
    // process input arguments
    if (argc != 6)
    {
        fprintf(stderr, "usage: \n \t imageIn n hom imageOut parameter \n");
        return 0;
    }
    
    char *filename_ImgIn = argv[1];
    int n = atoi(argv[2]);
    char *hom = argv[3];
    char *filename_ImgOut = argv[4];
    float delta = atof(argv[5]);
    
    char filename_hom[n*n][256];
    int nn = n*n;
    int w, h;
    float *im = iio_read_image_float(filename_ImgIn, &w, &h);
    
    float H[nn*9];
    
    float *coordx = malloc(n*n*sizeof(float));
    float *coordy = malloc(n*n*sizeof(float));
    coordinates_centers(w, h, n, coordx, coordy);
    
    
    for(int i = 0 ; i < nn ; i++)
    {
        sprintf(filename_hom[i], "%s_%d.txt", hom, i+1);
        
        FILE* f = NULL;
        f = fopen(filename_hom[i], "r");
        
        if (f != NULL)
        {
            fscanf(f, "%f %f %f %f %f %f %f %f %f", &H[9*i+0], &H[9*i+1], &H[9*i+2], &H[9*i+3], &H[9*i+4], &H[9*i+5], &H[9*i+6], &H[9*i+7], &H[9*i+8]);
            fclose(f);
        }
        else
        {
            printf("Impossible to open file");
        }
        
    }
    
    float *out = malloc(w*h*sizeof(float));
    float *Tx = malloc(w*h*sizeof(float));
    float *Ty = malloc(w*h*sizeof(float));
    delta = delta*(w+h)/(n+2);
    
    compute_transformation(n, coordx, coordy, H, Tx, Ty, w, h, delta);
    apply_transformation(im, Tx, Ty, w, h, out);
    
    iio_save_image_float(filename_ImgOut, out, w, h);
    
    free(out);
    free(Tx);
    free(Ty);
    free(coordx);
    free(coordy);
    
    return 0;
    
}











