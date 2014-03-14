//
//  irani.c
//  
//
//  Created by Clément on 07/03/2014.
//
//

// Program to register three images with the parametric method


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
static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return 0;
	else
		return x[j*w+i];
}

// subprogram to solve Ax=B

static int solvps(float *a,float *b,int n)
{ float *p,*q,*r,*s,t;
    int j,k;
    for(j=0,p=a; j<n ;++j,p+=n+1){
        for(q=a+j*n; q<p ;++q) *p-= *q* *q;
        if(*p<=0.) return -1;
        *p=sqrt(*p);
        for(k=j+1,q=p+n; k<n ;++k,q+=n){
            for(r=a+j*n,s=a+k*n,t=0.; r<p ;) t+= *r++ * *s++;
            *q-=t; *q/= *p;
        }
    }
    for(j=0,p=a; j<n ;++j,p+=n+1){
        for(k=0,q=a+j*n; k<j ;) b[j]-=b[k++]* *q++;
        b[j]/= *p;
    }
    for(j=n-1,p=a+n*n-1; j>=0 ;--j,p-=n+1){
        for(k=j+1,q=p+n; k<n ;q+=n) b[j]-=b[k++]* *q;
        b[j]/= *p;
    }
    return 0;
}

// solve Ax=b
int antislash_symmetric_positive_definite(float *x, float *A, float *b, int n)
{
	float *inA = malloc(n*n*sizeof*inA);
	float *inb = malloc(n*sizeof*inb);
	for (int i = 0; i < n*n; i++)
		inA[i] = A[i];
	for (int i = 0; i < n; i++)
		inb[i] = b[i];
	int r = solvps(inA, inb, n);
	for (int i = 0; i < n; i++)
		x[i] = inb[i];
	free(inA);
	free(inb);
	return r;
}

//create parametric transformation matrix (quadratic)
void fill_matrix_deg2(float *X, int i, int j)
{
    float tmp[] = { 1,i,j,0,0,0,i*i,i*j,0,0,0,1,i,j,i*j,j*j };
    for (int k = 0; k<16; k++)
    {
        X[k]=tmp[k];
    }
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

// compute the four directional derivatives
//
// in : input image
// w : width
// h : height
// out : four output images

void compute_directional_derivatives(float *im, int w, int h, float *out[4])
{
    float *dec[4];
    for (int k=0; k<4; k++)
    {
        out[k] = malloc(w*h*sizeof(float));
        dec[k] = malloc(w*h*sizeof(float));
    }


// create translated images
    apply_translation(dec[0], 1, 0, im, w, h);
    apply_translation(dec[1], 0, 1, im, w, h);
    apply_translation(dec[2], -1, 1, im, w, h);
    apply_translation(dec[3], -1, -1, im, w, h);

// compute derivatives
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            out[0][j*w+i] = (dec[0][j*w+i] - im[j*w+i])*(dec[1][j*w+i] - im[j*w+i]);
            out[1][j*w+i] = (dec[1][j*w+i] - im[j*w+i])*(dec[1][j*w+i] - im[j*w+i]);
            out[2][j*w+i] = (dec[2][j*w+i] - im[j*w+i])*(dec[2][j*w+i] - im[j*w+i])/2;
            out[3][j*w+i] = (dec[3][j*w+i] - im[j*w+i])*(dec[3][j*w+i] - im[j*w+i])/2;
        }
    
    free(dec[0]);
    free(dec[1]);
    free(dec[2]);
    free(dec[3]);
}

// subprogram to interpolate the pixels

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

// warp the image according to a 8-parameter transformation
//
// im : input image
// p : parameter
// w : width
// h : height
// out : warped image

void warping_8_parameters(float *im, float p[8], int w, int h, float *out)
{
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            float u;
            float v;
            u=p[0]+p[1]*i+p[2]*j+p[6]*i*i+p[7]*i*j;
            v=p[3]+p[4]*i+p[5]*j+p[6]*i*j+p[7]*j*j;
                                   
            out[j*w+i]=bicubic_interpolation(im, w, h, i+u, j+v);
                                   
        }
}

// compute the mean in a 11-pixel window centered on each pixel
//
// im : input image
// w : width
// h : height
// out : mean image
                                   
void mean_window_11(float *im, int w, int h, float *out)
{
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            float S=0;
        
            for(int k = -5; k <= 5 ; k++)
                for(int l = -5; l <= 5 ; l++)
                {
                    int ii=i+l;
                    int jj=j+k;
                    S = S + getpixel_0(im, w, h, ii, jj);
                }
            out[j*w+i]=S/121;
        }
}

//compute the norm 2 of the difference between a window and its mean for each pixel of an image
//
// im : input image
// mean : image of the means (obtained for example with mean_window_11)
// w : width
// h : height
//

void norm_2_window_11(float *im, float *mean, int w, int h, float *out)
{
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            float S=0;
            for(int k = -5; k <= 5 ; k++)
                for(int l = -5; l <= 5 ; l++)
                {
                    int ii=i+l;
                    int jj=j+k;
                    S = S + (getpixel_0(im, w, h, ii, jj)-getpixel_0(mean, w, h, i, j))*(getpixel_0(im, w, h, ii, jj)-getpixel_0(mean, w, h, i, j));
                }
            out[j*w+i]=sqrt(S);
        }
}

// compute the correlation surface between two images on every pixel for a given displacement
//
// im1, im2 : input images
// w : width
// h : height
// u,v : coordinates of the displacement
// out : image of the correlation surface

void correlation_surface(float *im1, float *im2, int w, int h, int u, int v, float *out)
{
    float *mean1 = malloc(w*h*sizeof(float));
    float *mean2 = malloc(w*h*sizeof(float));
    float *norm1 = malloc(w*h*sizeof(float));
    float *norm2 = malloc(w*h*sizeof(float));
    
    mean_window_11(im1, w, h, mean1);
    mean_window_11(im2, w, h, mean2);
    norm_2_window_11(im1, mean1, w, h, norm1);
    norm_2_window_11(im2, mean2, w, h, norm2);
    
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            float S=0;
            int ii=i+u;
            int jj=j+v;
            
            for(int k = -5; k <= 5 ; k++)
                for(int l = -5; l <= 5 ; l++)
                {
                    int i2=i+k;
                    int j2=j+l;
                    int ii2=ii+k;
                    int jj2=jj+l;
                    
                    S = S + (getpixel_0(im1, w, h, i2, j2)-getpixel_0(mean1, w, h, i, j))*(getpixel_0(im2, w, h, ii2, jj2)-getpixel_0(mean2, w, h, ii, jj));
                }
            
            out[j*w+i] = S/(norm1[j*w+i]*norm2[jj*w+ii]);
        }
    
    free(mean1);
    free(mean2);
    free(norm1);
    free(norm2);
}

//compute the gradient of the correlation surface in (0,0)
//
// im : input image
// w, h : dimensions
// out : gradient image (a 2x1 vector for each pixel)

void grad_correlation_surface(float *im1, float *im2, int w, int h, float *out[2])
{
    printf("Calcul du gradient... \n");
    
    float *surf00 = malloc(w*h*sizeof(float));
    float *surf01 = malloc(w*h*sizeof(float));
    float *surf10 = malloc(w*h*sizeof(float));
    
    correlation_surface(im1, im2, w, h, 0, 0, surf00);
    correlation_surface(im1, im2, w, h, 1, 0, surf10);
    correlation_surface(im1, im2, w, h, 0, 1, surf01);
    
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            out[0][j*w+i] = surf10[j*w+i]-surf00[j*w+i];
            out[1][j*w+i] = surf01[j*w+i]-surf00[j*w+i];
        }
    
    free(surf00);
    free(surf01);
    free(surf10);
}

// compute the hessian matrix of the correlation surface in (0,0)
//
// im : input image
// w, h : dimensions
// out : hessian image (a 2x2 matrix for each pixel)

void hessian_correlation_surface(float *im1, float *im2, int w, int h, float *out[4])
{
    printf("Calcul de la hessienne... \n");
    
    float *surf00 = malloc(w*h*sizeof(float));
    float *surf01 = malloc(w*h*sizeof(float));
    float *surf10 = malloc(w*h*sizeof(float));
    float *surf0_1 = malloc(w*h*sizeof(float));
    float *surf_10 = malloc(w*h*sizeof(float));
    float *surf11 = malloc(w*h*sizeof(float));
    
    correlation_surface(im1, im2, w, h, 0, 0, surf00);
    correlation_surface(im1, im2, w, h, 1, 0, surf10);
    correlation_surface(im1, im2, w, h, 0, 1, surf01);
    correlation_surface(im1, im2, w, h, -1, 0, surf_10);
    correlation_surface(im1, im2, w, h, 0, -1, surf0_1);
    correlation_surface(im1, im2, w, h, 1, 1, surf11);
    
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            out[0][j*w+i] = surf10[j*w+i]+surf_10[j*w+i]-2*surf00[j*w+i];
            out[1][j*w+i] = surf00[j*w+i]+surf11[j*w+i]-surf10[j*w+i]-surf01[j*w+i];
            out[2][j*w+i] = surf00[j*w+i]+surf11[j*w+i]-surf10[j*w+i]-surf01[j*w+i];
            out[3][j*w+i] = surf01[j*w+i]+surf0_1[j*w+i]-2*surf00[j*w+i];
        }
    
    free(surf00);
    free(surf01);
    free(surf10);
    free(surf0_1);
    free(surf_10);
    free(surf11);
    
}

// compute A (square matrix with length = number of parameters, here 8) with the hessian and the matrix X(i,j)
//

void A_matrix_8(float *im1, float *im2, int w, int h, float out[64])
{
    float *hess[4];
    
    for (int k=0; k<4; k++)
    {
        hess[k] = malloc(w*h*sizeof(float));
    }
    
    float X[16];
    
    hessian_correlation_surface(im1, im2, w, h, hess);
    
    for(int k = 0 ; k < 8 ; k++)
        for(int l = 0 ; l < 8 ; l++)
        {
            out[8*k+l] = 0;
        }
                
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            fill_matrix_deg2(X, i, j);
            
            for(int k = 0 ; k < 8 ; k++)
                for(int l = 0 ; l < 8 ; l++)
                {
                    out[8*k+l] = out[8*k+l] + X[k]*(X[l]*hess[0][j*w+i]+X[8+l]*hess[1][j*w+i])+X[8+k]*(X[l]*hess[2][j*w+i]+X[8+l]*hess[3][j*w+i]);
                }
        }
    
    free(hess[0]);
    free(hess[1]);
    free(hess[2]);
    free(hess[3]);
}


//compute B (vector with same length as the number of parameters, here 8) cith the gradient and the matrix X(i,j)
//

void B_matrix_8(float *im1, float *im2, int w, int h, float out[8])
{
    float *grad[2];
    grad[0] = malloc(w*h*sizeof(float));
    grad[1] = malloc(w*h*sizeof(float));
    float X[16];
    
    grad_correlation_surface(im1, im2, w, h, grad);
    
    for(int l = 0 ; l < 8 ; l++)
    {
        out[l] = 0;
    }
    
    for(int j = 0 ; j < h ; j++)
        for(int i = 0 ; i < w ; i++)
        {
            fill_matrix_deg2(X, i, j);
            
            for(int k = 0 ; k < 8 ; k++)
                {
                    out[k] = out[k] - grad[0][j*w+i]*X[k] - grad[1][j*w+i]*X[8+k];
                }
        }
    
    free(grad[0]);
    free(grad[1]);

}

//find displacement
//
// A, B : matrices computed with A_matrix and B_matrix

void delta_displacement_8(float *A, float *B, float *out)
{
    antislash_symmetric_positive_definite(out, A, B, 8);
}

// main program

int main(int argc, char **argv)
{
    // process input arguments
    if (argc != 4) {
        fprintf(stderr, "usage:\n \t imageIn1 imageIn2 imageWarped2 \n");
        //                         0     1       2        3
    }
    
    char *filename_ImgIn1 = argv[1];
    char *filename_ImgIn2 = argv[2];
    char *filename_ImgOut = argv[3];
    
    //read input images
    int w, h;
    float *im1 = iio_read_image_float(filename_ImgIn1, &w, &h);
    float *im2 = iio_read_image_float(filename_ImgIn2, &w, &h);
    
    //allocate space for derivatives and displacement
    float *der1[4];
    float *der2[4];
    
    for (int k=0; k<4; k++)
    {
        der1[k] = malloc(w*h*sizeof(float));
        der2[k] = malloc(w*h*sizeof(float));
    }
    
    float delta[8];
    
    printf("Calcul des derivees directionnelles... \n");
    
    //compute derivatives
    compute_directional_derivatives(im1, w, h, der1);
    compute_directional_derivatives(im2, w, h, der2);
    
    
    
    //compute A and B for each derivative and sum it
    float A[4][64];
    float B[4][8];
    float AA[64] ;
    float BB[8] ;
    
    printf("Calcul des matrices A et B... \n");
    
    for(int i=0 ; i < 4 ; i++)
    {
        A_matrix_8(der1[i], der2[i], w, h, A[i]);
        B_matrix_8(der1[i], der2[i], w, h, B[i]);
    }
    
    
    
    for(int j = 0 ; j < 8 ; j++)
        for(int i = 0 ; i < 8 ; i++)
        {
            AA[j*8+i]=A[0][j*8+i]+A[1][j*8+i]+A[2][j*8+i]+A[3][j*8+i];
           
        }
    for(int j=0 ; j<8 ; j++)
    {
        BB[j]=B[0][j]+B[1][j]+B[2][j]+B[3][j];
    }

    printf("Calcul du vecteur delta... \n");
    
    // compute displacement
    delta_displacement_8(AA, BB, delta);
    
    
    
    //allocate space for the output image
    float *out2 = malloc(w*h*sizeof(float));
    
    printf("Warping... \n");
    
    // warp image
    warping_8_parameters(im2, delta, w, h, out2);
    
    printf("Warping : done \n");
    
    //save the output image
    iio_save_image_float(filename_ImgOut, out2, w, h);
    
    //cleanup and exit
    free(out2);
    for (int k=0; k<4; k++)
    {
        free(der1[k]);
        free(der2[k]);
    }
    
    return 0;
}



