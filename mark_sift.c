//
//  mark_sift.c
//  
//
//  Created by Clément on 03/05/2014.
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

// function to put a black pixel in a white square somewhere in the image

void mark_pixel(float *out, int w, int h, int x, int y)
{
    out[y*w+x] = 255;
    out[(y-1)*w+x-1] = 0;
    out[(y-1)*w+x] = 0;
    out[(y-1)*w+x+1] = 0;
    out[y*w+x-1] = 0;
    out[y*w+x+1] = 0;
    out[(y+1)*w+x-1] = 0;
    out[(y+1)*w+x] = 0;
    out[(y+1)*w+x+1] = 0;
}

int count_line(FILE* fichier)
{
    int line = 0;
    int c=0;
    while ((c=fgetc(fichier)) != EOF )
    {
        if(c== '\n')
        {
            line++;
        } 
    } 
    return line; 
}


// function to read the file with the coordinates of the SIFT points

void read_coordinates(char *filename[], float coordx1[256], float coordy1[256], float coordx2[256], float coordy2[256], int size)
{
    FILE* file = NULL;
    
    file = fopen(*filename, "r");
    
    if (file != NULL)
    {
        int i = 0;
        size = 0;
        int c = fgetc(file);
        rewind(file);

        while(c != EOF)
        {
            fscanf(file, "%f %f %f %f", &coordx1[i], &coordy1[i], &coordx2[i], &coordy2[i]);
            
            c = fgetc(file);
            i=i+1;
            size = size+1;
        }
        
        fclose(file);
    
    }
    
    else
    {
        printf("Impossible to open file");
    }
}

// main fucntion

int main(int argc, char **argv)
{
    // process input arguments
    if (argc != 6)
    {
        fprintf(stderr, "usage: \n \t image1 image2 inliers out1 out2");
        //                         0  1      2      3       4    5
        return 0;
    }
    
    char *filename_img1 = argv[1];
    char *filename_img2 = argv[2];
    char *filename_inliers = argv[3];
    char *filename_out1 = argv[4];
    char *filename_out2 = argv[5];

    int w,h;
    float *out1 = iio_read_image_float(filename_img1, &w, &h);
    float *out2 = iio_read_image_float(filename_img2, &w, &h);
    
    // read coordinates
    
    float coordx1[256];
    float coordy1[256];
    float coordx2[256];
    float coordy2[256];
    int size;
    
    FILE* file = NULL;
    
    file = fopen(filename_inliers, "r");
    
    int n = count_line(file);
    
    fclose(file);
    
    read_coordinates(&filename_inliers, coordx1, coordy1, coordx2, coordy2, size);
    printf("size = %d \n", n);
    
    
    
    for (int i = 0 ; i < n ; i++)
    {
        mark_pixel(out1, w, h, floor(coordx1[i]), floor(coordy1[i]));
        mark_pixel(out2, w, h, floor(coordx2[i]), floor(coordy2[i]));
    }

    //save outputs
    
    iio_save_image_float(filename_out1, out1, w, h);
    iio_save_image_float(filename_out2, out2, w, h);
    
    return 0;
    
}













