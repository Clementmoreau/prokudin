#!/bin/sh

cut3 $1 blue.tiff green.tiff red.tiff

quantize 0 65000 blue.tiff blue.png
quantize 0 65000 red.tiff red.png
quantize 0 65000 green.tiff green.png

demo_orsa_homography red.png green.png all.txt in.txt redh.png greenh.png
demo_orsa_homography blue.png green.png all.txt in.txt blueh.png greenh.png

join3 redh.png green.png blueh.png $2