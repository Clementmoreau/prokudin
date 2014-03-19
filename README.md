prokudin
========
These programs are designed to register pictures of the Prokudin-Gorskii collection.
http://www.loc.gov/pictures/collection/prok

* gauss : apply a gaussian filter (sigma = 0.8) to an image.
   
usage : gauss imgIn.png imgOut.png
   
* cut3 : cut a gray image into three equal height vertical slices
   
usage : cut3 input.tiff out1.tiff out2.tiff out3.tiff


* join3 : join three gray images of the same size to form a color image
   
usage : join3 in_red.tiff in_green.tiff in_blue.tiff out_colorized.tiff


* quantize : quantize an image to 8 bits
   
usage : quantize 0 65000 in.tiff out.png


* translation : apply an integral displacement to an image (filling-in by zero)
   
usage : translation 10 -20 in.png out.png


* registration : naive algorithm to register two gray images
   
usage : registration A.tiff B.tiff Btranslated.tiff
   
* pyramide : return seven images, each one zommed out by a factor two with a gaussian filter
   
usage : pyramide imgIn.png pyr1.png pyr2.png pyr3.png pyr4.png pyr5.png pyr6.png pyr7.png
   
* demo_orsa_homography : find the best homography to fit two images with the SIFT matching.
  
This program has been adapted from http://www.ipol.im/pub/art/2012/mmm-oh/
   
usage : demo_orsa_homography imgInA.png imgInB.png allMatches.txt inMatches.txt imgAWarped.png imgBWarped.png
   
* main : compute the four directional derivatives of an image (squared)
   
usage : main imgIn.png der1.png der2.png der3.png der4.png

* irani : registration algorithm based on the Irani article (Robust Multi-Sensor Image Alignment)

usage : irani imgIn1.png imgIn2.png imgWarped2.png

* homographie.sh : script which takes the original glass negative and return the registered image (with orsa-homography method)

* irani1.sh : script which takes the original glass negative, does the homography registration, and the irani warping. This may be very long for big images.



