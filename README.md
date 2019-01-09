# SABERProbMapSegmentation
The main function is saberProbMapSegmentation.m

This takes class probability maps of nuclei contours derived from a U-Net architecture for pixel classification as input. The center of the nuclei contour is found and used as seeds for marker controlled watershed segmentation. Some default parameters include nuclei and cytoplasm channel as 2nd and 3rd channel respectively, and using the distance transform to approximate the cytoplasm.

Please download bioformats_package.jar and include in the Matlab path if you don't already have this. 

Developed by:
Clarence Yapp, 2018
