# (A) processData_CVDomes_Yang_LeastSquare.m 
## **Objective:** 
It generates the automobile frequency domain data, named "phInterp". 

## **Annotation:** 
- "phInterp" size: 128 x 128 x 128
- "phInterp" actual distance: 14.3m x 14.3m x 14.3m
- kx = 4 * (pi) * f / c * cos(azimuth_Center_angle) * cos(elevation_angle) 
- ky = 4 * (pi) * f / c * sin(azimuth_Center_angle) * cos(elevation_angle) 
- kz = 4 * (pi) * f / c * sin(elevation_angle)
- f means SAR frequency. c means light velocity (Unit: m/s). kz, ky, and kz means spatial frequency components along x, y, and z directions (Unit: cycle / m). 

# (B) plotStuff_yang.m 
## **Objective:** 
Extract the center 72 x 72 x 72 cube from automobile frequency domain data, sized 128 x 128 x 128 cube. 

## **Annotation:**  
- Since "phInterp" is centered at the center 72 x 72 x 72 cube, we extract the center cube to do K-SVD Dictionary Training. 

# (C) myfile_train_dictionary.m 
## **Objective:** 
Train the Dictionary based on "phInterp" by K-SVD Algortihm. 

## **Annotation:** 
- Currently, K-SVD is the most popular method to do Dictionary Training.
- It takes around 1 day to train 10 automobiles. 

# (D) Test_phInterp_D4_azCenter_137.m 
## **Objective:** 
Use "Basis Pursuit Denoise" to find Sparse Coefficient based on "phInterp" and the trained Dictionary by K-SVD Algorithm. 

## **Annotation:** 
- Advantage: It has the best performance compared to "Compressive Sampling Matching Pursuit" visually and quantitively (MSE = 0.01).
- Disadvantage: It takes around 10 hours to get the sparse coefficient. 

# (E) Test_Image_D4_azCenter_137_s.m 
## **Objective:** 
Prove Sparse Coefficient Accuracy visually and quantitively (MSE). 

## **Annotation:** 
- Currently, "Basis Pursuit Denoise" has the best performance. 
