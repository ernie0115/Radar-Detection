% In real life, we just only have P and D. (P: phInterp, D: Dksvd1, Dksvd2,Dksvd3, Dksvd4.) 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\Camry\FFT\Results_137.mat'); % 1 image of Camry, from "Results_137.mat"; 
% (F_3D^(-1))P = y; y - "X2_1"; P - "phInterp"; azCenter = 137; phInterp and X2_1 - 128*128*128; 
% Refer to the data from processData_CVDomes_Yang_LeastSquare.mat 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\T16\D4T16.mat'); % Dksvd4 has the size of 512*8192 when T = 16. (Trained Dictionary by K-SVD) 
addpath('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\compbox10'); 
addpath('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\spgl1-2.0'); 

x4 = zeros(size(Dksvd4,2),4096); % Ultimately, x4 should be 8192*4096. 
r4 = randi([1 8192],16,4096); % Select 16*4096 positions in x4 for 16*4096 nonzero entries. 
x4(r4) = 1; % 16*4096 nonzero entries in x4 are 1. (T = 16) 
x4 = x4(:); 
y = forward_operator(Dksvd4, x4,16,16,16,128,128,128,1); 

SNR = 30; % Signal to Noise Ratio 
New_x4 = norm(y,'fro'); 
sigma = New_x4*10^(-SNR/20); 

New_phInterp = reshape(phInterp, 2097152,1); 

A = @(x,mode)forward_operator(Dksvd4,x,16,16,16,128,128,128,mode); 
options = spgSetParms('isComplex','true'); 
[X4] = spg_bpdn(A,New_phInterp,sigma,options); 

save test_data_D4_Camry_SNR_20_azCenter_137.mat 