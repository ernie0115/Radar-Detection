% Train the dictionary by K-SVD Algorithm! 
% compbox10 contains the cksvd(). It can do dictionary training by K-SVD Algorithm. 
addpath('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\compbox10'); 

% Y is the training target of K-SVD Algorithm. (Y: 512*(729*357*10) - 512*2602530) 
% Size of the 8*8*8 cube is transformed to 512*1. 
% (72/8)*(72/8)*(72/8), 729, is the total number of 8*8*8 cube. 
% 357 is the number of different azimuth angle. 
% 10 means 10 different vehicles. 
Y = []; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\Camry\FFT\Results_Submatrix'); 
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\HondaCivic4dr\FFT\Results_Submatrix'); 
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\Jeep93\FFT\Results_Submatrix'); 
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\Jeep99\FFT\Results_Submatrix'); 
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\Maxima\FFT\Results_Submatrix'); 
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\MazdaMPV\FFT\Results_Submatrix'); 
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\Mitsubishi\FFT\Results_Submatrix'); 
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\Sentra\FFT\Results_Submatrix');
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\ToyotaAvalon\FFT\Results_Submatrix');
Y = [Y XX_New]; 
load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\CV Domes\ToyotaTacoma\FFT\Results_Submatrix'); 
Y = [Y XX_New]; 

% data is training target of K-SVD Algorithm. (Y: 512*729*357*10) 
% Tdata is the number of nonzero element in each column, 512*1. 
% Dictionary Size: 512*1024, 512*2048, 512*4096, 512*8192 (512*dictsize) 
params1.data = Y; params2.data = Y; params3.data = Y; params4.data = Y; 
params1.Tdata = 4; params2.Tdata = 4; params3.Tdata = 4; params4.Tdata = 4; 
params1.dictsize = 1024; params2.dictsize = 2048; params3.dictsize = 4096; params4.dictsize = 8192; 

[Dksvd1,g1,err1] = cksvd(params1,''); 
[Dksvd2,g2,err2] = cksvd(params2,''); 
[Dksvd3,g3,err3] = cksvd(params3,''); 
[Dksvd4,g4,err4] = cksvd(params4,''); 

save D1T4.mat Dksvd1; % Dksvd1 is the dictionary trained by K-SVD Algorithm. (Dksvd1 - 512*1024, T = 4) 
save x1T4.mat g1; 
save D2T4.mat Dksvd2; % Dksvd2 is the dictionary trained by K-SVD Algorithm. (Dksvd2 - 512*2048, T = 4) 
save x2T4.mat g2; 
save D3T4.mat Dksvd3; % Dksvd3 is the dictionary trained by K-SVD Algorithm. (Dksvd3 - 512*4096, T = 4) 
save x3T4.mat g3; 
save D4T4.mat Dksvd4; % Dksvd4 is the dictionary trained by K-SVD Algorithm. (Dksvd4 - 512*8192, T = 4) 
save x4T4.mat g4; 