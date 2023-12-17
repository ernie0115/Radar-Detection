clc; 
clear; 
close all; 
dirVehicle = 'Camry'; 
numRangeBins = 128; 
numHeightBins = 128; 
M_x = numRangeBins; 
M_y = numRangeBins; 
M_z = numHeightBins; 

N  = [M_x M_y,M_z]; % [128, 128, 128] 

D_x = 14.3; % D_x is 14.3 meters. 
D_y = 14.3; % D_y is 14.3 meters. 
D_z = 14.3; % D_z is 14.3 meters. 

Res_x = D_x/M_x; % Res_x = 14.3/128; Res_x is the resolution per unit along x axis. 
Res_y = D_y/M_y; % Res_y = 14.3/128; Res_y is the resolution per unit along y axis. 
Res_z = D_z/M_z; % Res_z = 14.3/128; Res_z is the resolution per unit along z axis. 
lsize = D_x; % lsize = 14.3 

xImage = (-M_x/2 + 1:M_x/2)*Res_x;  % ((-64 + 1):(64))*(14.3/128) 
yImage = (-M_y/2 + 1:M_y/2)*Res_y;  % ((-64 + 1):(64))*(14.3/128) 
zImage = -(-M_z/2 + 1:M_z/2)*Res_z; % ((-64 + 1):(64))*(14.3/128) 

[xxx,yyy] = ndgrid(xImage,yImage); 
[xGrid,yGrid,zGrid] = ndgrid(xImage,yImage,zImage); 
xx_1 = zeros(size(xGrid)); 

counts = 1; 

% Extract the 72*72*72 cube which is located at the center part of 128*128*128 cube. 
% 72*72*72 cube is where the vehicle located. 
idxCount = 1; 
XX_New = zeros(512,357*729); % XX_New has the size of 512*260253; 357*729 = 260253 

for azCenter = 2:358 
    
    load(sprintf('%s%sFFT%sResults_%d.mat',dirVehicle,filesep,filesep,azCenter)); 
    xx = fftshift(X2_1); % Data of the vehicle at the spatial domain 
    
    % Segmentation of 128*128*128 voxel 
    center_x = 64; 
    center_y = 64; 
    center_z = 64; 
    cropSize = 72; 
    
    xx_Cropped = center_x - cropSize/2:center_x + cropSize/2 - 1; % (64 - 36):(64 + 36 - 1) 
    yy_Cropped = center_y - cropSize/2:center_y + cropSize/2 - 1; % (64 - 36):(64 + 36 - 1) 
    zz_Cropped = center_z - cropSize/2:center_z + cropSize/2 - 1; % (64 - 36):(64 + 36 - 1) 
    XX = xx(xx_Cropped, yy_Cropped, zz_Cropped); 
    
    % Segmentation of 72*72*72 voxel 
    
    for i = 1:9 
        for j = 1:9 
            for k = 1:9 
                
                XX_New(:,idxCount) = reshape(XX((i-1)*8 + 1:(i)*8,(j-1)*8 + 1:(j)*8,(k-1)*8 + 1:(k)*8),512,1); 
                idxCount = idxCount + 1; 
                
            end 
        end 
    end 
    
end 

% XX_New has the size of 512*(729*357), 512*260253. 
% Size of the 8*8*8 cube is transformed to 512*1. 
% (72/8)*(72/8)*(72/8), 729, is the total number of 8*8*8 cube. 
% 357 is the number of different azimuth angle. 
save(sprintf('%s%sFFT%sResults_Submatrix.mat',dirVehicle,filesep,filesep),'XX_New','-v7.3'); 









% [aa,bb,cc] = ind2sub(size(xx),find(20*log10(abs(xx)) > max(20*log10(abs(xx(:)))) - 30)); 
% coordinates1 = []; 
% 
% coordinates1(1,:) = arrayfun(@(i)xGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
% coordinates1(2,:) = arrayfun(@(i)yGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
% coordinates1(3,:) = arrayfun(@(i)zGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
% amps = 20*log10(abs(arrayfun(@(i)xx(aa(i),bb(i),cc(i)),(1:length(aa)).'))); 
% 
% 
% figure;scatter3(coordinates1(1,:),coordinates1(2,:),coordinates1(3,:),30,'g','filled'); 
% xlim([-4,4]) 
% ylim([-4,4]) 
% zlim([-4,4]) 
% axis square 
%
% title(sprintf('Vehicle = %s, Viewing angle = %d',dirVehicle,azCenter)); 
% xlabel('X'); 
% ylabel('Y'); 
% zlabel('Z'); 
% 