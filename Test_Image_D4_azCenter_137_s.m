load('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\test_data_D4_Camry_SNR_20_azCenter_137.mat'); % X4 - 33554432*1 (8192*4096) 

% Original data (xx) 
xx = fftshift(X2_1); 
% X2_1_abs = abs(xx); 

x = linspace(-7.15,7.15,size(xx, 1)); 
y = linspace(-7.15,7.15,size(xx, 1)); 
z = linspace(-7.15,7.15,size(xx, 1)); 

% Create the coordinates of each voxel and store in variable xGrid, yGrid and zGrid 

[xGrid, yGrid, zGrid] = ndgrid(x, y, z); 
SNR_threshold = 30; 

% Threshold the image to get a sparser and cleaner image 
% We get the index locations of each dominant point in spatial domain. aa, bb and cc 

[aa,bb,cc] = ind2sub(size(xx), find(20*log10(abs(xx)) > max(20*log10(abs(xx(:)))) - SNR_threshold)); 
coordinates1 = []; 

% Looping through the detected points to obtain the corresponding coordinates of those dominant points 

coordinates1(1,:) = arrayfun(@(i)xGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
coordinates1(2,:) = arrayfun(@(i)yGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
coordinates1(3,:) = arrayfun(@(i)zGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
amps = 20*log10(abs(arrayfun(@(i)xx(aa(i),bb(i),cc(i)),(1:length(aa)).'))); 

figure(1);scatter3(coordinates1(1,:), coordinates1(2,:), coordinates1(3,:), 30, amps, 'filled'); 
colorbar;colormap jet; 

% Original data (X2_1_abs) 
X2_1_abs = abs(fftshift(X2_1)); 

x = linspace(-7.15,7.15,size(X2_1_abs, 1)); 
y = linspace(-7.15,7.15,size(X2_1_abs, 1)); 
z = linspace(-7.15,7.15,size(X2_1_abs, 1)); 

% Create the coordinates of each voxel and store in variable xGrid, yGrid and zGrid 

[xGrid, yGrid, zGrid] = ndgrid(x, y, z); 
SNR_threshold = 30; 

% Threshold the image to get a sparser and cleaner image 
% We get the index locations of each dominant point in spatial domain. aa, bb and cc 

[aa,bb,cc] = ind2sub(size(X2_1_abs), find(20*log10(abs(X2_1_abs)) > max(20*log10(abs(X2_1_abs(:)))) - SNR_threshold)); 
coordinates2 = []; 

% Looping through the detected points to obtain the corresponding coordinates of those dominant points 

coordinates2(1,:) = arrayfun(@(i)xGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
coordinates2(2,:) = arrayfun(@(i)yGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
coordinates2(3,:) = arrayfun(@(i)zGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
amps2 = 20*log10(abs(arrayfun(@(i)X2_1_abs(aa(i),bb(i),cc(i)),(1:length(aa)).'))); 

figure(2);scatter3(coordinates2(1,:), coordinates2(2,:), coordinates2(3,:), 30, amps2, 'filled'); 
colorbar;colormap jet; 

% My trained data 
SNR_threshold = 30; 

X4 = reshape(X4, 8192, 4096); 
Image4 = Dksvd4 * X4; 

Y4 = zeros(128, 128, 128); 

for i = 1:16 
    for j = 1:16 
        for k = 1:16 
            
            Y4(8*(i - 1) + 1:8*i,8*(j - 1) + 1:8*j,8*(k - 1) + 1:8*k) = reshape(Image4(:,256*(i - 1) + 16*(j - 1) + k), 8, 8, 8); 
            
        end 
    end 
end 

New_Image4 = abs(fftshift(Y4)); 

[aa,bb,cc] = ind2sub(size(New_Image4),find(20*log10(abs(New_Image4)) > max(20*log10(abs(New_Image4(:)))) - SNR_threshold)); 
coordinates4_image4 = []; 

% Looping through the detected points to obtain the corresponding coordinates of those dominant points 

coordinates4_image4(1,:) = arrayfun(@(i)xGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
coordinates4_image4(2,:) = arrayfun(@(i)yGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
coordinates4_image4(3,:) = arrayfun(@(i)zGrid(aa(i),bb(i),cc(i)),(1:length(aa)).'); 
amps_image4 = 20*log10(abs(arrayfun(@(i)New_Image4(aa(i),bb(i),cc(i)),(1:length(aa)).'))); 

figure(3);scatter3(coordinates4_image4(1,:), coordinates4_image4(2,:), coordinates4_image4(3,:), 30, amps_image4, 'filled'); 
colorbar;colormap jet; 

% Assuming you have two cubes called cube1 and cube2 
cube1 = xx; 
cube2 = New_Image4; 
cube3 = X2_1_abs; 

% Step 1: Ensure both cubes have the same size 
assert(isequal(size(cube1), size(cube2)), 'Cubes must have the same size.'); 
assert(isequal(size(cube3), size(cube2)), 'Cubes must have the same size.'); 

% Step 2: Calculate the difference cube 
diff_cube1_2 = cube1 - cube2; 
diff_cube3_2 = cube3 - cube2; 

% Step 3: Square each element of the difference cube 
diff_cube_squared1_2 = diff_cube1_2.^2; 
diff_cube_squared3_2 = diff_cube3_2.^2; 

% Step 4: Sum all the squared differences 
sum_squared_diff1_2 = sum(diff_cube_squared1_2(:)); 
sum_squared_diff3_2 = sum(diff_cube_squared3_2(:)); 

% Step 5: Calculate the mean square error 
mse1_2 = sum_squared_diff1_2 / numel(cube1); 
mse3_2 = sum_squared_diff3_2 / numel(cube3); 

% Display the result 
fprintf('The Mean Square Error between xx and New_Image4 is: %f.\n', mse1_2) 
fprintf('The Mean Square Error between X2_1_abs and New_Image4 is: %f.\n', mse3_2) 