clc; 
clear; 

addpath('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\spgl1-2.0'); % directory path of spgl1 
addpath('C:\Users\USER\Downloads\Ohio State University\OSU 課程\Spring 2022\ECE 6193\NUFFT_code'); % directory path of NUFFT 
dirVehicle = 'Camry'; % Name of the directory of the vehicle 
elev = [41:0.0625:44.9375];  % The total coverage of measurements in elevation domain (degrees) [64 different elevations] 
cspeed = 299792458; % Velocity of light (c is approximately 2.99752458*10^8 m/s) 
bandwidth = 670e6;  % Bandwidth of the signal (BW = 0.67G Hz) 
azSpan = 4.1;  % The total coverage of measurements in azimuth domain (degrees) 

% Create a 128 x 128 x 128 3D image of the object. 

numRangeBins = 128; numHeightBins = 128; M_x = numRangeBins; M_y = numRangeBins; M_z = numHeightBins; 
N  = [M_x M_y,M_z]; % [128, 128, 128] 

D_x = 14.3; % Size of the scene in X-axis (14.3 m) (128 pixels) 
D_y = 14.3; % Size of the scene in Y-axis (14.3 m) (128 pixels) 
D_z = 14.3; % Size of the scene in Z-axis (14.3 m) (128 pixels) 

Res_x = D_x/M_x; % Resolution of the object in X-axis (14.3/128), 0.112 
Res_y = D_y/M_y; % Resolution of the object in Y-axis (14.3/128), 0.112 
Res_z = D_z/M_z; % Resolution of the object in Z-axis (14.3/128), 0.112 

lsize = D_x; 

xImage = (-M_x/2 + 1:M_x/2)*Res_x; % ((-64 + 1):(64))*(14.3/128) 
yImage = (-M_y/2 + 1:M_y/2)*Res_y; % ((-64 + 1):(64))*(14.3/128) 
zImage = (-M_z/2 + 1:M_z/2)*Res_z; % ((-64 + 1):(64))*(14.3/128) 

numElevationPasses = 64; % There are 64 different elevations. 
fMin = 9.2e9; % Minimum Frequency is 9.2G Hz. 
fMax = fMin + bandwidth; % Maximum Frequency is 9.2G + 0.67G = 9.87G Hz. 

[xGrid,yGrid,zGrid] = ndgrid(xImage,yImage,zImage); % Not Necessary 
[zGrid1,~,~] = ndgrid(zImage,yImage,xImage); % Not Necessary 
[xGrid1,yGrid1] = meshgrid(xImage,yImage); % Not Necessary 

kx_grid = linspace(-1/2/Res_x, 1/2/Res_x, M_x + 1); 
kx_grid(end) = []; 
% Create 129 evenly spaced points between (-1/2/(14.3/128)), -4.475, and (1/2/(14.3/128)), 4.475. 
% And then, delete the 129th element. So kx_grid contains 128 evenly spaced points between -4.475 to 4.4056. 

ky_grid = linspace(-1/2/Res_y, 1/2/Res_y, M_y + 1); 
ky_grid(end) = []; 
% Create 129 evenly spaced points between (-1/2/(14.3/128)), -4.475, and (1/2/(14.3/128)), 4.475. 
% And then, delete the 129th element. So ky_grid contains 128 evenly spaced points between -4.475 to 4.4056. 

kz_grid = linspace(-1/2/Res_z, 1/2/Res_z, M_z + 1); 
kz_grid(end) = []; 
% Create 129 evenly spaced points between (-1/2/(14.3/128)), -4.475, and (1/2/(14.3/128)), 4.475. 
% And then, delete the 129th element. So kz_grid contains 128 evenly spaced points between -4.475 to 4.4056. 

% Load the measurements from the directory. 
for i = 1:numElevationPasses 
    
    passes(i) = load(sprintf('%s%s%s_el%.4f.mat',dirVehicle,filesep,dirVehicle,elev(i))); 
    
end 

% Looping over the azimuth angles from 2 degrees to 358 degrees 
for azCenter = 2:358 
    
    k_x_total = []; 
    k_y_total = []; 
    k_z_total = []; 
    phTotal1 = []; 
    samplesIndexPass = 1; 
    
    % Loading the measurements spanning the required azimuth and elevation angles 
    for i = 1:numElevationPasses 
        
        maxAz = azCenter + azSpan/2; % Maximum Azimuth Angle inside a certain azCenter 
        minAz = azCenter - azSpan/2; % Minimum Azimuth Angle inside a certain azCenter       
        idxAz =  find((double(passes(i).data.azim) < maxAz) & (double(passes(i).data.azim) >= minAz)); 
        % Find the angles between minAz and maxAz. 
        % The size of idxAz is 1*65. 
        [azimuthVals1{i},idxSorted] = sort(double(passes(i).data.azim(idxAz))); 
        % Sort all columns in passes(i).data.azim(idxAz). 
        % azimuthVals1{i} is the sorted version of passes(i).data.azim(idxAz). 
        % idxSorted is the index of elements in azimuthVals1{i} according to passes(i).data.azim(idxAz). 
        azimuthVals1{i} = azimuthVals1{i} - mean(azimuthVals1{i}); 
        idxAz = idxAz(idxSorted); 
        
        elev1{i} = elev(i); % Elevation angle of each pass 
        
        ff = double(passes(i).data.FGHz)*1e9; 
        f1{i} = ff((ff > fMin) & (ff < fMax)); % Frequency used in each pass (9.2GB < f1 < 9.87GB) 
        idxF = find((ff > fMin & ff < fMax)); % The size of idxF is 64*1. (219 - 282) 
        
        numFreqs(i) = length(f1{i}); 
        numPulses(i) = length(idxAz); 
        
        pp = (double(passes(i).data.hh(idxF,idxAz))); 
        phTotal1(:,:,numElevationPasses - i + 1) = pp; % Size: 64*65*64 
        samplesIndexPass(i + 1) = length(f1{i})*length(azimuthVals1{i}); 
        meanF(i) = mean(f1{i}); 
        meanElev(i) = mean(elev1{i}); 
        
    end 
    
    fixedCenterFreq = mean(meanF); 
    fixedElev = mean(meanElev); 
    azCenterRadarCentric = 0; 
    phaseRampCenter_y = 2*cosd(fixedElev)*sind(azCenterRadarCentric)*fixedCenterFreq/cspeed; 
    phaseRampCenter_x = 2*cosd(fixedElev)*cosd(azCenterRadarCentric)*fixedCenterFreq/cspeed; 
    phaseRampCenter_z = 2*sind(fixedElev)*fixedCenterFreq/cspeed; 
    
    % Generating the spatial frequency points 
    for i = 1:numElevationPasses 
        
        frep = repmat(f1{i},1,length(azimuthVals1{i})); 
        % azimuthVals1{i}: 1*65; length(azimuthVals1{i}) = 65 
        
        azrep = repmat(((azimuthVals1{i})),length(f1{i}),1); 
        elrep = repmat(((elev1{i})),length(f1{i}),length(azimuthVals1{i})); 
        
        k_y = 2/cspeed*sind(azrep).*(frep).*cosd(elrep); 
        k_y_c = phaseRampCenter_y; 
        k_y = k_y - k_y_c; 
        
        k_x = 2/cspeed*cosd(azrep).*(frep).*cosd(elrep); 
        k_x_c = phaseRampCenter_x; 
        k_x = k_x - k_x_c; 
  
        k_z = 2/cspeed*(frep).*sind(elrep); 
        k_z_c = phaseRampCenter_z; 
        k_z = k_z - k_z_c; 
        
        k_x_total(:,:,numElevationPasses - i + 1) = k_x; 
        k_y_total(:,:,numElevationPasses - i + 1) = k_y; 
        k_z_total(:,:,numElevationPasses - i + 1) = k_z; 
        
    end 
    
    [kx1,ky1,kz1] = ndgrid(kx_grid,ky_grid,kz_grid); 
    F = scatteredInterpolant(k_x_total(:),k_y_total(:),k_z_total(:),phTotal1(:),'linear','none'); 
    phInterp = F(kx1,ky1,kz1); 
    phInterp(isnan(phInterp)) = 0; 
    
    X2_1 =  ifftn((phInterp)); 
    
    % Saving 3D representation 
    save(sprintf('%s%sFFT%sResults_%d',dirVehicle,filesep,filesep,... 
        azCenter),'X2_1','azCenter','phInterp'); 
    
end 
