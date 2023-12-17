% 
% Inputs 
% D - Dictionary matrix 
% x - Spatial domain Dictionary coefficients if mode == 1 and Fourier domain 
%    measurement if mode == 2 
% nCubes_x - number of small cubes in x dimension 
% nCubes_y - number of small cubes in y dimension 
% nCubes_z - number of small cubes in z dimension 
% nX - size of x dimension of Image 
% nX - size of y dimension of Image 
% nZ - size of z dimension of Image 
% mode  - 1 : Forward operator and x is in spatial domain 
%       - 2 : Adjoint of Forward operator and x is in Frequecy domain 

function y = forward_operator(D,x,nCubes_x,nCubes_y,nCubes_z,nX,nY,nZ,mode) 
 
if mode == 1 
    
    % Compute the Forward operation which takes the 3D FFT 
    
    x = reshape(x,size(D,2),nCubes_x*nCubes_y*nCubes_z); % 8192*4096 

    Dx = D*x; 
    smallCubeSize = 8; 
    
    yTemp = zeros(nX,nY,nZ); % yTemp has the size of 128*128*128. 
    
    for i = 1:nCubes_x 
        for j = 1:nCubes_y 
            for k = 1:nCubes_z 
                
                yTemp((i-1)*smallCubeSize + 1:i*smallCubeSize,...
                    (j-1)*smallCubeSize + 1:j*smallCubeSize,...
                    (k-1)*smallCubeSize + 1:k*smallCubeSize) = ...
                    reshape(Dx(:,256*(i - 1) + 16*(j - 1) + k),smallCubeSize,...
                    smallCubeSize,smallCubeSize); 
                
            end 
        end 
    end 
      
    y = yTemp; 
    
    y = 1/sqrt(nX*nY*nZ)*fftn(yTemp); % y means the Fourier Domain of the automobile data. 
    
    y = y(:); % y will be changed to 2097152*1 from 128*128*128. 
    
else 
    
    % Compute the Adjoint of the Forward operator which returns the spatial domain estimate 
    
    smallCubeSize = 8; 
    x = reshape(x,nX,nY,nZ); 
    y = zeros(size(D,2),nCubes_x*nCubes_y*nCubes_z); % 8192*4096 
    yTemp = sqrt(nX*nY*nZ)*ifftn(x); 
    
    for i = 1:nCubes_x 
        for j = 1:nCubes_y 
            for k = 1:nCubes_z 
                
                yy = yTemp((i-1)*smallCubeSize+1:i*smallCubeSize,...
                    (j-1)*smallCubeSize+1:j*smallCubeSize,...
                    (k-1)*smallCubeSize+1:k*smallCubeSize); 
                yy = yy(:); 
                y(:,256*(i - 1) + 16*(j - 1) + k) = D'*yy; 
                
            end 
        end 
    end 
    
    y = y(:); % y will be changed to 33554432*1 from 8192*4096. 
    
end 

end 