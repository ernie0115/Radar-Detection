function [x_hat] = omp(A, y, ks) 
    % OMP algorithm for sparse signal recovery 
    % A: Sensing matrix; y: Measurements; k: Sparsity level 

    smallCubeSize = 8; nCubes_x = 16; nCubes_y = 16; nCubes_z = 16; nX = 128; nY = 128; nZ = 128; 
    y_spatial = sqrt(nX * nY * nZ) * ifftn(y); % size(y_spatial) = 128 * 128 * 128 
    y_s_2D = zeros(smallCubeSize * smallCubeSize * smallCubeSize, nCubes_x * nCubes_y * nCubes_z); % size(y_s_2D) = 512 * 4096 

    for i = 1:nCubes_x 
        for j = 1:nCubes_y 
            for k = 1:nCubes_z 
                
                y_s_2D(:, nCubes_y * nCubes_z * (i - 1) + nCubes_z * (j - 1) + k) = ... 
                reshape(y_spatial((i-1)*smallCubeSize+1:i*smallCubeSize,...
                    (j-1)*smallCubeSize+1:j*smallCubeSize,...
                    (k-1)*smallCubeSize+1:k*smallCubeSize), 512, 1); 

            end 
        end 
    end 

    r = y_s_2D; % Residual initialization (Size: 512 * 4096) 
    Omega = []; % Support set initialization 
    x_hat = zeros(size(A, 2), nCubes_x * nCubes_y * nCubes_z); % Sparse coefficient initialization (Size: 8192 * 4096) 

    % whos y_s_2D; disp(size(y_s_2D)); whos A; whos r; disp(size(A)); disp(size(r)); 

    for iter = 1:ks 
        % if ~ismatrix(A') 
        %     error('ctranspose(A) is not a matrix.'); 
        % end 
        % 
        % if ~ismatrix(r) 
        %     error('r is not a matrix.'); 
        % end 
        % 
        % if ~ismatrix(y_s_2D) 
        %     error('y_s_2D is not a matrix.'); 
        % end 
        % 
        % if size(A', 2) ~= size(r, 1) 
        %     error('Input dimensions are not compatible for matrix multiplication'); 
        % end 

        % Compute projection of residual onto columns of A 
        proj = abs(A' * r); % size(proj) = 8192 * 4096 

        % Find index of maximum projection (Find max number in each column.) 
        [~, idx] = max(proj); 
        
        % Update support set 
        Omega = [Omega, idx]; % size(Omega) = 1 * 4096, 1 * 8192, 1 * 12288, 1 * 16384 
        
        % Construct estimate of the support 
        Psi = A(:, Omega); % size(Psi) = 512 * 4096, 512 * 8192, 512 * 12288, 512 * 16384 

        % Solve least squares problem
        x_tilde = pinv(Psi) * y_s_2D; % size(x_tilde) = 4096 * 4096, 8192 * 4096, 12288 * 4096, 16384 * 4096 
        
        % Update sparse coefficient estimate
        x_hat(Omega, :) = x_tilde; 
        
        Dx = A * x_hat; 
        yTemp = zeros(nX,nY,nZ); % yTemp has the size of 128 * 128 * 128. 
    
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

        % Update residual 
        r_3D = y - 1/sqrt(nX * nY * nZ) * fftn(yTemp); 

        for i = 1:nCubes_x 
            for j = 1:nCubes_y 
                for k = 1:nCubes_z 
                    
                    r(:, nCubes_y * nCubes_z * (i - 1) + nCubes_z * (j - 1) + k) = ... 
                    reshape(r_3D((i-1)*smallCubeSize+1:i*smallCubeSize,...
                        (j-1)*smallCubeSize+1:j*smallCubeSize,...
                        (k-1)*smallCubeSize+1:k*smallCubeSize), 512, 1); 
    
                end 
            end 
        end 

    end 

    % disp(size(find(Omega - ones(1, 4096 * 4)))); % 16384 out of 16384 elements aren't 1. 
    % disp(size(unique(Omega))); % 2299 different atoms out of 8192 atoms are chosen from Dksvd4. 

end 