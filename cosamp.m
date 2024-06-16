function [x_hat] = cosamp(A, y, ks, max_iter) 
    % CoSaMP algorithm for sparse signal recovery 
    % A: Sensing matrix 
    % y: Measurements 
    % k: Sparsity level 
    % max_iter: Maximum number of iterations 
    
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

    r = y_s_2D; % Residual initialization 
    x_hat = zeros(size(A, 2), nCubes_x * nCubes_y * nCubes_z); % Sparse coefficient initialization (Size: 8192 * 4096) 

    for iter = 1:max_iter 

        % Estimate support set 
        [~, idx] = sort(abs(A' * r), 'descend'); 
        Omega = idx(1:2*ks); 
        
        % Construct estimate of the support 
        Psi = A(:, Omega); 
        
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
                    
                    r_temp((i-1)*smallCubeSize + 1:i*smallCubeSize,...
                        (j-1)*smallCubeSize + 1:j*smallCubeSize,...
                        (k-1)*smallCubeSize + 1:k*smallCubeSize) = ...
                        reshape(r(:,256*(i - 1) + 16*(j - 1) + k),smallCubeSize,...
                        smallCubeSize,smallCubeSize); 
                    
                end 
            end 
        end 

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

        disp(iter); 
        % Calculate the Frobenius of "old r"! (3D) 
        frobenius_norm_old_r = norm(r_temp, 'fro'); 
        disp(frobenius_norm_old_r); 

        % Compare the Frobenius of "old r and new r"! (3D) 
        frobenius_norm_delta_r = norm(r_temp - r_3D, 'fro'); 
        disp(frobenius_norm_delta_r); 

    end 
end 