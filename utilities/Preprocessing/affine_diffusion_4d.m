function [smoothed_array] = affine_diffusion_4d(original_array,t_tot,dt,aff_flag,timer,padding)
% Created by Xinan Chen on 12/06/2022 based on affine_diffusion_3d.m
% affine_diffusion_4d smooth an image using the affine-invariant
% mean_curvature flow.
% Modified by Kaiming Xu: replace loops by matrix calculation to reduce
% computing time.
%
%   t_tot    = total evolution time, in "seconds"
%   dt       = numerical evolution parameter
%   aff_flag = 1 if affine-invariant flow, 0 for linear smoothing
%   timer  = 1 if step counter and total runtime should be printed
%       (probably leave these flags both 1)
%   padding = 'zeros' or 'replicate' or 'none' (padding in 4th dimension, the time dimension, by one voxel)
%
% Example usage:
% C = affine_diffusion_4d(A, .5, 0.02, 1, 1);

arguments
    original_array 
    t_tot = 1
    dt = 0.1
    aff_flag = 1
    timer = 1
    padding = 'replicate'
end

if nargin < 6
    fprintf('In affine_diffusion_4d.m: No padding option was loaded. Then use default padding = %s\n',padding);
end

phi = double(original_array);    % cast to double to avoid rounding

switch padding
    case 'none'

    case 'zeros'
        [N1,N2,N3,N4] = size(phi);
        tmp = zeros(N1,N2,N3,N4+2);
        tmp(:,:,:,2:N4+1) = phi;
        phi = tmp;
    case 'replicate'
        [N1,N2,N3,N4] = size(phi);
        tmp = zeros(N1,N2,N3,N4+2);
        tmp(:,:,:,2:N4+1) = phi;
        tmp(:,:,:,1) = phi(:,:,:,1);
        tmp(:,:,:,end) = phi(:,:,:,end);
        phi = tmp;
    otherwise
        error('In affine_diffusion_4d.m: padding option is only limited to none, zeros or replicate.');
end

%% Set up parameters
if ndims(phi) ~= 4 
    error('Error: takes a 4D array.')
end
[Nx,Ny,Nz,Nt] = size(phi);


n_t = round(t_tot/dt);
%% Evolve the image
if timer
    tic
end

for t = 1:n_t

    % first order derivatives
    I1 = zeros(size(phi)); I2 = zeros(size(phi)); I3 = zeros(size(phi)); I4 = zeros(size(phi)); 
    % second order derivatives
    I11 = zeros(size(phi)); I22 = zeros(size(phi)); I33 = zeros(size(phi)); I44 = zeros(size(phi));
    if(aff_flag)
        I12 = zeros(size(phi)); I13 = zeros(size(phi)); I14 = zeros(size(phi));
        I23 = zeros(size(phi)); I24 = zeros(size(phi)); I34 = zeros(size(phi));
    end
    
    % Compute derivatives for all voxels
    I1(2:end-1,2:end-1,2:end-1,2:end-1) = 0.5*(phi(3:end,2:end-1,2:end-1,2:end-1) - phi(1:end-2,2:end-1,2:end-1,2:end-1));
    I2(2:end-1,2:end-1,2:end-1,2:end-1) = 0.5*(phi(2:end-1,3:end,2:end-1,2:end-1) - phi(2:end-1,1:end-2,2:end-1,2:end-1));
    I3(2:end-1,2:end-1,2:end-1,2:end-1) = 0.5*(phi(2:end-1,2:end-1,3:end,2:end-1) - phi(2:end-1,2:end-1,1:end-2,2:end-1));
    I4(2:end-1,2:end-1,2:end-1,2:end-1) = 0.5*(phi(2:end-1,2:end-1,2:end-1,3:end) - phi(2:end-1,2:end-1,2:end-1,1:end-2));

    I11(2:end-1,2:end-1,2:end-1,2:end-1) = phi(3:end,2:end-1,2:end-1,2:end-1) - 2*phi(2:end-1,2:end-1,2:end-1,2:end-1) + phi(1:end-2,2:end-1,2:end-1,2:end-1);
    I22(2:end-1,2:end-1,2:end-1,2:end-1) = phi(2:end-1,3:end,2:end-1,2:end-1) - 2*phi(2:end-1,2:end-1,2:end-1,2:end-1) + phi(2:end-1,1:end-2,2:end-1,2:end-1);
    I33(2:end-1,2:end-1,2:end-1,2:end-1) = phi(2:end-1,2:end-1,3:end,2:end-1) - 2*phi(2:end-1,2:end-1,2:end-1,2:end-1) + phi(2:end-1,2:end-1,1:end-2,2:end-1);
    I44(2:end-1,2:end-1,2:end-1,2:end-1) = phi(2:end-1,2:end-1,2:end-1,3:end) - 2*phi(2:end-1,2:end-1,2:end-1,2:end-1) + phi(2:end-1,2:end-1,2:end-1,1:end-2);
    if aff_flag
        I12(2:end-1,2:end-1,2:end-1,2:end-1) = 0.25*(phi(3:end,3:end,2:end-1,2:end-1) + phi(1:end-2,1:end-2,2:end-1,2:end-1) - phi(3:end,1:end-2,2:end-1,2:end-1) - phi(1:end-2,3:end,2:end-1,2:end-1));
        I23(2:end-1,2:end-1,2:end-1,2:end-1) = 0.25*(phi(2:end-1,3:end,3:end,2:end-1) + phi(2:end-1,1:end-2,1:end-2,2:end-1) - phi(2:end-1,3:end,1:end-2,2:end-1) - phi(2:end-1,1:end-2,3:end,2:end-1));
        I13(2:end-1,2:end-1,2:end-1,2:end-1) = 0.25*(phi(3:end,2:end-1,3:end,2:end-1) + phi(1:end-2,2:end-1,1:end-2,2:end-1) - phi(3:end,2:end-1,1:end-2,2:end-1) - phi(1:end-2,2:end-1,3:end,2:end-1));
        I24(2:end-1,2:end-1,2:end-1,2:end-1) = 0.25*(phi(2:end-1,3:end,2:end-1,3:end) + phi(2:end-1,1:end-2,2:end-1,1:end-2) - phi(2:end-1,3:end,2:end-1,1:end-2) - phi(2:end-1,1:end-2,2:end-1,3:end));
        I14(2:end-1,2:end-1,2:end-1,2:end-1) = 0.25*(phi(3:end,2:end-1,2:end-1,3:end) + phi(1:end-2,2:end-1,2:end-1,1:end-2) - phi(3:end,2:end-1,2:end-1,1:end-2) - phi(1:end-2,2:end-1,2:end-1,3:end));
        I34(2:end-1,2:end-1,2:end-1,2:end-1) = 0.25*(phi(2:end-1,2:end-1,3:end,3:end) + phi(2:end-1,2:end-1,1:end-2,1:end-2) - phi(2:end-1,2:end-1,3:end,1:end-2) - phi(2:end-1,2:end-1,1:end-2,3:end));
    end

%     if t==1
% %         figure, montageArray(squeeze(J(1,:,:,:,:))), axis image, title('first order derivative, x axis')
% %         figure, montageArray(squeeze(J(2,:,:,:,:))), axis image, title('first order derivative, y axis')
% %         figure, montageArray(squeeze(J(3,:,:,:,:))), axis image, title('first order derivative, z axis')
% %         figure, montageArray(squeeze(J(4,:,:,:,:))), axis image, title('first order derivative, t axis')
% 
%         figure, montageArray(squeeze(I(1,1,:,:,:,:))), axis image, title('second order derivative, 11')
%         figure, montageArray(squeeze(I(2,1,:,:,:,:))), axis image, title('second order derivative, 21')
%         figure, montageArray(squeeze(I(1,2,:,:,:,:))), axis image, title('second order derivative, 12')
%         figure, montageArray(squeeze(I(2,2,:,:,:,:))), axis image, title('second order derivative, 22')
%     end
    
    if aff_flag
        % Compute curvature quantities
        
        meanCurvNum_part1 = (I1.^2+I2.^2+I3.^2+I4.^2).*(I11+I22+I33+I44);
        
        meanCurvNum_part2 = I1.^2.*I11 + I2.^2.*I22 + I3.^2.*I33 + I4.^2.*I44 + ...
            2*(I1.*I2.*I12 + I1.*I3.*I13 + I1.*I4.*I14 + I3.*I2.*I23 + I4.*I2.*I24 + I3.*I4.*I34);
        meanCurvNum = meanCurvNum_part1- meanCurvNum_part2;

        gausCurvNum = -(I1.^2.*(-I22.*I34.^2 - I33.*I24.^2 - I44.*I23.^2 + I22.*I33.*I44 + 2*I23.*I24.*I34)...
            + I2.^2.*(-I11.*I34.^2 - I33.*I14.^2 - I44.*I13.^2 + I11.*I33.*I44 + 2*I13.*I14.*I34)...
            + I3.^2.*(-I11.*I24.^2 - I22.*I14.^2 - I44.*I12.^2 + I11.*I22.*I44 + 2*I12.*I14.*I24)...
            + I4.^2.*(-I11.*I23.^2 - I22.*I13.^2 - I33.*I12.^2 + I11.*I22.*I33 + 2*I12.*I13.*I23)...
            - 2*I1.*I2.*(I12.*(I33.*I44 - I34.^2) - I33.*I14.*I24 - I44.*I13.*I23 + I34.*(I13.*I24+I14.*I23))...
            - 2*I1.*I3.*(I13.*(I22.*I44 - I24.^2) - I22.*I14.*I34 - I44.*I12.*I23 + I24.*(I12.*I34+I14.*I23))...
            - 2*I1.*I4.*(I14.*(I22.*I33 - I23.^2) - I22.*I13.*I34 - I33.*I12.*I24 + I23.*(I12.*I34+I13.*I24))...
            - 2*I2.*I3.*(I23.*(I11.*I44 - I14.^2) - I11.*I24.*I34 - I44.*I12.*I13 + I14.*(I12.*I34+I13.*I24))...
            - 2*I2.*I4.*(I24.*(I11.*I33 - I13.^2) - I11.*I23.*I34 - I33.*I12.*I14 + I13.*(I12.*I34+I14.*I23))...
            - 2*I3.*I4.*(I34.*(I11.*I22 - I12.^2) - I11.*I23.*I24 - I22.*I13.*I14 + I12.*(I13.*I24+I14.*I23)));

        % Update phi according to mean curvature flow
        %phi = phi + dt*sign(meanCurvNum).*nthroot(gausCurvNum,5);
        phi = phi + dt*sign(meanCurvNum).*nthroot(max(0,gausCurvNum),5);
        %phi = max(zeros(size(phi)),phi + dt*sign(meanCurvNum).*nthroot(gausCurvNum,5));

%         if t==1
%             figure, montageArray(gausCurvNum), axis image, title('gausCurvNum 4d')
%             %figure, histogram(gausCurvNum)
%             figure, montageArray(meanCurvNum), axis image, title('meanCurvNum 4d')
%             %figure, histogram(meanCurvNum)
%             %figure, montageArray(squeeze(I(1,3,:,:,:,:))), axis image, title('I(1,3)')
%         end

    else % just do linear
        phi = phi + dt*(I11+I22+I33+I44);
    end
    
    if timer
        fprintf('Step %d out of %d.\n', t,n_t);
    end
end

if timer
    rt = toc;
    fprintf('\n Evolution took %f seconds.\n', rt);
end

% return
%smoothed_array = phi;
switch padding
    case 'none'
        smoothed_array = phi;
    case {'replicate','zeros'}
        smoothed_array = phi(:,:,:,2:Nt-1);
end

end
