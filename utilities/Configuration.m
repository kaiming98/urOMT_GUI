classdef Configuration < handle & matlab.mixin.CustomDisplay
    %CONFIGURATION Superclass that contains the MRI data images and the
    %functions necessary to run the regularized optimal mass transport
    %(rOMT) analysis. A separate superclass and script are used for
    %post-processing and visualization. See "Postprocessing.m"
    %

    properties
        densityTarget % true final density (from next image in OMT loop)
        alpha = 0 % penalty constant for source term
        alphaFunctionalType = 'constant'
        pathAlpha
        pathAnatomyImage
        mask Mask
    end

    properties (SetAccess = protected)
        B % placeholder for a matrix used in AdvecDiffusion calculation
        Xc
        Yc
        Grad % gradient operator matrix
        hd % volume of one voxel
        Data %struct containing 3D image data
        A % A and Abig are helper matrixes used in GaussNewton
        Abig
        nFiles {mustBeInteger}
        pathOutput
    end

    properties (SetAccess = immutable)
        n % domain size
        dimension
        pathMask
        pathData
        pathOutputFolder = './output'
        isDataPacked = 0
        startpointsMaskName
        pathStartpointsMask
        pathMaxDensity
        extension % dataset type
        xRange
        yRange
        zRange
        boundaryConditions
    end

    properties (SetAccess = private)
        reinitializeVelocityFlag = 1
        reinitializeSourceFlag = 1
        %from JSON configuration file (constants)
        tag string = ''
        dataset_name
        mask_number = 1 % identifier for which mask
        do_resize = 0
        sizeFactor = 1
        smooth = 1
        preSmooth = 'None'
        reinitializeDensityFlag % reinitiliaze input data
        dilate = 0
        data_start
        data_end
        boundaryConditionOpen = 1
        maxIterationsVelocity = 6
        timeInitial 
        timeJump = 2
        timeFinal
        sigma = 0.002
        convert2concentration = 0; %do the concentration images need to be converted to signal (0 no, 1 yes)
        isMaskFilled = 1; %1 mask is filled; 0 mask is not filled and needs further filling.
        dt = 0.4% numerical time difference between each adjacent time steps.
        nt  = 10 % sets temporal discretization: # of numerical time steps
        % between source image and target image in one loop.
        gamma = 80
        beta  = 10000 %penalty constant for final frame matching error
        maxIterations_pcg % iterations to run PCG algorithm
        domainSize
        trueSize
        rawSize
        addSource = 0
        diffusionCoefficientType = 'constant'
        sigmaFunctionType = 1
        sigmaFunctionConstant = 1
        KList
    end

    methods

        function obj = Configuration(tag)
            %Configuration Construct an instance of this class
            %   Detailed explanation goes here
            % load json struct fields as objects
            arguments
                tag string
            end

            pathConfiguration = sprintf("%s.json",tag);
            obj.tag = tag;
            if exist(pathConfiguration,'file')
                json = jsondecode(fileread(pathConfiguration));
            else
                error("missing or unrecognized configuration file ")
            end

            for fn = fieldnames(json)' %enumerate fields

                try
                    obj.(fn{1}) = json.(fn{1}); %and copy
                catch
                    warning('Could not copy field %s', fn{1});
                end

            end
            obj.timeFinal = obj.timeFinal - obj.timeJump;
            % obj.do_resize = 1;
            % obj.sizeFactor = 0.5;
            %assert(strcmp(obj.tag,tag),"configuration file name and tag do not match");
            if ~isempty(obj.pathMask)
                [obj.xRange, obj.yRange, obj.zRange, obj.rawSize] = find_range(obj.pathMask,obj.dilate);
                [~, ~, obj.extension] = fileparts(obj.pathMask);
            else
                [~, ~, obj.extension] = fileparts(obj.pathData);
                if obj.isDataPacked
                    rawImage = load_untouch_nii(obj.pathData);

                else
                    rawImage = nii2mat(sprintf('%s%02d%s', obj.pathData, ...
                        timeInitial, '.nii'));
                end
                obj.xRange = 1:rawImage.hdr.dime.dim(2);
                obj.yRange = 1:rawImage.hdr.dime.dim(3);
                obj.zRange = 1:rawImage.hdr.dime.dim(4);
                obj.rawSize = rawImage.hdr.dime.dim(2:4);                  
            end
            obj.domainSize = [length(obj.xRange), length(obj.yRange), length(obj.zRange)];

            obj.trueSize = round(obj.sizeFactor * obj.domainSize);
            
            %outputPath containing frames chosen
            if strcmp(obj.preSmooth,'None')
                pathOutput = sprintf('%s/%s_E%d_E%d_addSource_%d',obj.pathOutputFolder,obj.tag,obj.timeInitial,obj.timeFinal+obj.timeJump,obj.addSource);
            else
                pathOutput = sprintf('%s/%s_E%d_E%d_presmooth_%s_addSource_%d',obj.pathOutputFolder,obj.tag,obj.timeInitial,obj.timeFinal+obj.timeJump,obj.preSmooth,obj.addSource);
            end
            if obj.beta~=10000
                pathOutput = sprintf('%s_beta_%d',pathOutput,obj.beta);
            end
            
            if obj.addSource
                switch obj.alphaFunctionalType
                    case 'constant'
                        pathOutput = strcat(pathOutput,sprintf('_a=%0.3f',obj.alpha));
                    case 'auto'
                        pathOutput = strcat(pathOutput,'_a=auto');
                    case 'manual'
                        pathOutput = strcat(pathOutput,'_',obj.pathAlpha(1:end-4));
                end
            end
            if obj.mask_number~=1
                pathOutput = strcat(pathOutput,sprintf('_mask%d',obj.mask_number));
            end
            if obj.dt~=0.4
                pathOutput = strcat(pathOutput,sprintf('_dt=%0.1f',obj.dt));
            end
            if obj.sigma~=0.002
                pathOutput = strcat(pathOutput,sprintf('_sigma=%0.5f',obj.sigma));
            end
            if strcmp(obj.diffusionCoefficientType,'anisotropic')
                pathOutput = strcat(pathOutput,sprintf('_%s_sType=%d_K=%d',obj.diffusionCoefficientType,obj.sigmaFunctionType,obj.sigmaFunctionConstant));
            end
            if strcmp(obj.diffusionCoefficientType,'autoAnisotropic')
                pathOutput = strcat(pathOutput,sprintf('_%s_sType=%d',obj.diffusionCoefficientType,obj.sigmaFunctionType));
            end
            if strcmp(obj.diffusionCoefficientType,'brightness')
                pathOutput = strcat(pathOutput,sprintf('_%s_sType=%d_K=%d',obj.diffusionCoefficientType,obj.sigmaFunctionType,obj.sigmaFunctionConstant));
            end
            obj.pathOutput = pathOutput;           
            
            if obj.boundaryConditionOpen
                obj.boundaryConditions = 'open';
            else
                obj.boundaryConditions = 'closed';
            end
            
            assert(isinteger(int8((obj.timeFinal-obj.timeInitial)/obj.timeJump)),...
                "nFiles must be an integer. Adjust timeInitial, timeFinal and timeJump in configuration .json file");
            obj.nFiles = length(obj.timeInitial:obj.timeJump:obj.timeFinal);

            obj.n = obj.trueSize';
            obj.dimension = length(obj.n);

            % Preallocate data struct
            obj.Data = struct('volume', cell(1, obj.nFiles));

        end

        function loadVolume(obj,timeInitial,timeJump,timeFinal)
            if nargin == 1
                timeInitial = obj.timeInitial;
                timeJump = obj.timeJump;
                timeFinal = obj.timeFinal;
            end
            % load volumetric data
            if isempty(obj.pathData)
                error("need path to data")
            end

            if isempty(obj.mask) 
                obj.mask = Mask(obj.pathMask,obj.isMaskFilled,obj.xRange,obj.yRange,obj.zRange);
                if obj.do_resize
                    obj.mask = obj.mask.resize(obj.sizeFactor);
                end

                if obj.dilate>0
                    obj.mask = obj.mask.dilate(obj.dilate);
                end
            end
            
            % load data 
            if strcmp(obj.extension,'.gz') || obj.isDataPacked % packed 

                X = nii2mat(obj.pathData, obj.xRange, obj.yRange, obj.zRange);

                if obj.sizeFactor && ~(obj.sizeFactor == 1)
                    X = resizeMatrix(X, round(obj.sizeFactor .* size(X)), 'linear');
                end

                for kk = timeInitial:timeJump:(timeFinal+timeJump)
                    volume = X(:,:,:,kk);
                    if obj.smooth > 0
                        volume = affine_diffusion_3d(volume, obj.smooth, 0.1, 1);
                    end
                    volume(~obj.mask.contents) = 0;
                    % calculate contrast concentration from image intensity
                    if obj.convert2concentration
                        volume = obj.sig2conc(volume);
                    end
                    obj.Data((kk-timeInitial)/timeJump+1).volume = volume;
                end

            elseif strcmp(obj.extension,'.nii') && ~obj.isDataPacked % unpacked
                nFiles = length(timeInitial:timeJump:timeFinal)+1;
                for kk = 1:nFiles
                    X = nii2mat(sprintf('%s%02d%s', obj.pathData, ...
                        timeInitial + (kk - 1) * timeJump, '.nii'), ...
                        obj.xRange, obj.yRange, obj.zRange);
    
                    if obj.sizeFactor && ~(obj.sizeFactor == 1)
                        X = resizeMatrix(X, round(obj.sizeFactor .* size(X)), 'linear');
                    end
    
                    if obj.smooth > 0
                        X = affine_diffusion_3d(X, obj.smooth, 0.1, 1);
                    end
    
                    X(~obj.mask.contents) = 0;
                    obj.Data(kk).volume = X;
                end
            else
                error("Data must be packed ('*nii.gz') or unpacked ('*.nii') NIfTI files")
            end

        end

        function saveInitialDensity(obj)
            
            % create an output directory
            if ~exist(obj.pathOutput, 'dir')
                mkdir(obj.pathOutput)
            end

            if isempty(obj.Data(1).volume)
                error("nothing to save. object has not been initialized wiht data")
            end

            fileName = sprintf('%s/rho_%s_%d_t_0.txt', obj.pathOutput, obj.tag, obj.timeInitial);

            if exist(fileName, 'file')
                delete(fileName);
            end
            writematrix(obj.Data(1).volume(:), fileName);
                

        end

        function initializeGrid(obj)

            obj.checkDimension(obj.dimension);

            obj.getCellCenteredData;

            % Precompute operators for later use in the gradient descent
            % for the velocity (GNblock_u.m)

            % calculate B, A, Abig
            DiffusionOperator = -obj.sigma * (obj.Grad') * obj.Grad; % sigma * laplacian
            I = speye(prod(obj.n));
            obj.B = I - obj.dt * DiffusionOperator;

            obj.A = kron(ones(1, obj.dimension), speye(prod(obj.n)));
            obj.Abig = kron(speye(obj.nt), obj.A);

        end

        function velocity = getVelocity(obj, Rho_i, velocity)
            %%GETVELOCITY
            % Finds optimal velocity using Gauss Newton optimization
            % algorithm.

            arguments
                obj Configuration
                Rho_i (:, 1) double % initial (source) density
                velocity double % velocity
            end

            % calculates initial phi (energy functional, Gamma or phi)
            phi = obj.getPhi(Rho_i, velocity);

            if ~strcmp(obj.diffusionCoefficientType,'constant')
                % calculate derivative of diffusion term w.r.t initial density
                dD{1} = obj.calculateDiffusionJacobian(Rho_i);
            end

            for i = 1:obj.maxIterationsVelocity 

                Rho = obj.advectionDiffusion(Rho_i, velocity);

                RHO0 = [Rho_i, Rho(:, 1:end - 1)];

                U = reshape(velocity, obj.dimension * prod(obj.n), []);

                M = obj.calculateM(RHO0, U);
                if ~strcmp(obj.diffusionCoefficientType,'constant')
                    %calculate the derivative of diffusion term w.r.t density
                    for iTime = 2:obj.nt
                        dD{iTime} = obj.calculateDiffusionJacobian(RHO0(:,iTime));
                    end
                    M.dD = dD;
                end
                %% Calculate derivatives of MK functional
        
                dMongeKantorovichFunctional = zeros(obj.dimension * prod(obj.n), obj.nt);

                for j = 1:obj.nt
                    dMongeKantorovichFunctional(:, 1:j) = dMongeKantorovichFunctional(:, 1:j) + reshape(obj.dt * ...
                        get_drNduT(M, j, obj.dt, obj, obj.A * (U(:, j) .* U(:, j))),...
                        obj.dimension * prod(obj.n), j);
                end
                
                %% Calculate gradients (g) of total functional 

                g = (2 * obj.dt * Rho(:)' * obj.Abig * sdiag(velocity(:)))' ...
                    + dMongeKantorovichFunctional(:) ...
                    + obj.beta * get_drNduT(M, obj.nt, obj.dt, obj, Rho(:, end) - obj.densityTarget) ...
                    + get_dRudu(velocity, obj.nt, obj)' * obj.gamma * obj.dt;

                fprintf('%3d.%d\t      %3.5e \t     ||g|| = %3.5e\n', i, 0, phi, norm(g));

                %% Calculate the Hessian (H) of total functional
                H = @(x) 2 * obj.dt * sdiag(Rho(:)' * obj.Abig) * x + ...
                obj.gamma * obj.dt .* kron(speye(obj.nt * obj.dimension), (obj.Grad)' * obj.Grad) * x + ...
                    obj.beta * get_drNduTdrNdu(M, obj.nt, obj.dt, obj, x);

                %% Gauss-Newton optimization algorithm using pcg for  
                % Solving for s in the linear equation Hs = -g 

                [s, pcgflag, relres, iter] = pcg(H, -g, 0.01, obj.maxIterations_pcg);

                if pcgflag ~= 0
                    fprintf(['\nUPDATE: pcg output is approximated in case %s (pcgflag = %d).\n' ...
                        'Relative residual = %3.2e after %d/%d pcg iterations' ...
                        ' and %d/%d getVelocity main loop iterations.\n\n']...
                        ,obj.tag, pcgflag, relres, iter, obj.maxIterations_pcg, i, obj.maxIterationsVelocity);
                end

                %% line search
                % Find the appropriate length of the gradient descent step. 

                lengthLineSearch = 0.7; % initialization (typically 1)

                iLineSearch = 1;

                while 1
                    velocityTemp = velocity(:) + lengthLineSearch * s;

                    phiTemp = obj.getPhi(Rho_i, velocityTemp);

                    fprintf('%3d.%d\t      %3.5e \t     phiTemp  = %3.5e        %s\n',...
                        i, iLineSearch, phi, phiTemp, obj.tag);
                    
                    phiThreshold = phi + 1e-8 * s' * g;
                    
                    uApproachingZero = norm(velocity(:))<1e-5;


                    if phiTemp < phiThreshold
                        break;
                    end

                    lengthLineSearch = lengthLineSearch / 2;
                    iLineSearch = iLineSearch + 1;

                    if iLineSearch > 20
                        %velocity = u;
                        fprintf('LineSearch breaks');

                        if i<= 2 
                            fprintf('\n image pair may be too similar\n');
                        else
                            fprintf('\n');
                        end

                        return;

                    end

                end

                velocity = velocityTemp;
                phi = phiTemp;

            end

        end
        
        function [velocity,source] = getVelocityAndSource(obj, Rho_i, velocity, source)
            %%GETVELOCITYANDSOURCE
            % Finds optimal velocity and source using Gauss Newton optimization
            % algorithm.

            arguments
                obj Configuration
                Rho_i (:, 1) double % initial density
                velocity double % velocity
                source double % source
            end

            % calculates initial phi (energy functional, Gamma or phi)
            phi = obj.getPhi(Rho_i, velocity, source);
            
            if ~strcmp(obj.diffusionCoefficientType,'constant')
                % calculate derivative of diffusion term w.r.t initial density
                dD{1} = obj.calculateDiffusionJacobian(Rho_i);
            end




            for i = 1:obj.maxIterationsVelocity 

                Rho = obj.advectionDiffusionSource(Rho_i, velocity, source);

                RHO0 = [Rho_i, Rho(:, 1:end - 1)];

                U = reshape(velocity, obj.dimension * prod(obj.n), []);
                R = reshape(source, prod(obj.n), []);

                M = obj.calculateM(RHO0, U, R);
                
                if ~strcmp(obj.diffusionCoefficientType,'constant')
                    %calculate the derivative of diffusion term w.r.t density
                    for iTime = 2:obj.nt
                        dD{iTime} = obj.calculateDiffusionJacobian(RHO0(:,iTime));
                    end
                    M.dD = dD;
                end
                %% Calculate gradients (g) of total functional 
                dMatchingTerm_u = 2 * obj.beta * get_drNduT(M,obj.nt,obj.dt,obj,Rho(:, end) - obj.densityTarget);
                dMatchingTerm_r = 2 * obj.beta * get_drNdrT(M,obj.nt,obj.dt,obj,Rho(:, end) - obj.densityTarget);
                
                dKineticEnergy_u1 = 2 * (Rho(:)' * obj.Abig * sdiag(velocity(:)))';
                dSourceEnergy_r1 = 2 * obj.alpha * (Rho(:)' * sdiag(source(:)))';
                              
                dKineticEnergy_u2 = zeros(obj.dimension * prod(obj.n), obj.nt);
                dSourceEnergy_u = zeros(obj.dimension * prod(obj.n), obj.nt);
                dKineticEnergy_r = zeros(prod(obj.n), obj.nt);
                dSourceEnergy_r2 = zeros(prod(obj.n), obj.nt);

                for j = 1:obj.nt
                    dKineticEnergy_u2(:, 1:j) = dKineticEnergy_u2(:, 1:j) + ...
                        reshape(get_drNduT(M, j, obj.dt, obj, obj.A * (U(:, j) .* U(:, j))),...
                        obj.dimension * prod(obj.n), j);
                    dSourceEnergy_u(:, 1:j) = dSourceEnergy_u(:, 1:j) + ...
                        reshape(get_drNduT(M, j, obj.dt, obj, R(:, j) .* R(:, j)),...
                        obj.dimension * prod(obj.n), j);
                    dKineticEnergy_r(:, 1:j) = dKineticEnergy_r(:, 1:j) + ...
                        reshape(get_drNdrT(M, j, obj.dt, obj, obj.A * (U(:, j) .* U(:, j))),...
                        prod(obj.n), j);
                    dSourceEnergy_r2(:, 1:j) = dSourceEnergy_r2(:, 1:j) + ...
                        reshape(get_drNdrT(M, j, obj.dt, obj, R(:, j) .* R(:, j)),...
                        prod(obj.n), j);
                end
                dKineticEnergy_u2 = dKineticEnergy_u2(:);
                dSourceEnergy_u = obj.alpha * dSourceEnergy_u(:);
                dKineticEnergy_r = dKineticEnergy_r(:);
                dSourceEnergy_r2 = obj.alpha * dSourceEnergy_r2(:);
                
                g = [obj.dt * (dKineticEnergy_u1 + dKineticEnergy_u2 + dSourceEnergy_u) + dMatchingTerm_u; ...
                    obj.dt * (dKineticEnergy_r + dSourceEnergy_r2 + dSourceEnergy_r1) + dMatchingTerm_r];

                fprintf('%3d.%d\t      %3.5e \t     ||g|| = %3.5e\n', i, 0, phi, norm(g));
    
                %% Calculate the Hessian (H) of total functional
                
                H = @(x) [2*obj.dt*sdiag(Rho(:)'*obj.Abig)*x(1:obj.dimension*obj.nt*prod(obj.n)) + ...
                    2*obj.beta*get_drNduTdrNdu(M,obj.nt,obj.dt,obj,x(1:obj.dimension*obj.nt*prod(obj.n))) + ...
                    2*obj.beta*get_drNduTdrNdr(M,obj.nt,obj.dt,obj,x(obj.dimension*obj.nt*prod(obj.n)+1:end)); ...
                    2*obj.beta*get_drNdrTdrNdu(M,obj.nt,obj.dt,obj,x(1:obj.dimension*obj.nt*prod(obj.n))) + ...
                    2*obj.beta*get_drNdrTdrNdr(M,obj.nt,obj.dt,obj,x(obj.dimension*obj.nt*prod(obj.n)+1:end)) + ...
                    2*obj.alpha*obj.dt*sdiag(Rho(:))*x(obj.dimension*obj.nt*prod(obj.n)+1:end)];
                
                %% Gauss-Newton optimization algorithm using pcg for  
                % Solving for s in the linear equation Hs = -g 

                [s, pcgflag, relres, iter] = pcg(H, -g, 0.01, obj.maxIterations_pcg);

                if pcgflag ~= 0
                    fprintf(['\nUPDATE: pcg output is approximated in case %s (pcgflag = %d).\n' ...
                        'Relative residual = %3.2e after %d/%d pcg iterations' ...
                        ' and %d/%d getVelocity main loop iterations.\n\n']...
                        ,obj.tag, pcgflag, relres, iter, obj.maxIterations_pcg, i, obj.maxIterationsVelocity);
                end

                %% line search
                % Find the appropriate length of the gradient descent step. 

                lengthLineSearch = 0.7; % initialization (typically 1)

                iLineSearch = 1;

                while 1
                    velocityTemp = velocity(:) + lengthLineSearch * s(1:obj.dimension*prod(obj.n)*obj.nt);
                    sourceTemp = source(:) + lengthLineSearch * s(obj.dimension*prod(obj.n)*obj.nt+1:end);

                    phiTemp = obj.getPhi(Rho_i, velocityTemp, sourceTemp);

                    fprintf('%3d.%d\t      %3.5e \t     phiTemp  = %3.5e        %s\n',...
                        i, iLineSearch, phi, phiTemp, obj.tag);
                    
                    phiThreshold = phi + 1e-8 * s' * g;

                    if phiTemp < phiThreshold
                        break;
                    end

                    lengthLineSearch = lengthLineSearch / 2;
                    iLineSearch = iLineSearch + 1;

                    if iLineSearch > 4
                        
                        fprintf('LineSearch breaks');

                        if i<= 2 
                            fprintf('\n image pair may be too similar\n');
                        else
                            fprintf('\n');
                        end

                        return;

                    end

                end

                velocity = velocityTemp;
                source = sourceTemp;
                phi = phiTemp;

            end

        end

        function [phi, mk, phiN, rho, Ru, sourceEnergy] = getPhi(obj, rho0, u, r)
            %% GETPHI Calculates energy functional
            %   rho0 - initial image data (flattened vector of size n1*n2*n3*1)
            %   u    - optimal velocity field
            %   r    - source term
            
            if obj.addSource && nargin<4
                error("missing source term");
            end
            if obj.addSource && isempty(obj.alpha)
                error("missing source term penalty constant");
            end
            if obj.addSource
                rho = obj.advectionDiffusionSource(rho0, u, r);
                r = reshape(r,prod(obj.n),obj.nt);
                sourceEnergy = obj.dt*sum(rho.*(r.*r),'all');
            else
                rho = obj.advectionDiffusion(rho0, u);
                sourceEnergy = 0;
            end

            mk = obj.getMongeKantorovichFunctional(rho, u); %=hd*dt*rho'*||v||^2
            
            %calculating matching term in the energy functional. (second term)
            phiN = 0.5 * norm(rho(:, end) - obj.densityTarget)^2;

            %% smoothing deformation field
            Ru = 0;

            if obj.gamma ~= 0
                uvec = reshape(u(:), [], obj.dimension * obj.nt);

                for ii = 1:obj.dimension * obj.nt
                    Ru = Ru + 0.5 * obj.dt * (norm(obj.Grad * uvec(:, ii))^2);
                end

            end

            phi = mk + obj.beta * phiN + obj.gamma * Ru + obj.alpha * sourceEnergy;
        end

        function saveRecord(obj, matrixToSave)

            % create an output directory
            if ~exist(obj.pathOutput, 'dir')
                mkdir(obj.pathOutput)
            end

            % create record file
            recordPath = sprintf('%s/record.txt', obj.pathOutput);

            if ~exist(recordPath, 'file')
                if ~obj.addSource
                    obj.createRecordFile(recordPath, [0 0 0 0 0 0 0 0 0], ...
                        {'time-ind', 'ti', 'tf', ...
                            'phi', 'mk', 'Ru', ...
                            'phiN', 'max(u)', 'toc'});
                else
                    obj.createRecordFile(recordPath, [0 0 0 0 0 0 0 0 0 0 0], ...
                        {'time-ind', 'ti', 'tf', ...
                            'phi', 'mk', 'Ru', ...
                            'phiN', 'se','max(u)', 'max(r)', 'toc'});
                end
            end

            dlmwrite(recordPath, matrixToSave, '-append');

        end

        function saveVariable(obj, timeOrigin, timeTarget, iTime, variable)
            %SAVE Saves the variable calculated following a rOMT loop
            %   timeOrigin  - the time point of the source image
            %   timeTarget  - the time point of the target image
            %   iTime       - index counter between first and last time points
            %   variable    - what will be saved

            arguments
                obj
                timeOrigin (1, 1) double
                timeTarget (1, 1) double
                iTime (1, 1) double
                variable double
            end

            variableName = inputname(5);
            save(sprintf('%s/%s_%s_%d_%d_t_%d.mat', obj.pathOutput, obj.tag, variableName(1), ...
                timeOrigin, timeTarget, iTime), 'variable');

        end
        
        function rho = advectionDiffusion(obj, rho0, u)
            
            if isempty(obj.B)
                obj.initializeGrid;
            end
            
            rho = zeros(prod(obj.n), obj.nt + 1);
            rho(:, 1) = rho0;

            u = reshape(u, 2 * prod(obj.n), obj.nt);
            for i = 2:obj.nt + 1

                U1 = reshape(u(1:prod(obj.n), i - 1), obj.n');
                U2 = reshape(u(prod(obj.n) + 1:end, i - 1), obj.n');

                S = dTrilinears2d(rho(:, i - 1), obj.Xc + obj.dt * U1, obj.Yc + obj.dt * U2, ...
                    1, 1, obj.boundaryConditions);
                
                switch obj.diffusionCoefficientType
                    case 'constant'
                        rho(:, i) = S * rho(:, i - 1);
                        [rho(:, i), pcgflag] = pcg(obj.B, rho(:, i));
                        if pcgflag ~= 0
                            warning('MATLAB:pcgExitFlag', 'Warning: advecDiff.m >>> while finding rho(:,%d), pcg exit flag = %d', i, pcgflag)
                        end
                    case 'anisotropic'
                        rho(:, i) = S * rho(:, i - 1) + ...
                            obj.dt*obj.calculateAnisotropicDiffusion(rho(:, i - 1));
                    case 'autoAnisotropic'
                        rho(:, i) = S * rho(:, i - 1) + ...
                            obj.dt*obj.calculateAnisotropicDiffusion(rho(:, i - 1));
                    case 'brightness'
                        diffusionResult = obj.calculateBrightnessDiffusion(rho(:, i - 1));
                        rho(:, i) = S * rho(:, i - 1) + ...
                            obj.dt*diffusionResult; 
                end

            end

            rho = rho(:, 2:end);
        end 
        
        function rho = advectionDiffusionSource(obj, rho0, u, r)
            if isempty(obj.B)
                obj.initializeGrid;
            end
            
            rho = zeros(prod(obj.n), obj.nt + 1);
            rho(:, 1) = rho0;

            u = reshape(u, 2 * prod(obj.n), obj.nt);
            r = reshape(r, prod(obj.n), obj.nt);
            for i = 2:obj.nt + 1
                
                rhoSource = (1 + obj.dt * r(:, i - 1)) .* rho(:, i - 1);
                
                U1 = reshape(u(1:prod(obj.n), i - 1), obj.n');
                U2 = reshape(u(prod(obj.n) + 1:end, i - 1), obj.n');

                S = dTrilinears2d(rhoSource, obj.Xc + obj.dt * U1, obj.Yc + obj.dt * U2, ...
                    1, 1, obj.boundaryConditions);
                
                switch obj.diffusionCoefficientType
                    case 'constant'
                        rho(:, i) = S * rhoSource;
                        [rho(:, i), pcgflag] = pcg(obj.B, rho(:, i));
                        if pcgflag ~= 0
                            warning('MATLAB:pcgExitFlag', 'Warning: advecDiff.m >>> while finding rho(:,%d), pcg exit flag = %d', i, pcgflag)
                        end
                    case 'anisotropic'
                        rho(:, i) = S * rhoSource + ...
                            obj.dt*obj.calculateAnisotropicDiffusion(rho(:, i - 1));
                    case 'brightness'
                        diffusionResult = obj.calculateBrightnessDiffusion(rho(:, i - 1));
                        rho(:, i) = S * rhoSource + ...
                            obj.dt*diffusionResult; 
                end

            end

            rho = rho(:, 2:end);
        end
        
        function matrixOut = sigmaFunctional(obj,matrixIn,K)
            %SIGMAFUNCTIONAL calculate \sigma(x)
            %x could be \rho and |\nabla\rho|
            %case |\nabla\rho|
            %obj.sigmaFunctionType=1 \sigma(x) = \sigma_0/(1+(x/K)^2)
            %obj.sigmaFunctionType=2 \sigma(x) = \sigma_0*exp(-(x/K)^2)
            %To reduce calculation cost, we always assume \sigma(x) is a
            %function of x^2 and matrixIn=x^2.
            %case \rho
            %obj.sigmaFunctionType=0 \sigma(x) = \sigma_0*x/K
            %obj.sigmaFunctionType=1 \sigma(x) = \sigma_0/(1+x/K)
            if nargin==2
                K = obj.sigmaFunctionConstant;
            end
            switch obj.diffusionCoefficientType
                case {'anisotropic', 'autoAnisotropic'}
                    switch obj.sigmaFunctionType
                        case 1
                            matrixOut = obj.sigma./(1+matrixIn/K^2);
                        case 2
                            matrixOut = obj.sigma.*exp(-matrixIn/K^2);
                    end
                case 'brightness'
                    switch obj.sigmaFunctionType
                        case 0
                            matrixOut = obj.sigma.*matrixIn/K;
                        case 1
                            matrixOut = obj.sigma./(1+matrixIn/K);
                    end
            end
        end
        
        function plotVolume(obj,dataInitial,dataFinal)          
            if nargin == 1
                dataInitial = obj.timeInitial;
                dataFinal = obj.timeFinal;
            end
            timeInitial = obj.timeInitial;
            timeFinal = obj.timeFinal;
            timeInitial = timeInitial-dataInitial+1;
            timeFinal = timeFinal-dataInitial+1;
            obj.loadVolume(dataInitial,1,dataFinal-1);
            dataFinal = dataFinal-dataInitial+1;
            dataInitial = 1;
            nFiles = length(dataInitial:dataFinal);
            T = dataInitial:dataFinal;
            Rho = zeros(1,nFiles);
            mask_brain = Mask(obj.pathMask, obj.isMaskFilled, obj.xRange,...
                        obj.yRange, obj.zRange);
            if obj.do_resize
                mask_brain = mask_brain.resize(obj.sizeFactor);
            end
            for i = 1:nFiles
                A = obj.Data(i).volume;
                A(mask_brain.contents == 0) = 0;
                Rho(i) = sum(A,'all');
            end
            N = 10;
            Rhotmp = [Rho(1)*ones(1,N/2),Rho,Rho(end)*ones(1,N/2)];
            B = 1/N*ones(N,1);
            RhoMovingAverage10 = filter(B,1,Rhotmp);
            Rho = RhoMovingAverage10(N+1:end);
            [rho_max, t_max] = max(Rho);
            T_max = t_max + dataInitial-1;
%            T_diff = timeInitial+obj.timeJump/2:obj.timeJump:timeFinal+obj.timeJump/2;
%            Rho_diff = abs(Rho(timeInitial+obj.timeJump:obj.timeJump:timeFinal+obj.timeJump) - ...
%                Rho(timeInitial:obj.timeJump:timeFinal))./Rho(timeInitial:obj.timeJump:timeFinal);
            if obj.addSource
                switch obj.alphaFunctionalType
                    case 'auto'
                        T_alpha = T_diff;
                        alpha = (Rho_diff>0.5)*25000+(Rho_diff<=0.5).*(75000./(4*Rho_diff+1));
                    case 'manual'
                        T_alpha = timeInitial:timeFinal+obj.timeJump;
                        load(obj.pathAlpha);
                    otherwise
                        T_alpha = timeInitial:timeFinal+obj.timeJump;
                        alpha = ones(size(T_alpha))*obj.alpha;
                end
            end
            figure 
            plot(T,Rho, T_max, rho_max, 'ro');
%             hold on
%             p2 = plot(T_diff, Rho_diff*1e6);
            hold on
            plot([T_max, T_max],[0 rho_max],':m')
            xticks(T);
            hold on
            anArrow = annotation('doublearrow') ;
            anArrow.Parent = gca;
            anArrow.Position = [timeInitial, rho_max/2, timeFinal-timeInitial+obj.timeJump, 0];
            line1 = annotation('line') ;
            line1.Parent = gca;
            line1.Position = [timeInitial, 0, 0, Rho(timeInitial)] ;
            line2 = annotation('line') ;
            line2.Parent = gca;
            line2.Position = [timeFinal+obj.timeJump, 0, 0, Rho(timeFinal+obj.timeJump)] ;
            %annotation('doublearrow',[(timeInitial-1)/(timeFinal-1) (timeFinal+1)/(timeFinal-1)],[1/2 1/2])
            %annotation('doublearrow',[0 1],[0 1])
%             yyaxis right
%             p3 = plot(T_alpha,alpha);
%             text(T_alpha(1),alpha(1)+5000,num2str(alpha(1)));
%             text(T_alpha(end),alpha(end)+5000,num2str(alpha(end)));
%             ylim([0 1e5])
%             legend([p2,p3],{'relative absolute difference * 10^6','\alpha'});
%            title(sprintf("%s:total density over time",obj.tag));
            title('Total density over time');
            out_dir = sprintf('./output/%s_%s',obj.dataset_name,obj.preSmooth);
            if ~exist(out_dir, 'dir')
                mkdir(out_dir)
            end
            saveas(gcf,sprintf('%s/frameChose+alpha_%s_%s_mask_%d.png',out_dir,obj.tag,obj.preSmooth,obj.mask_number));
            close
        end
        
        function exportVolumeToVtk(obj,timeInitial,timeFinal)
            if nargin == 1
                timeInitial = obj.timeInitial;
                timeFinal = obj.timeFinal;
            end
            nFiles = length(timeInitial:timeFinal);
            obj.loadVolume(timeInitial,1,timeFinal);
            T = timeInitial:timeFinal;
            mask_brain = Mask(obj.pathMask, obj.isMaskFilled, obj.xRange,...
                        obj.yRange, obj.zRange);
            if obj.do_resize
                mask_brain = mask_brain.resize(obj.sizeFactor);
            end
            out_dir = sprintf('./output/%s', obj.tag);
            if ~exist(out_dir, 'dir')
                mkdir(out_dir)
            end
            for i = T
                A = obj.Data(i-timeInitial+1).volume;
                A(mask_brain.contents == 0) = 0;
                vtkwrite(sprintf('%s/data%d.vtk', out_dir,i), 'structured_points',...
                    'mask', A);
            end
        end
        
        function calculateK(obj,rho)
            rho = reshape(rho,obj.n');
            [g1,g2,g3] = gradient(rho);
            g = sqrt(g1.^2+g2.^2+g3.^2);
            h = histogram(g);
            l = length(h.Values);
            multipleValue = zeros(1,l);
            for i=1:l
                multipleValue(i) = (i-1)*h.Values(i);
            end
            valueSum = cumsum(multipleValue);
            for i=1:l-1
                if valueSum(i)<=0.9*valueSum(l) && valueSum(i+1)>=0.9*valueSum(l)
                    break;
                end
            end
            obj.sigmaFunctionConstant = i;
            if isempty(obj.KList)
                obj.KList = [i];
            else
                obj.KList(end+1)=i;
            end
        end

        function visualizeData(obj)
            x = 1:obj.trueSize(1);
            y = 1:obj.trueSize(2);
            z = 1:obj.trueSize(3);
            obj.loadVolume;
            density = zeros(obj.trueSize(1),obj.trueSize(2),obj.trueSize(3),length(obj.Data));
            for i=1:length(obj.Data)
                density(:,:,:,i) = obj.Data(i).volume;
            end
            density1 = density(density>0);
            maxset = maxk(density1(:),round(0.05*length(density1(:))));
            densityMax = maxset(end);
            for i = obj.timeInitial:obj.timeJump:obj.timeFinal+obj.timeJump
                figure,
                tmp = obj.Data((i-obj.timeInitial)/obj.timeJump+1).volume;
                hs=slice(y,x,z,tmp,y,x,z); 
                set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
                alpha('color'),alphamap(gca,linspace(0,1,100))
                xticks(0:10:obj.trueSize(2)); yticks(0:10:obj.trueSize(1)); zticks(0:10:obj.trueSize(3));
                ax = gca; ax.FontSize = 16; 
                xlim([0 obj.trueSize(1)]); ylim([0 obj.trueSize(2)]); zlim([0 obj.trueSize(3)])
                clim([0,densityMax])
                colormap(gca,jet)
                view([242.1011   14.4475])
                set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
                grid on, axis image
            
                %saveas(gcf,sprintf('%s/images_%d_tind_%d.png',cfg.out_dir,i,(i-cfg.first_time)/cfg.time_jump+1))
            end
        end
        
    end

    methods (Access = protected)

        function footer = getFooter(obj)
            % sets a custom footer in the display of an object
            footer = [newline, ...
                '   Show <a href="matlab: methods(''', class(obj), '''); ">all methods', ...
                    newline, '</a>'];
        end

        function propertyGroup = getPropertyGroups(obj)
            %defines a group of properties to display
            propertyList = ["tag", "dataset_name", "sigma", "beta", "nt" ...
                        , "nIterations_pcg", "trueSize", "do_resize", "sizeFactor" ...
                            , "first_time", "last_time", "time_jump"];
            propertyGroup = matlab.mixin.util.PropertyGroup(propertyList);
        end

        function M = calculateM(obj, RHO0, U, R)
        % calculates intermediary matrix M used in getVelocity()
        % 
        % M.S is a matrix, that when multiplied against a density
        % distribution, returns a new, advected density disitrubution. 
        %
        % M.Tx, M.Ty are spatial derivatives (should be renamed
        % because T refers to density, which is the variable used in the
        % dTrilinears3d function but not in the rest of this code. 
        % perhaps rename to getAdvectionIntermediaryMatrices(). 
        %
        % M.Tr is the derivative with respect to the source term.

            for k = 1:obj.nt
                U1 = reshape(U(1:prod(obj.n), k), obj.n');
                U2 = reshape(U(prod(obj.n) + 1:end, k), obj.n');
                if ~obj.addSource
                    [M.S{k}, M.Tx{k}, M.Ty{k}] = dTrilinears2d(RHO0(:, k), ...
                        obj.Xc + obj.dt * U1, obj.Yc + obj.dt * U2, ...
                        1, 1, obj.boundaryConditions);
                else
                    M.R{k} = 1 + obj.dt * R(:, k);
                    M.Rho{k} = RHO0(:, k);
                    [M.S{k}, M.Tx{k}, M.Ty{k}] = dTrilinears2d(M.R{k}.*RHO0(:, k), ...
                        obj.Xc + obj.dt * U1, obj.Yc + obj.dt * U2, ...
                        1, 1, obj.boundaryConditions);                    
                end
            end

        end

        function mk = getMongeKantorovichFunctional(obj, rho, u)

            % Vectorization of velocity and density (capitalized to notate vecotorization)
            U = reshape(u, obj.dimension * prod(obj.n), obj.nt);
            Density = reshape(rho, prod(obj.n), obj.nt);

            mk = 0;

            for i = 1:obj.nt
                mk = mk + obj.dt * Density(:, i)' * obj.A * (U(:, i) .* U(:, i));
            end

        end

        function getCellCenteredData(obj)
            % 2D version
            CellDimension1 = ones(obj.trueSize(1), 1);
            CellDimension2 = ones(obj.trueSize(2), 1);
            [obj.Xc, obj.Yc] = getCellCenteredGrid(CellDimension1, CellDimension2);
            BC = {'ccn' 'ccn'};
            obj.Grad = getCellCenteredGradMatrix(BC, CellDimension1, CellDimension2);
        end
        
        function columnVector = splitVelocityAndSource(obj,matrix)
            %SPLITVELOCITYANDSOURCE tranform [u;r] matrix into [u(:);r(:)]
            %column vector. Here u and r are two matrices.
            matrixVelocity = matrix(1:obj.dimension*prod(obj.n),:);
            matrixSource = matrix(obj.dimension*prod(obj.n)+1:end,:);
            columnVector = [matrixVelocity(:);matrixSource(:)];
        end
        
        function diffusionResult = calculateAnisotropicDiffusion(obj,rho)
            % CALCULATEANISOTROPICDIFFUSION calculate the diffusion term
            % \nabla\cdot(\sigma(abs(\nabla rho))\nabla rho).
            % Here we only consider cell-centered neumann boundary
            % conditions(ccn) and \sigma(x)=\sigma_0/(1+|x|^2).
            Rho = reshape(rho,obj.domainSize);
            % boundary condition
            bc = {'ccn','ccn'};
            % spatial grid size, h = [dx,dy].
            h = {ones(1,obj.domainSize(2)-1), ones(1,obj.domainSize(1)-1)}; 
            differenceXMinus = obj.doXMinusDifference(Rho,bc{1},h{1});
            differenceXPlus = obj.doXPlusDifference(Rho,bc{1},h{1});
            differenceYMinus = obj.doYMinusDifference(Rho,bc{2},h{2});
            differenceYPlus = obj.doYPlusDifference(Rho,bc{2},h{2});
            diffusionResult = obj.doXMinusDifference(obj.sigmaFunctional(differenceXPlus.^2 + ...
                obj.minmod(differenceYMinus,differenceYPlus).^2).*differenceXPlus,bc{1},h{1}) + ...
                obj.doYMinusDifference(obj.sigmaFunctional(differenceYPlus.^2 + ...
                obj.minmod(differenceXMinus,differenceXPlus).^2).*differenceYPlus,bc{2},h{2});
            diffusionResult = diffusionResult(:);
        end
        
        function diffusionJacobian = calculateDiffusionJacobian(obj,rho)
            %CALCULATEDIFFUSIONJACOBIAN calculate the jacobian matrix of
            %anisotropic diffuion term D w.r.t the density rho.
            % rho(i,j) i index for y axis; j index for x axis.
            Rho = reshape(rho,obj.domainSize);
            % boundary condition
            bc = {'ccn','ccn'};
            % spatial grid size, h = [dx,dy].
            h = {ones(1,obj.domainSize(2)-1), ones(1,obj.domainSize(1)-1)}; 
            differenceXMinus = obj.doXMinusDifference(Rho,bc{1},h{1});
            differenceXPlus = obj.doXPlusDifference(Rho,bc{1},h{1});
            differenceYMinus = obj.doYMinusDifference(Rho,bc{2},h{2});
            differenceYPlus = obj.doYPlusDifference(Rho,bc{2},h{2});
            
            % gradient rho at different place
            gradientXPlus = differenceXPlus.^2 + obj.minmod(differenceYMinus,differenceYPlus).^2;
            gradientXMinus = zeros(size(gradientXPlus)); 
            gradientXMinus(:,2:end) = gradientXPlus(:,1:end-1);
            gradientXMinus(:,1) = obj.minmod(differenceYMinus(:,1),differenceYPlus(:,1)).^2;
            
            gradientYPlus = differenceYPlus.^2 + obj.minmod(differenceXMinus,differenceXPlus).^2;
            gradientYMinus = zeros(size(gradientYPlus)); 
            gradientYMinus(2:end,:) = gradientYPlus(1:end-1,:);
            gradientYMinus(1,:) = obj.minmod(differenceXMinus(1,:),differenceXPlus(1,:)).^2;
            
            % diffusion coefficient at different place
            sigmaXPlus = obj.sigmaFunctional(gradientXPlus);
            sigmaXMinus = obj.sigmaFunctional(gradientXMinus);            
            sigmaYPlus = obj.sigmaFunctional(gradientYPlus);
            sigmaYMinus = obj.sigmaFunctional(gradientYMinus);
            
            % diffusion coefficient derivative at different place
            dsigmaXPlus = obj.sigmaFunctionalDerivative(gradientXPlus);
            dsigmaXMinus = obj.sigmaFunctionalDerivative(gradientXMinus);
            dsigmaYPlus = obj.sigmaFunctionalDerivative(gradientYPlus);
            dsigmaYMinus = obj.sigmaFunctionalDerivative(gradientYMinus);
            
            RhoGhost = obj.addGhostPoints(Rho);
            % \partial D_{i,j}/\partial \rho_{i,j}
            fijRhoij = -sigmaXPlus - sigmaXMinus + differenceXPlus.*dsigmaXPlus.*(-2*differenceXPlus + ...
                obj.minmod2Derivative(RhoGhost(3:end,2:end-1),Rho,RhoGhost(1:end-2,2:end-1),2)) - ...
                2*differenceXMinus.^2.*dsigmaXMinus - sigmaYPlus - sigmaYMinus + ...
                differenceYPlus.*dsigmaYPlus.*(-2*differenceYPlus + obj.minmod2Derivative(RhoGhost(2:end-1,3:end),Rho,RhoGhost(2:end-1,1:end-2),2)) - ...
                2*differenceYMinus.^2.*dsigmaYMinus;
            % \partial D_{i,j}/\partial \rho_{i,j+1}
            fijRhoijPlus = sigmaXPlus + 2*(differenceXPlus.^2).*dsigmaXPlus + ...
                differenceYPlus.*obj.minmod2Derivative(RhoGhost(2:end-1,3:end),Rho,RhoGhost(2:end-1,1:end-2),1).*dsigmaYPlus;
            % \partial D_{i,j}/\partial \rho_{i,j-1}
            fijRhoijMinus = sigmaXMinus - differenceXMinus.*dsigmaXMinus.*(-2*differenceXMinus +  ...
                obj.minmod2Derivative(RhoGhost(3:end,1:end-2),RhoGhost(2:end-1,1:end-2),RhoGhost(1:end-2,1:end-2),2)) + ...
                differenceYPlus.*obj.minmod2Derivative(RhoGhost(2:end-1,3:end),Rho,RhoGhost(2:end-1,1:end-2),3).*dsigmaYPlus;
            % \partial D_{i,j}/\partial \rho_{i+1,j}
            fijRhoiPlusj = sigmaYPlus + 2*(differenceYPlus.^2).*dsigmaYPlus + differenceXPlus.* ...
                obj.minmod2Derivative(RhoGhost(3:end,2:end-1),Rho,RhoGhost(1:end-2,2:end-1),1).*dsigmaXPlus;
            % \partial D_{i,j}/\partial \rho_{i-1,j}
            fijRhoiMinusj = sigmaYMinus - differenceYMinus.*dsigmaYMinus.*(-2*differenceYMinus + ...
                obj.minmod2Derivative(RhoGhost(1:end-2,3:end),RhoGhost(1:end-2,2:end-1),RhoGhost(1:end-2,1:end-2),2)) + ...
                differenceXPlus.*obj.minmod2Derivative(RhoGhost(3:end,2:end-1),Rho,RhoGhost(1:end-2,2:end-1),3).*dsigmaXPlus;
            % \partial D_{i,j}/\partial \rho_{i-1,j+1}
            fijRhoiMinusjPlus = -differenceYMinus.*obj.minmod2Derivative(RhoGhost(1:end-2,3:end),RhoGhost(1:end-2,2:end-1),RhoGhost(1:end-2,1:end-2),1).*dsigmaYMinus;
            % \partial D_{i,j}/\partial \rho_{i+1,j-1}
            fijRhoiPlusjMinus = -differenceXMinus.*obj.minmod2Derivative(RhoGhost(3:end,1:end-2),RhoGhost(2:end-1,1:end-2),RhoGhost(1:end-2,1:end-2),1).*dsigmaXMinus;
            % \partial D_{i,j}/\partial \rho_{i-1,j-1}
            fijRhoiMinusjMinus = -differenceYMinus.*obj.minmod2Derivative(RhoGhost(1:end-2,3:end),RhoGhost(1:end-2,2:end-1),RhoGhost(1:end-2,1:end-2),3).*dsigmaYMinus - ...
                differenceXMinus.*obj.minmod2Derivative(RhoGhost(3:end,1:end-2),RhoGhost(2:end-1,1:end-2),RhoGhost(1:end-2,1:end-2),3).*dsigmaXMinus;
            
            %Put information into the Jacobian sparse matrix jacobian(D)
            nRho = prod(obj.domainSize);
            ii1 = 1:nRho; jj1 = 1:nRho; uu1 = fijRhoij(:);
            A1 = sparse(ii1,jj1,uu1);
            
            [X,Y] = meshgrid(1:obj.domainSize(2),1:obj.domainSize(1)-1);
            ind2 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii2 = ii1(ind2); jj2 = ii2+1; uu2 = fijRhoiPlusj(:); uu2 = uu2(ind2);
            A2 = sparse(ii2,jj2,uu2);
            
            [X,Y] = meshgrid(1:obj.domainSize(2),2:obj.domainSize(1));
            ind3 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii3 = ii1(ind3); jj3 = ii3-1; uu3 = fijRhoiMinusj(:); uu3 = uu3(ind3);
            A3 = sparse(ii3,jj3,uu3);
            
            [X,Y] = meshgrid(1:obj.domainSize(2)-1,1:obj.domainSize(1));
            ind4 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii4 = ii1(ind4); jj4 = ii4+obj.domainSize(1); uu4 = fijRhoijPlus(:); uu4 = uu4(ind4);
            A4 = sparse(ii4,jj4,uu4);
            
            [X,Y] = meshgrid(2:obj.domainSize(2),1:obj.domainSize(1));
            ind5 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii5 = ii1(ind5); jj5 = ii5-obj.domainSize(1); uu5 = fijRhoijMinus(:); uu5 = uu5(ind5);
            A5 = sparse(ii5,jj5,uu5);
            
            [X,Y] = meshgrid(2:obj.domainSize(2),1:obj.domainSize(1)-1);
            ind6 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii6 = ii1(ind6); jj6 = ii6-obj.domainSize(1)+1; uu6 = fijRhoiPlusjMinus(:); uu6 = uu6(ind6);
            A6 = sparse(ii6,jj6,uu6);
            
            [X,Y] = meshgrid(1:obj.domainSize(2)-1,2:obj.domainSize(1));
            ind7 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii7 = ii1(ind7); jj7 = ii7+obj.domainSize(1)-1; uu7 = fijRhoiMinusjPlus(:); uu7 = uu7(ind7);
            A7 = sparse(ii7,jj7,uu7);
            
            [X,Y] = meshgrid(2:obj.domainSize(2),2:obj.domainSize(1));
            ind8 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii8 = ii1(ind8); jj8 = ii8-obj.domainSize(1)-1; uu8 = fijRhoiMinusjMinus(:); uu8 = uu8(ind8);
            A8 = sparse(ii8,jj8,uu8);
            
            diffusionJacobian = A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8;            
        end    
        
        function matrixOut = sigmaFunctionalDerivative(obj,matrixIn,K)
            %SIGMAFUNCTIONALDERIVATIVE returns 1st order derivative of
            %\sigma(x)
            %x could be \rho or |\nabla\rho|
            %case |\nabla\rho|
            %obj.sigmaFunctionType=1 \sigma'(x) =
            %-2\sigma_0*x/(K^2(1+(x/K)^2)^2)
            %obj.sigmaFunctionType=2 \sigma'(x) =
            %-2\sigma_0*x*exp(-(x/K)^2)/K^2
            %To reduce calculation cost, we always assume \sigma(x) is a
            %function of x^2 and matrixIn=x^2. Since 2x will be eliminated
            %after taking \partial |\nabla rho|/\partial \rho, the result
            %here is for example matrixOut=-\sigma_0/(K^2(1+matrixIn/K^2)^2)
            %or matrixOut = -\sigma_0*exp(-matrixIn/K^2)/K^2
            %case \rho
            %obj.sigmaFunctionType=0 \sigma'(x) = \sigma_0/K
            %obj.sigmaFunctionType=2 \sigma'(x) = -\sigma_0/(K*(1+x/K)^2)
            if nargin==2
                K = obj.sigmaFunctionConstant;
            end
            switch obj.diffusionCoefficientType
                case {'anisotropic', 'autoAnisotropic'}
                    switch obj.sigmaFunctionType
                        case 1
                            matrixOut = -obj.sigma./(K^2 * ...
                                (1+matrixIn/K^2).^2);
                        case 2
                            matrixOut = -obj.sigma*exp(-matrixIn/K^2)/K^2;
                    end
                case 'brightness'
                    switch obj.sigmaFunctionType
                        case 0
                            matrixOut = obj.sigma/K;
                        case 1
                            matrixOut = -obj.sigma./(K*(1+matrixIn/K).^2);
                    end
            end
        end
        
        function diffusionResult = calculateBrightnessDiffusion(obj,rho)
            % CALCULATEBRIGHTNESSDIFFUSION calculate the brightness
            % diffusion term, i.e. \nabla\cdot(\rho\nabla\rho)
            % Here we only consider cell-centered neumann boundary
            % conditions(ccn).
            Rho = reshape(rho,obj.domainSize);
            % boundary condition
            bc = {'ccn','ccn'}; 
            h = {ones(1,obj.domainSize(2)-1), ones(1,obj.domainSize(1)-1)}; 
            differenceXPlus = obj.doXPlusDifference(Rho,bc{1},h{1});
            differenceYPlus = obj.doYPlusDifference(Rho,bc{2},h{2}); 
            RhoGhost = obj.addGhostPoints(Rho);
            RhoXPlus = RhoGhost(2:end-1,3:end);
            RhoYPlus = RhoGhost(3:end,2:end-1);
            diffusionResult = obj.doXMinusDifference(obj.sigmaFunctional((Rho+RhoXPlus)/2).*differenceXPlus,bc{1},h{1}) + ...
                obj.doYMinusDifference(obj.sigmaFunctional((Rho+RhoYPlus)/2).*differenceYPlus,bc{2},h{2});
            diffusionResult = diffusionResult(:);
        end
        
        function diffusionJacobian = calculateBDiffusionJacobian(obj,rho)
            % CALCULATEBDIFFUSIONJACOBIAN calculate jacobian of brightness
            % diffusion w.r.t density rho.
            Rho = reshape(rho,obj.domainSize);
            % boundary condition
            bc = {'ccn','ccn'}; 
            h = {ones(1,obj.domainSize(2)-1), ones(1,obj.domainSize(1)-1)}; 
            differenceXPlus = obj.doXPlusDifference(Rho,bc{1},h{1});
            differenceYPlus = obj.doYPlusDifference(Rho,bc{2},h{2});
            differenceXMinus = obj.doXMinusDifference(Rho,bc{1},h{1});
            differenceYMinus = obj.doYMinusDifference(Rho,bc{2},h{2});
            RhoGhost = obj.addGhostPoints(Rho);
            meanRhoXPlus = (RhoGhost(2:end-1,3:end)+Rho)/2;
            meanRhoYPlus = (RhoGhost(3:end,2:end-1)+Rho)/2;
            meanRhoXMinus = (RhoGhost(2:end-1,1:end-2)+Rho)/2;
            meanRhoYMinus = (RhoGhost(1:end-2,2:end-1)+Rho)/2;
            
            
            DijRhoij = obj.sigmaFunctionalDerivative(meanRhoXPlus).*differenceXPlus/2 - ...
                obj.sigmaFunctional(meanRhoXPlus) - obj.sigmaFunctional(meanRhoXMinus) - ...
                obj.sigmaFunctionalDerivative(meanRhoXMinus).*differenceXMinus/2 + ...
                obj.sigmaFunctionalDerivative(meanRhoYPlus).*differenceYPlus/2 - ...
                obj.sigmaFunctional(meanRhoYPlus) - obj.sigmaFunctional(meanRhoYMinus) - ...
                obj.sigmaFunctionalDerivative(meanRhoYMinus).*differenceYMinus/2;
            
            DijRhoiPlusj = obj.sigmaFunctionalDerivative(meanRhoYPlus).*differenceYPlus/2 + ...
                obj.sigmaFunctional(meanRhoYPlus);           
            DijRhoiMinusj = -obj.sigmaFunctionalDerivative(meanRhoYMinus).*differenceYMinus/2 + ...
                obj.sigmaFunctional(meanRhoYMinus);            
            DijRhoijPlus = obj.sigmaFunctionalDerivative(meanRhoXPlus).*differenceXPlus/2 + ...
                obj.sigmaFunctional(meanRhoXPlus);            
            DijRhoijMinus = -obj.sigmaFunctionalDerivative(meanRhoXMinus).*differenceXMinus/2 + ...
                obj.sigmaFunctional(meanRhoXMinus);
                       
            nRho = prod(obj.domainSize);
            ii1 = 1:nRho; jj1 = 1:nRho; uu1 = DijRhoij(:);
            A1 = sparse(ii1,jj1,uu1,nRho,nRho);
            
            
            [X,Y] = meshgrid(1:obj.domainSize(2),1:obj.domainSize(1)-1);
            ind2 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii2 = ii1(ind2); jj2 = ii2+1; uu2 = DijRhoiPlusj(:); uu2 = uu2(ind2);
            A2 = sparse(ii2,jj2,uu2,nRho,nRho);
            
            [X,Y] = meshgrid(1:obj.domainSize(2),2:obj.domainSize(1));
            ind3 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii3 = ii1(ind3); jj3 = ii3-1; uu3 = DijRhoiMinusj(:); uu3 = uu3(ind3);
            A3 = sparse(ii3,jj3,uu3,nRho,nRho);
            
            [X,Y] = meshgrid(1:obj.domainSize(2)-1,1:obj.domainSize(1));
            ind4 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii4 = ii1(ind4); jj4 = ii4+obj.domainSize(1); uu4 = DijRhoijPlus(:); uu4 = uu4(ind4);
            A4 = sparse(ii4,jj4,uu4,nRho,nRho);
            
            [X,Y] = meshgrid(2:obj.domainSize(2),1:obj.domainSize(1));
            ind5 = sub2ind(obj.domainSize,Y(:),X(:));           
            ii5 = ii1(ind5); jj5 = ii5-obj.domainSize(1); uu5 = DijRhoijMinus(:); uu5 = uu5(ind5);
            A5 = sparse(ii5,jj5,uu5,nRho,nRho);
            
            diffusionJacobian = A1 + A2 + A3 + A4 + A5;  
        end
        
    end

    methods (Static)

        function checkDimension(dimension)

            arguments
                dimension (1, 1) double
            end

            if ~(dimension == 2)
                error("Dimension of data should be 2. \n For 3D data use MRI3D subclass")
            end

        end

        function createRecordFile(filename, m, headers, r, c)
            % This function functions like the build in MATLAB function csvwrite but
            % allows a row of headers to be easily inserted
            %
            % known limitations
            % 	The same limitation that apply to the data structure that exist with
            %   csvwrite apply in this function, notably:
            %       m must not be a cell array
            %
            % Inputs
            %
            %   filename    - Output filename
            %   m           - array of data
            %   headers     - a cell array of strings containing the column headers.
            %                 The length must be the same as the number of columns in m.
            %   r           - row offset of the data (optional parameter)
            %   c           - column offset of the data (optional parameter)
            %
            %
            % Outputs
            %   None

            arguments
                filename char
                m
                headers char
                r = 0
                c = 0
            end

            if length(headers) ~= size(m, 2)
                error('number of header entries must match the number of columns in the data')
            end

            %% write the header string to the file

            %turn the headers into a single comma seperated string
            header_string = join(string(headers), ',');

            %write the string to a file
            fid = fopen(filename, 'w');
            fprintf(fid, '%s\r\n', header_string);
            fclose(fid);

            %% write the append the data to the file

            %
            % Call dlmwrite with a comma as the delimiter
            %
            dlmwrite(filename, m, '-append', 'delimiter', ',', 'roffset', r, 'coffset', c);
        end
        
        function c = minmod(a,b)
            % MINMOD find the second large value among a,b,0
            % a,b need to be the same size.
            c = ((a>0) & (b>0)).*min(a,b) + ((a<0) & (b<0)).*max(a,b);
        end
        
        function matrixOut = doXPlusDifference(matrixIn,BoundaryConditionType,h)
            [~,h2] = size(matrixIn);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(:,2:end-1);
                matrixNPlus = matrixIn(:,3:end);
                matrixOut = (matrixNPlus-matrixN)*diag(1./h(2:end));
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(:,1:end-1);
                matrixNPlus = matrixIn(:,2:end);
                matrixOut = (matrixNPlus-matrixN)*diag(1./h);
                matrixOut(:,h2) = 0;
            end
        end
        
        function matrixOut = doXMinusDifference(matrixIn,BoundaryConditionType,h)
            [~,h2] = size(matrixIn);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(:,2:end-1);
                matrixNMinus = matrixIn(:,1:end-2);
                matrixOut = (matrixN - matrixNMinus)*diag(1./h(1:end-1));
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(:,2:end);
                matrixNMinus = matrixIn(:,1:end-1);
                matrixOut(:,2:h2) = (matrixN - matrixNMinus)*diag(1./h);
                matrixOut(:,1) = 0;
            end
        end
        
        function matrixOut = doYPlusDifference(matrixIn,BoundaryConditionType,h)
            [h1,~] = size(matrixIn);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(2:end-1,:);
                matrixNPlus = matrixIn(3:end,:);
                matrixOut = diag(1./h(2:end))*(matrixNPlus-matrixN);
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(1:end-1,:);
                matrixNPlus = matrixIn(2:end,:);
                matrixOut = diag(1./h)*(matrixNPlus-matrixN);
                matrixOut(h1,:) = 0;
            end
        end
        
        function matrixOut = doYMinusDifference(matrixIn,BoundaryConditionType,h)
            [h1,~] = size(matrixIn);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(2:end-1,:);
                matrixNMinus = matrixIn(1:end-2,:);
                matrixOut = diag(1./h(1:end-1))*(matrixN - matrixNMinus);
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(2:end,:);
                matrixNMinus = matrixIn(1:end-1,:);
                matrixOut(2:h1,:) = (1./h)*(matrixN - matrixNMinus);
                matrixOut(1,:) = 0;
            end
        end
        
        function d = minmod2Derivative(a,b,c,variableIndex)
            %MINMOD2DERIVATIVE calculate the first order of derivative of
            %function minmod(a-b,b-c)^2. This derivative is taken w.r.t
            %variableIndex variable.
            
            switch variableIndex
                case 1
                    d = 2*((a-b>0 & b-c>a-b) | (a-b<0 & b-c<a-b)).*(a-b);
                case 2
                    d = 2*(((a-b>0 & b-c>a-b) | (a-b<0 & b-c<a-b)).*(b-a) + ...
                        ((b-c>0 & b-c<a-b) | (b-c<0 & b-c>a-b)).*(b-c));
                case 3
                    d = 2*((b-c>0 & b-c<a-b) | (b-c<0 & b-c>a-b)).*(c-b);
            end
        end
        
        function rhoGhost = addGhostPoints(rho)
            rhoGhost = zeros(size(rho)+2);
            rhoGhost(2:end-1,2:end-1) = rho;
            rhoGhost(1,:) = rhoGhost(2,:);
            rhoGhost(end,:) = rhoGhost(end-1,:);
            rhoGhost(:,1) = rhoGhost(:,2);
            rhoGhost(:,end) = rhoGhost(:,end-1);
        end
        
        function C = sig2conc(img)
            % functon to calculate contrast concentration from image intensity
            % equation: 
            % R1(t)=-1/TR*ln(1-((S(t)-S(0))/S0*sin(alpha))+(1-m)/(1-m*cos(alpha)))
            %over 1-cos(alpha)*((S(t)-S(0))/S0*sin(alpha))+(1-m)/(1-m*sin(alpha)))
            % where m=exp(-R10*TR)
            % then C(t)=(R1(t)-R1(0))/r1
            
            % Yi Guo, 06/12/2014
            % No mask comparing to Marc's version, otherwise allmost the same
            % some simplification, pay attention to R1, and alpha unit!!!
            % See source: https://github.com/shwetha0312/DCE_MRI/blob/c99f907db74206f4f8a2aae6c38df4d73e9288cb/conc2sig.m
            
            trueSize = size(img);
            imgB=img(:,:,:,1);
            
            alpha = pi*15/180; %angle
            TR = 0.006; 
            M0 = 5*ones(trueSize(1:2));
            R10 = ones(trueSize(1:2));
            
            Rcs=4.39; % contrast relaxivity
            nt=size(img,4); % temporal dimension size
            
            m=exp(-repmat(R10,[1 1 1 nt])*TR);
            par2=(1-m)./(1-m*cos(alpha));
            
            par1=(img-repmat(imgB,[1 1 1 nt]))./repmat(M0+eps,[1 1 1 nt])/sin(alpha);
            
            B=(par1+par2);
            
            Rt=-1/TR*real(log((1-B)./(1-cos(alpha)*B+eps)));
            
            R1B=Rt(:,:,:,1); % baseline R1 is equal to R10 if img(1)==imgB
            %R1B=R10;
            C=Rt-repmat(R1B,[1 1 1 nt]);
            
            C=C/Rcs;
            
            C(C<0) = 0;  % get rid of outliers
            C(C>180)=0;
        end
    end

end
