function [obj, flag] = runRegularizedOMT(tag, dimension)
%%RUNROMT Runs regularized optimal mass transport algorithm (rOMT)
% MRI image series containing dynamic tracking of injected contrast
% agent are processed for extracting estimates of contrast flow, measured
% in a velocity vector field. 
% 
% Results are saved in an output folder, to be further postprocessed and     
% visualized with figures in the Postprocessing.m script. 
%
%   tag         - name of the experiment case
%   dimension   - 2 (2D data) or 3 (3D)

    arguments
        tag string
        dimension (1,1) double = 3
    end
    
    addpath('./Sensitivities', './analyzeFlows', genpath('./utilities'))

    % construct class instance based on dimension
    if dimension == 3
        obj = Configuration3D(tag);
    elseif dimension == 2
        obj = Configuration(tag);
    else
        warning("dimension must be 2 or 3")
    end

    % load volumetric data
    obj.loadVolume; %includes preprocessing (3D affine diffusion smoothing)    

    % initializes and saves density at initial time point
    densityCached = obj.Data(1).volume(:); % first image (flattened vector of the image matrix)
    obj.saveInitialDensity;

    fprintf('\n =============== rOMT Starts ===============\n')
    fprintf('______________________________________________\n\n')
    disp(obj); %displays important properties of the object in Configuration class

    %profile on
    clear T;
    T = 0;

    % initialize cell centered grid and gradients
    obj.initializeGrid;

    %calculate relative total intensity change
    if obj.addSource
        switch obj.alphaFunctionalType
            case 'auto'
                Rho = zeros(1,obj.nFiles);
                mask_brain = Mask(obj.pathMask, obj.isMaskFilled, obj.xRange,...
                            obj.yRange, obj.zRange);
                if obj.do_resize
                    mask_brain = mask_brain.resize(obj.sizeFactor);
                end
                for i = 1:obj.nFiles+1
                    A = obj.Data(i).volume;
                    A(mask_brain.contents == 0) = 0;
                    Rho(i) = sum(A,'all');
                end
                Rho_diff = abs(Rho(2:end)-Rho(1:end-1))./Rho(1:end-1); 
            case 'manual'
                alphaNum = load(obj.pathAlpha);
        end
    end
    
    % Iterates over time series of images
    nTimePoints = length(obj.timeInitial:obj.timeJump:obj.timeFinal);    
    for iTime = 1:nTimePoints
        tic
        % Reinitialization decouples the image pairs from their temporal
        % relationship.
        if obj.reinitializeDensityFlag %
            densitySource = obj.Data(iTime).volume(:);
        else
            densitySource = densityCached(:);
        end
        if strcmp(obj.diffusionCoefficientType, 'autoAnisotropic')
            obj.calculateK(densitySource);
        end
        % true final density from target image in the loop
        obj.densityTarget = obj.Data(iTime + 1).volume(:);
        if ~obj.addSource
            % case of not adding source term
            % initial guess for velocity (all zero)
            if iTime == 1 || obj.reinitializeVelocityFlag
                velocity = zeros(length(obj.trueSize) * prod(obj.trueSize), obj.nt);
            end

            %% Descent for velocity
            fprintf('\n =============== Descent on u ===============\n')
            fprintf('time index = %d / %d \n', iTime, nTimePoints)
            fprintf('______________________________________________\n\n')
            fprintf('iLineSearch \tphi    \t      descent output\n')
            fprintf('________    ___________     __________________\n')

            velocity = obj.getVelocity(densitySource, velocity);

            %Calculate energy functional (Phi)
            [phi, mk, phiN, Rho, Ru] = obj.getPhi(densitySource, velocity);

            % get final interpolated image from nt # of rhos
            densityCached = Rho(:, end); % Rho = [rho_0, rho2,...,rho_nt]
            btoc = toc;
            T = T + btoc;

            %% Save results for each image pair

            % save record of parameters 
            timeOfDensityOrigin = obj.timeInitial + (iTime - 1) * obj.timeJump;
            timeOfDensityTarget = timeOfDensityOrigin + obj.timeJump;
            obj.saveRecord([iTime, timeOfDensityOrigin, timeOfDensityTarget, ...
                            phi, mk, Ru, phiN, max(velocity(:)), btoc]);

            % save optimal velocity and density matrices
            obj.saveVariable(timeOfDensityOrigin, timeOfDensityTarget, iTime, velocity)
            obj.saveVariable(timeOfDensityOrigin, timeOfDensityTarget, iTime, densityCached)

            fprintf('time index = %d, max(u) = %5.4f\n', iTime, max(velocity));
        else
            % case of adding source term
            % set alpha value
            switch obj.alphaFunctionalType
                case 'auto'
                    if Rho_diff(iTime)>=0.5
                        obj.alpha = 25000;
                    else
                        obj.alpha = 75000/(4*Rho_diff(iTime)+1);
                    end
                case 'manual'
                    obj.alpha = (alphaNum.alpha(1+(iTime - 1) * obj.timeJump) + alphaNum.alpha(1 + iTime * obj.timeJump))/2;
            end
            % initial guess for velocity (all zero)           
            if iTime == 1 || obj.reinitializeVelocityFlag
                velocity = zeros(length(obj.trueSize) * prod(obj.trueSize), obj.nt);
            end
            % initial guess for source (all zero)
            if iTime == 1 || obj.reinitializeSourceFlag
                source = zeros(prod(obj.trueSize), obj.nt);
            end
            %% Descent for velocity and source
            fprintf('\n =============== Descent on u and r ===============\n')
            fprintf('time index = %d / %d \n', iTime, nTimePoints)
            fprintf('______________________________________________\n\n')
            fprintf('iLineSearch \tphi    \t      descent output\n')
            fprintf('________    ___________     __________________\n')

            [velocity, source] = obj.getVelocityAndSource(densitySource, velocity, source);

            %Calculate energy functional (Phi)
            [phi, mk, phiN, Rho, Ru, sourceEnergy] = obj.getPhi(densitySource, velocity, source);

            % get final interpolated image from nt # of rhos
            densityCached = Rho(:, end); % Rho = [rho_0, rho2,...,rho_nt]
            btoc = toc;
            T = T + btoc;

            %% Save results for each image pair

            % save record of parameters 
            timeOfDensityOrigin = obj.timeInitial + (iTime - 1) * obj.timeJump;
            timeOfDensityTarget = timeOfDensityOrigin + obj.timeJump;
            obj.saveRecord([iTime, timeOfDensityOrigin, timeOfDensityTarget, ...
                            phi, mk, Ru, phiN, sourceEnergy,max(velocity(:)), ...
                            max(source(:)), btoc]);

            % save optimal velocity and density matrices
            obj.saveVariable(timeOfDensityOrigin, timeOfDensityTarget, iTime, velocity)
            obj.saveVariable(timeOfDensityOrigin, timeOfDensityTarget, iTime, densityCached)
            obj.saveVariable(timeOfDensityOrigin, timeOfDensityTarget, iTime, source)

            fprintf('time index = %d, max(u) = %5.4f, max(r) = %5.4f\n', iTime, max(velocity), max(source));
        end
    end
    if strcmp(obj.diffusionCoefficientType,'autoAnisotropic')
        KList = obj.KList;
        save(sprintf('%s_KList.mat',obj.pathOutput),'KList');
    end
    fprintf('\n =============== rOMT Ends ===============\n')
    fprintf('\n Elapsed Time: %s\n', datestr(seconds(T), 'HH:MM:SS'))

    flag = 1;
end
