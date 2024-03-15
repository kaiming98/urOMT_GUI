classdef LagrangianSpecification < handle
    %LAGRANGIANSPECIFICATION Lagrangian specification of flow field
    %   Looking motion of an individual partical over time
    %   
    %   Properties defined in Lagrangian coordinates, i.e., in a coordinate
    %   system that follows trajectories of particles from a starting
    %   position (referred here as StartingPoints). 
    %   
    %   The below scalar properties are defined by a starting point (x0 and
    %   time, t): X(x0,t) where X is a scalar property. 
    %
    %   The dimensions of each property X is therefore defined by the # of
    %   startingpoints, Np, and the number of timesteps, Nt, namely [Np,Nt]
    %   matrices. 
    %
    %   The position matrix is defined spatially, so [Np, (x,y,z),Nt] or
    %   [Np, 3, Nt]

    properties
        position % dimension: L
        rho % density, dimension: L^-3M
        speed % dimension: T^-1L
        diffusiveSpeed
        augmentedSpeed
        timeStep % dimension: T
        flux % dimension: T^-1L^-2M
        drho_dt %drho/dt, dimension T^-1L^-3M
        drho_dx %|\nabla rho|
        source %r: relative source
        sFlux %\rho*r: source flux
        diffusionCoefficient
        startPoints
        experimentMask
        analysisConstants
        cutoffThresholds
    end

    properties (Access = protected)
        step
        addSource
    end

    methods

        function obj = LagrangianSpecification(tag,analysisConstants)
            %LAGRANGIANSPECIFICATION Construct an instance of this class
            
            configuration = Configuration3D(tag);
            obj.analysisConstants = analysisConstants;
            obj.addSource = configuration.addSource;
    
            % create starting points instance used in velocity visualization
            obj.startPoints = OrderedStartingPoints(configuration.pathStartpointsMask,...
                                                    configuration.startpointsMaskName);
            
            obj.startPoints.updateJump(analysisConstants.jump);
            % set cutoff threshold
            switch analysisConstants.cutoffFlag
                case 'min'
                    obj.cutoffThresholds.density = 0.0001;
                    obj.cutoffThresholds.front = 0;
                    obj.cutoffThresholds.flow = 0;
                    obj.cutoffThresholds.speed = 0.0001;
                case 'max'
                    obj.cutoffThresholds.density = .0001;
                    obj.cutoffThresholds.front = 0;
                    obj.cutoffThresholds.flow = 0;
                    obj.cutoffThresholds.speed = 0;
                case 'mean'
                    obj.cutoffThresholds.density = .00001;
                    obj.cutoffThresholds.front = 0;
                    obj.cutoffThresholds.flow = 0.0001;
                    obj.cutoffThresholds.speed = 0.0001;
            end
            % create mask
            obj.getmasks(configuration);
            
            % starting points thresholding:Only using the starts points from max(dpsnrv)>startpointsThresh
            obj.startPoints.doThresholding(tag);
            
            % select starting points based on some laws (ordered or uniform).
            obj.startPoints.select(tag);  
            
            % set useful variables
            h1 = 1; h2 = 1; h3 = 1;
            ti = configuration.timeInitial;
            tj = configuration.timeJump;
            tf = configuration.timeFinal;
            nt = configuration.nt;

            %convert from matlab grid to cell-centered grid:
            s1 = (obj.startPoints.y - 0.5) .* h1; %i/y-axis
            s2 = (obj.startPoints.x - 0.5) .* h2; %j/x-axis
            s3 = (obj.startPoints.z - 0.5) .* h3; %k/z-axis
            sp_123 = [s1, s2, s3];
            %current point i.e. list of current location in each
            %streamline that hasn't been terminated
            pointsCurrent = sp_123;
            %keep track of the # of streamlines that have not yet been terminated
            npoints = length(pointsCurrent);

            obj.initialize(npoints, ti, tj, tf, nt, analysisConstants.nEulerianSteps);

            obj.step = 1;
            xt = pointsCurrent; %current position in rcl orientation
            xt = max([h1 * 0.5, h2 * 0.5, h3 * 0.5], ...
                min(xt, [h1 * (configuration.trueSize(1) - .5001), ...
                    h2 * (configuration.trueSize(2) - .5001), h3 * (configuration.trueSize(3) - .5001)]));

            obj.calculate(configuration, xt, analysisConstants)
        end

        function calculate(obj, configuration, xt, analysisConstants)
            %CALCULATE calculate properties of flow field in lagrangian
            %coordinates

            % set useful varibles
            ti = configuration.timeInitial;
            tj = configuration.timeJump;
            tf = configuration.timeFinal;
            nt = configuration.nt;
            [x, y, z] = meshgrid(1:configuration.trueSize(2), 1:configuration.trueSize(1), ...
                1:configuration.trueSize(3));
            if isempty(configuration.Grad)
                configuration.initializeGrid;
            end
            Mdis = -configuration.sigma * configuration.Grad' * configuration.Grad;

            for t1 = ti:tj:tf

                [RHO_t, U, R] = obj.loadRomtResult(t1, configuration, obj.analysisConstants.isInterpolated);
                
                for t2 = 1:nt
                    TInd = ((t1 - ti) / tj) * nt + t2;

                    T = t1 + (t2 - 1) * (tj / nt);
                    fprintf('t = %d (t1 = %d, t2 = %d -> T = %.3f)\n', TInd, t1, t2, T);

                    if obj.analysisConstants.isInterpolated
                        d = reshape(RHO_t(:, t2), configuration.trueSize);
                    else
                        d = getData(tag, round(T), 'none');

                        if configuration.smooth > 0
                            d = affine_diffusion_3d(d, configuration.smooth, 0.1, 1, 1);
                        end

                    end

                    if analysisConstants.intensityRangeFlag
                        d = d - min(d(:));
                    end
                    %make sure density is non-negative
                    d(d<0)=0;
                    
                    u = reshape(U(:, t2), [], 3);
                    v1 = reshape(u(:, 1), configuration.trueSize);
                    v2 = reshape(u(:, 2), configuration.trueSize);
                    v3 = reshape(u(:, 3), configuration.trueSize);

                    %norm of gradient rho
                    [gradientRho2, gradientRho1, gradientRho3] = gradient(d);
                    gradientRhoNorm = gradientRho1.^2+gradientRho2.^2+gradientRho3.^2;
                    %diffusive speed
                    switch configuration.diffusionCoefficientType
                        case 'constant'
                            %add eps to d so can take log(d) and not log(0)
                            [w2, w1, w3] = gradient(log(d + 2 * eps));
                            diffusionCoeff = configuration.sigma .*ones(size(d));
                            w1 = configuration.sigma .*w1;
                            w2 = configuration.sigma .*w2;
                            w3 = configuration.sigma .*w3;
                            du = [w1(:), w2(:), w3(:)];
                        case 'anisotropic'
                            [w2, w1, w3] = gradient(log(d + 2 * eps));
                            sigmaGradientRho = configuration.sigmaFunctional(gradientRhoNorm);
                            diffusionCoeff = sigmaGradientRho;
                            w1 = sigmaGradientRho.*w1;
                            w2 = sigmaGradientRho.*w2;
                            w3 = sigmaGradientRho.*w3;
                            du = [w1(:), w2(:), w3(:)];
                        case 'autoAnisotropic'
                            load(sprintf('%s_KList.mat',configuration.pathOutput));
                            K = KList((t1 - ti) / tj + 1);
                            [w2, w1, w3] = gradient(log(d + 2 * eps));
                            sigmaGradientRho = configuration.sigmaFunctional(gradientRhoNorm,K);
                            diffusionCoeff = sigmaGradientRho;
                            w1 = sigmaGradientRho.*w1;
                            w2 = sigmaGradientRho.*w2;
                            w3 = sigmaGradientRho.*w3;
                            du = [w1(:), w2(:), w3(:)];
                        case 'brightness' %du = \rho \nabla log(\rho)
                            [w2, w1, w3] = gradient(log(d + 2 * eps));
                            diffusionCoeff = configuration.sigmaFunctional(d);
                            w1 = configuration.sigmaFunctional(d).*w1;
                            w2 = configuration.sigmaFunctional(d).*w2;
                            w3 = configuration.sigmaFunctional(d).*w3;
                            du = [w1(:), w2(:), w3(:)];
                    end
                    
                    u1 = v1 - w1;
                    u2 = v2 - w2;
                    u3 = v3 - w3;
                   
                    r = reshape(R(:, t2), configuration.trueSize);
                    
                    switch analysisConstants.flowType
                        case 'vel'
                            a1 = u1;
                            a2 = u2;
                            a3 = u3;
                        case 'flw'
                            a1 = u1 .* d;
                            a2 = u2 .* d;
                            a3 = u3 .* d;
                    end

                    [Gx, Gy, Gz] = gradient(d);
                    speedTmp = reshape(sqrt(sum(u.^2, 2)), configuration.trueSize);
                    diffusiveSpeedTmp = reshape(sqrt(sum(du.^2, 2)), configuration.trueSize);
                    speedAug = sqrt(u1.^2 + u2.^2 + u3.^2);
                    img_flow = speedTmp .* d;
                    src_flux = r .* d;
                    drdt_dif = Mdis * d(:);
                    drdt_ad1 = Gx .* reshape(u(:, 2), configuration.trueSize) + ...
                        Gy .* reshape(u(:, 1), configuration.trueSize) + ...
                        Gz .* reshape(u(:, 3), configuration.trueSize);
                    DIVu = divergence(x, y, z, reshape(u(:, 2), configuration.trueSize), ...
                        reshape(u(:, 1), configuration.trueSize), reshape(u(:, 3), configuration.trueSize));
                    drdt_ad2 = d(:) .* DIVu(:);
                    drdt = reshape(drdt_dif - (drdt_ad1(:) + drdt_ad2), configuration.trueSize);

                    %update first step of density and speed
                    if obj.step == 1
                        obj.doFirstStep(configuration.trueSize, obj.startPoints.x,...
                            obj.startPoints.y, obj.startPoints.z, xt, d, sqrt(gradientRhoNorm),...
                            speedTmp, diffusiveSpeedTmp, speedAug, img_flow, drdt, r, src_flux);
                    end

                    switch analysisConstants.cutoffFlag
                        case 'min'
                            conf.density = 1;
                            conf.speed = 1;
                        case 'max'
                            conf.density = mean(d(d > 0)) + std(d(d > 0));
                            conf.speed = mean(speedTmp(speedTmp > 0)) + std(speedTmp(speedTmp > 0));
                        case 'mean'
                            conf.density = mean(d(d > 0));
                            conf.speed = mean(speedTmp(speedTmp > 0));
                    end

                    %vector field to be integrated in order to compute streamlines
                    V = [a1(:), a2(:), a3(:)];
                    V(obj.experimentMask.contents == 0, :) = 0; %don't want to move outside of the masked ROI
                    u(obj.experimentMask.contents == 0, :) = 0;

                    xt = obj.update(xt, V, u, du, d, sqrt(gradientRhoNorm), img_flow, drdt, r, src_flux, diffusionCoeff, conf, configuration.trueSize, analysisConstants);
                end

            end

        end
                        
        function getmasks(obj,configuration)          
            % GETMASKS get experiment mask with dilate 0 and the mask used 
            % for starting points.
            obj.experimentMask = Mask(configuration.pathMask,configuration.isMaskFilled,...
                configuration.xRange, configuration.yRange, configuration.zRange);

            if configuration.do_resize
                obj.experimentMask = obj.experimentMask.resize(configuration.sizeFactor);
            end
            
            if configuration.dilate > 0
                obj.experimentMask = obj.experimentMask.dilate(configuration.dilate);
            end
            
            if ~isempty(obj.startPoints)
                obj.startPoints.loadMask(obj.experimentMask,configuration);
            else
                error("StartingPoints instance is not exist");
            end
        end

    end

    methods (Access = protected)

        function initialize(obj, npoints, ti, tj, tf, nt, nstep)
            % initialize pathlines
            % npoints: number of starting points
            % ti,tj,tf,nt,nstep:first time, time jump, final time
            % discrete time number, number of eulerian step
            obj.position = NaN(npoints, 3, length(ti:tj:tf) * nt * nstep + 1);
            obj.rho = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            obj.drho_dx = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            obj.speed = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            obj.diffusiveSpeed = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            obj.augmentedSpeed = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            obj.timeStep = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            obj.flux = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            obj.drho_dt = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            obj.diffusionCoefficient = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            if obj.addSource
                obj.source = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
                obj.sFlux = NaN(npoints, length(ti:tj:tf) * nt * nstep + 1);
            end
        end

        function doFirstStep(obj, n, sx, sy, sz, xt, d, drhodx, speed, diffusiveSpeed, augmentedSpeed, imageFlow, drhodt, source, sourceFlux)
            %update first time step of properties
            obj.position(:, :, 1) = xt;
            obj.rho(:, 1) = d(sub2ind(n, sy, sx, sz));
            obj.drho_dx(:, 1) = drhodx(sub2ind(n, sy, sx, sz));
            obj.speed(:, 1) = speed(sub2ind(n, sy, sx, sz));
            obj.diffusiveSpeed(:, 1) = diffusiveSpeed(sub2ind(n, sy, sx, sz));
            obj.augmentedSpeed(:, 1) = augmentedSpeed(sub2ind(n, sy, sx, sz));
            obj.timeStep(:, 1) = obj.step;
            obj.flux(:, 1) = imageFlow(sub2ind(n, sy, sx, sz));
            obj.drho_dt(:, 1) = drhodt(sub2ind(n, sy, sx, sz));
            if obj.addSource
                obj.source(:, 1) = source(sub2ind(n,sy,sx,sz));
                obj.sFlux(:, 1) = sourceFlux(sub2ind(n,sy,sx,sz));
            end
        end

        function xt = update(obj, xt, augmentedVelocity, velocity, diffusiveVelocity, density,...
                gradientdensity, imageFlow, drhodt, source, sourceFlux, diffusionCoefficient, conf, n, flowField)
            % UPDATE update one time step for all properties
            %xt: location
            %conf: cutoff threshold parameter.
            h1 = 1; h2 = 1; h3 = 1;

            for Estep = 1:flowField.nEulerianSteps
                obj.step = obj.step + 1;
                augmentedVelocityInterpolated = interp_vel(augmentedVelocity, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);
                velocityInterpolated = interp_vel(velocity, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);
                diffusiveVelocityInterpolated = interp_vel(diffusiveVelocity, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);

                densityInterpolated = interpF(density, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);
                gradientDensityInterpolated = interpF(gradientdensity, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);
                augmentedspeedInterpolated = sqrt(sum(augmentedVelocityInterpolated.^2, 2));
                speedInterpolated = sqrt(sum(velocityInterpolated.^2, 2));
                diffusiveSpeedInterpolated = sqrt(sum(diffusiveVelocityInterpolated.^2, 2));
                fluxInterpolated = interpF(imageFlow, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);
                drhodtInterpolated = interpF(drhodt, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);
                diffusionCoefficientInterpolated = interpF(diffusionCoefficient, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);

                if obj.addSource
                    sourceInterpolated = interpF(source, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);
                    sFluxInterpolated = interpF(sourceFlux, n, xt(:, 1), xt(:, 2), xt(:, 3), [h1, h2, h3]);
                end

                densityConfInterpolated = densityInterpolated ./ conf.density;
                speedConfInterpolated = speedInterpolated ./ conf.speed;

                indexThresholded = find(densityConfInterpolated > obj.cutoffThresholds.density & speedConfInterpolated > obj.cutoffThresholds.speed);

                if isempty(indexThresholded)
                    break
                else

                    switch flowField.distanceFlag
                        case 'T'
                            xt(indexThresholded, :) = xt(indexThresholded, :) + ...
                            flowField.nTimeInterval .* augmentedVelocityInterpolated(indexThresholded, :);
                        case 'X'
                            xt(indexThresholded, :) = xt(indexThresholded, :) + ...
                            flowField.nTimeInterval .* (augmentedVelocityInterpolated(indexThresholded, :) ./ speedInterpolated);
                    end

                    %make sure it stays in bounds:
                    xt = max([h1 * 0.5, h2 * 0.5, h3 * 0.5], ...
                    min(xt, [h1 * (n(1) - .5001), h2 * (n(2) - .5001), h3 * (n(3) - .5001)]));
                    obj.position(indexThresholded, :, obj.step) = xt(indexThresholded, :);
                    obj.rho(indexThresholded, obj.step) = densityInterpolated(indexThresholded);
                    obj.drho_dx(indexThresholded, obj.step) = gradientDensityInterpolated(indexThresholded);
                    obj.speed(indexThresholded, obj.step) = speedInterpolated(indexThresholded);
                    obj.diffusiveSpeed(indexThresholded, obj.step) = diffusiveSpeedInterpolated(indexThresholded);
                    obj.augmentedSpeed(indexThresholded, obj.step) = augmentedspeedInterpolated(indexThresholded);
                    obj.timeStep(indexThresholded, obj.step) = obj.step;
                    obj.flux(indexThresholded, obj.step) = fluxInterpolated(indexThresholded);
                    obj.drho_dt(indexThresholded, obj.step) = drhodtInterpolated(indexThresholded);
                    obj.diffusionCoefficient(indexThresholded, obj.step) = diffusionCoefficientInterpolated(indexThresholded);
                    if obj.addSource
                        obj.source(indexThresholded, obj.step) = sourceInterpolated(indexThresholded);
                        obj.sFlux(indexThresholded, obj.step) = sFluxInterpolated(indexThresholded);
                    end
                end

            end

        end
        
    end

    methods (Static)

        function [RHO_t, U, R] = loadRomtResult(t1, configuration, isInterpolated)
            % LOADROMTRESULT load romt density and velocity results
            ti = configuration.timeInitial;
            tj = configuration.timeJump;
            nt = configuration.nt;

            U = importdata(sprintf('%s/%s_v_%d_%d_t_%d.mat', ...
                configuration.pathOutput, configuration.tag, t1, ...
                t1 + tj, (t1 - ti) / tj + 1));
            U = reshape(U, [], nt);
            
            if configuration.addSource
                R = importdata(sprintf('%s/%s_s_%d_%d_t_%d.mat', ...
                    configuration.pathOutput, configuration.tag, t1, ...
                    t1 + tj, (t1 - ti) / tj + 1));
                R = reshape(R, [], nt);
            else
                R = zeros(prod(configuration.trueSize),nt);
            end

            if isInterpolated

                if t1 == ti
                    RHO = readmatrix(sprintf('%s/rho_%s_%d_t_0.txt', ...
                        configuration.pathOutput, configuration.tag, ti));
                else
                    RHO = importdata(sprintf('%s/%s_d_%d_%d_t_%d.mat', ...
                        configuration.pathOutput, configuration.tag, ...
                        t1 - tj, t1, (t1 - ti) ./ tj));
                end

                if configuration.nt > 1
                    if ~configuration.addSource
                        RHO_t = [RHO, configuration.advectionDiffusion(RHO, U(:))];
                    else
                        RHO_t = [RHO, configuration.advectionDiffusionSource(RHO, U(:), R(:))];
                    end
                else
                    RHO_t = RHO;
                end

            end

        end

    end

end
