classdef EulerianSpecification < FlowFieldCharacteristic
    %EULERIANSPECIFICATION Eulerian specification of flow field
    %   Looking the snapshot of flow properties
    %   
    %   Properties defined in Eulerian coordinates i.e. in a coordinates
    %   of location and time.
    %
    %   Example(3D): density defined over domain and time, rho =
    %   rho(x,y,z,t). Each scalar property will be stored in an nX*nY*nZ*nTime
    %   matrix. Each vector property will be stored in an
    %   nPoints*dimension*nTime
    
    properties
        velocity
        augmentedVelocity
        diffusiveVelocity
        speed
        augmentedSpeed
        diffusiveSpeed
        Pe
        source
        sFlux % source flux
        diffusionCoefficient
        points
        nData
    end
    
    properties (Access = private)
        anatomyImage
        experimentMaskNifti
        addSource
        indexROI
        isInterpolated = 1
    end
    
    methods
        function obj = EulerianSpecification(tag)
            %EULERIANSPECIFICATION construct an instance of this class
            %   calculate eulerian results of each property.
            configuration = Configuration3D(tag);
            obj = obj@FlowFieldCharacteristic(configuration);
            obj.jsonFileName = tag;
            % set useful varibles
            ti = configuration.timeInitial;
            tj = configuration.timeJump;
            tf = configuration.timeFinal;
            nt = configuration.nt;
            obj.addSource = configuration.addSource;
            if isempty(configuration.Grad)
                configuration.initializeGrid;
            end
            obj.nData = floor((tf - ti) / tj + 1);
            count = 1;
            obj.velocity = zeros(configuration.dimension * prod(configuration.n),obj.nData);
            obj.augmentedVelocity = zeros(configuration.dimension * prod(configuration.n),obj.nData);
            obj.diffusiveVelocity = zeros(configuration.dimension * prod(configuration.n),obj.nData);
            obj.speed = zeros(prod(configuration.n),obj.nData);
            obj.augmentedSpeed = zeros(prod(configuration.n),obj.nData);
            obj.diffusiveSpeed = zeros(prod(configuration.n),obj.nData);
            obj.Pe = zeros(prod(configuration.n),obj.nData);
            obj.source = zeros(prod(configuration.n),obj.nData);
            obj.sFlux = zeros(prod(configuration.n),obj.nData);           
            obj.diffusionCoefficient = zeros(prod(configuration.n),obj.nData);
            obj.points = 1:prod(configuration.n);         
            for t1 = ti:tj:tf
                
                [RHO_t, U, R] = obj.loadRomtResult(t1, configuration, obj.isInterpolated);
                obj.velocity(:, count) = mean(U,2);
                obj.source(:, count) = mean(R,2);
                for t2 = 1:nt
                    TInd = ((t1 - ti) / tj) * nt + t2;

                    T = t1 + (t2 - 1) * (tj / nt);
                    fprintf('t = %d (t1 = %d, t2 = %d -> T = %.3f)\n', TInd, t1, t2, T);

                    if obj.isInterpolated
                        d = reshape(RHO_t(:, t2), configuration.trueSize);
                    else
                        d = getData(tag, round(T), 'none');

                        if configuration.smooth > 0
                            d = affine_diffusion_3d(d, configuration.smooth, 0.1, 1, 1);
                        end

                    end
                    
                    %make sure density is non-negative
                    d(d<0)=0;
                    
                    u = reshape(U(:, t2), [], 3);
                    obj.speed(:, count) = obj.speed(:, count) + sqrt(sum(u.^2, 2));
                    obj.sFlux(:, count) = obj.sFlux(:, count) + R(:, t2).*d(:);
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
                            obj.diffusionCoefficient(:, count) = obj.diffusionCoefficient(:, count) + configuration.sigma*ones(numel(d),1);
                            w1 = configuration.sigma .*w1;
                            w2 = configuration.sigma .*w2;
                            w3 = configuration.sigma .*w3;
                            du = [w1(:), w2(:), w3(:)];                           
                        case 'anisotropic'
                            [w2, w1, w3] = gradient(log(d + 2 * eps));
                            sigmaGradientRho = configuration.sigmaFunctional(gradientRhoNorm);
                            obj.diffusionCoefficient(:, count) = obj.diffusionCoefficient(:, count) + sigmaGradientRho(:);
                            w1 = sigmaGradientRho.*w1;
                            w2 = sigmaGradientRho.*w2;
                            w3 = sigmaGradientRho.*w3;
                            du = [w1(:), w2(:), w3(:)];
                        case 'autoAnisotropic'
                            load(sprintf('%s_KList.mat',configuration.pathOutput));
                            K = KList((t1 - ti) / tj + 1);
                            [w2, w1, w3] = gradient(log(d + 2 * eps));
                            sigmaGradientRho = configuration.sigmaFunctional(gradientRhoNorm,K);
                            obj.diffusionCoefficient(:, count) = obj.diffusionCoefficient(:, count) + sigmaGradientRho(:);
                            w1 = sigmaGradientRho.*w1;
                            w2 = sigmaGradientRho.*w2;
                            w3 = sigmaGradientRho.*w3;
                            du = [w1(:), w2(:), w3(:)];
                        case 'brightness'
                            [w2, w1, w3] = gradient(log(d + 2 * eps));
                            dC = configuration.sigmaFunctional(d);
                            obj.diffusionCoefficient(:, count) = obj.diffusionCoefficient(:, count) + dC(:);
                            w1 = dC.*w1;
                            w2 = dC.*w2;
                            w3 = dC.*w3;
                            du = [w1(:), w2(:), w3(:)];
                    end
                    obj.diffusiveVelocity(:, count) = obj.diffusiveVelocity(:, count) + du(:);
                    obj.diffusiveSpeed(:, count) = obj.diffusiveSpeed(:, count) + sqrt(sum(du.^2, 2));
                    u1 = v1 - w1;
                    u2 = v2 - w2;
                    u3 = v3 - w3;
                    
                    au = [u1(:), u2(:), u3(:)];
                    obj.augmentedVelocity(:, count) = obj.augmentedVelocity(:, count) + au(:);
                    obj.augmentedSpeed(:, count) = obj.augmentedSpeed(:, count) + sqrt(sum(au.^2, 2));  
                end                
                obj.speed(:, count) = obj.speed(:, count)/nt;
                obj.diffusiveVelocity(:, count) = obj.diffusiveVelocity(:, count)/nt;
                obj.diffusiveSpeed(:, count) = obj.diffusiveSpeed(:, count)/nt;
                obj.augmentedVelocity(:, count) = obj.augmentedVelocity(:, count)/nt;
                obj.augmentedSpeed(:, count) = obj.augmentedSpeed(:, count)/nt;
                obj.Pe(:, count) = obj.Pe(:, count) + obj.speed(:, count)./(obj.diffusiveSpeed(:, count) + eps);
                obj.diffusionCoefficient(:, count) = obj.diffusionCoefficient(:, count)/nt;
                obj.sFlux(:, count) = obj.sFlux(:,  count)/nt;
                count = count + 1;
            end
            %reshape each properties to (x,y,z,t) form
            scalarSize = [configuration.trueSize(1), configuration.trueSize(2), ...
                configuration.trueSize(3), obj.nData];
            vectorSize = [prod(configuration.trueSize), 3, obj.nData];
            obj.speed  = reshape(obj.speed, scalarSize);
            obj.diffusiveSpeed  = reshape(obj.diffusiveSpeed, scalarSize);
            obj.augmentedSpeed  = reshape(obj.augmentedSpeed, scalarSize);
            obj.Pe  = reshape(obj.Pe, scalarSize);
            obj.source = reshape(obj.source, scalarSize);
            obj.sFlux = reshape(obj.sFlux, scalarSize);
            obj.velocity  = reshape(obj.velocity, vectorSize);
            obj.diffusiveVelocity  = reshape(obj.diffusiveVelocity, vectorSize);
            obj.augmentedVelocity  = reshape(obj.augmentedVelocity, vectorSize);
            obj.diffusionCoefficient = reshape(obj.diffusionCoefficient, scalarSize);
        end
        
        function exportAverageProperty(obj)
            %EXPORTAVERAGEPROPERTY export average properties on time to
            %vtk format.
            configuration = Configuration3D(obj.tag);
            obj.initializeOutDir(3,configuration);
            obj.createOutDir(obj.pathOutput,obj.tag);
            obj.getPropertiesInMask;
            %speed
            obj.speedAverage = mean(obj.speed,2);
            vtkwrite(sprintf('%s/%s/%s_speedAverage_%s.vtk', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion), ...
                    'structured_points', 'mask', obj.speedAverage);
        end

        function exportToMat(obj)
            configuration = Configuration3D(obj.jsonFileName);
            %create output directory
            obj.initializeOutDir(3,configuration);
            obj.createOutDir(obj.pathOutput,obj.tag);

            save(sprintf('%s/eulerian.mat',obj.pathOutput),'obj');
        end

        function exportToVtk(obj)
            %EXPORTTOVTK save Eulerian specification of each property in VTK
            %format
            configuration = Configuration3D(obj.jsonFileName);
            obj.initializeOutDir(3,configuration);
            obj.createOutDir(obj.pathOutput,obj.tag);
            obj.getPropertiesInMask;
            vtkwrite(sprintf('%s/%s/%s_anato_%s.vtk', obj.pathOutput,...
                obj.outdir, obj.tag, obj.outversion), 'structured_points',...
                'mask', obj.anatomyImage);
            for i = 1:obj.nData
                
                vtkwrite(sprintf('%s/%s/%s_speed_%s_%d.vtk', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion, i), ...
                    'structured_points', 'mask', obj.speed(:,:,:,i));
                vtkwrite(sprintf('%s/%s/%s_Pe_%s_%d.vtk', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion, i), ...
                    'structured_points', 'mask', obj.Pe(:,:,:,i));
                
                if obj.addSource
                    vtkwrite(sprintf('%s/%s/%s_source_%s_%d.vtk', obj.pathOutput,...
                        obj.outdir, obj.tag, obj.outversion, i), ...
                        'structured_points', 'mask', obj.source(:,:,:,i));
                    vtkwrite(sprintf('%s/%s/%s_sourceFlux_%s_%d.vtk', obj.pathOutput,...
                        obj.outdir, obj.tag, obj.outversion, i), ...
                        'structured_points', 'mask', obj.sFlux(:,:,:,i));
                end
                %{
                vtkwrite(sprintf('%s/%s/%s_diffCoeff_%s_%d.vtk', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion, i), ...
                    'structured_points', 'mask', obj.diffusionCoefficient(:,:,:,i), 'precision', 6);
                %}
                [indX,indY,indZ] = ind2sub(configuration.n,obj.points);
                vtkwrite(sprintf('%s/%s/%s_velocity_%s_%d.vtk',...
                    obj.pathOutput, obj.outdir, obj.tag, obj.outversion, i),...
                    'structured_grid', indX(obj.indexROI)',indY(obj.indexROI)', indZ(obj.indexROI)',  ...
                    'vectors', 'vector_field', obj.velocity(obj.indexROI, 1, i),...
                    obj.velocity(obj.indexROI, 2, i), obj.velocity(obj.indexROI, 3, i));
            end
        end
        
        function exportToNifti(obj)
            %EXPORTTONIFTI save Eulerian specification of each property in
            %NIFTI format

            configuration = Configuration3D(obj.jsonFileName);
            ti = configuration.timeInitial;
            tj = configuration.timeJump;
            tf = configuration.timeFinal;
            obj.initializeOutDir(3,configuration);
            obj.createOutDir(obj.pathOutput,obj.tag);
            obj.getPropertiesInMask;
%            vtkwrite(sprintf('%s/%s/%s_anato_%s.vtk', obj.pathOutput,...
%                obj.outdir, obj.tag, obj.outversion), 'structured_points',...
%                'mask', obj.anatomyImage);
            % create Nifti form mask
            obj.experimentMaskNifti = obj.getNiftiMask(obj.pathMask,configuration.do_resize,...
                configuration.rawSize,configuration.sizeFactor);
  
            for i = 1:obj.nData
                speedNiftiForm = obj.applyNiftiMask(obj.experimentMaskNifti, configuration.do_resize,...
                    configuration.sizeFactor, configuration.xRange, configuration.yRange, configuration.zRange, obj.speed(:,:,:,i));
                save_untouch_nii(speedNiftiForm, sprintf('%s/%s/%s_speed_%s_E%d_%d_T_%d_%d.nii', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion, ti, tf+tj, ti+(i-1)*tj, ti+i*tj));
                PeNiftiForm = obj.applyNiftiMask(obj.experimentMaskNifti, configuration.do_resize,...
                    configuration.sizeFactor, configuration.xRange, configuration.yRange, configuration.zRange, obj.Pe(:,:,:,i));
                save_untouch_nii(PeNiftiForm, sprintf('%s/%s/%s_Pe_%s_E%d_%d_T_%d_%d.nii', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion, ti, tf+tj, ti+(i-1)*tj, ti+i*tj));
                
                if obj.addSource
                    srcNiftiForm = obj.applyNiftiMask(obj.experimentMaskNifti, configuration.do_resize,...
                        configuration.sizeFactor, configuration.xRange, configuration.yRange, configuration.zRange, obj.source(:,:,:,i));
                    save_untouch_nii(srcNiftiForm, sprintf('%s/%s/%s_source_%s_E%d_%d_T_%d_%d.nii', obj.pathOutput,...
                        obj.outdir, obj.tag, obj.outversion, ti, tf+tj, ti+(i-1)*tj, ti+i*tj));
                    sFluxNiftiForm = obj.applyNiftiMask(obj.experimentMaskNifti, configuration.do_resize,...
                        configuration.sizeFactor, configuration.xRange, configuration.yRange, configuration.zRange, obj.sFlux(:,:,:,i));
                    save_untouch_nii(sFluxNiftiForm, sprintf('%s/%s/%s_sourceFlux_%s_E%d_%d_T_%d_%d.nii', obj.pathOutput,...
                        obj.outdir, obj.tag, obj.outversion, ti, tf+tj, ti+(i-1)*tj, ti+i*tj));
                end
                %{
                vtkwrite(sprintf('%s/%s/%s_diffCoeff_%s_%d.vtk', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion, i), ...
                    'structured_points', 'mask', obj.diffusionCoefficient(:,:,:,i), 'precision', 6);
                [indX,indY,indZ] = ind2sub(configuration.n,obj.points);
                vtkwrite(sprintf('%s/%s/%s_velocity_%s_%d.vtk',...
                    obj.pathOutput, obj.outdir, obj.tag, obj.outversion, i),...
                    'structured_grid', indX(obj.indexROI)',indY(obj.indexROI)', indZ(obj.indexROI)',  ...
                    'vectors', 'vector_field', obj.velocity(obj.indexROI, 1, i),...
                    obj.velocity(obj.indexROI, 2, i), obj.velocity(obj.indexROI, 3, i));
               %}     
            end
        end
        
        function getPropertiesInMask(obj)
            configuration = Configuration3D(obj.jsonFileName);
            if ~isempty(configuration.pathAnatomyImage)
                obj.anatomyImage = load_untouch_nii(configuration.pathAnatomyImage);

                obj.anatomyImage = obj.anatomyImage.img(configuration.xRange,...
                                    configuration.yRange, configuration.zRange);
                if configuration.do_resize
                     obj.anatomyImage = resizeMatrix(obj.anatomyImage, ...
                        round(configuration.sizeFactor .* size(obj.anatomyImage)), 'linear');
                end    
            else
                anatomyImage = Mask(configuration.pathMask, configuration.isMaskFilled,...
                    configuration.xRange, configuration.yRange, configuration.zRange);
                if configuration.do_resize
                     anatomyImage = anatomyImage.resize(configuration.sizeFactor);
                end    
                if configuration.dilate > 0
                    anatomyImage = anatomyImage.dilate(configuration.dilate);
                end
                obj.anatomyImage = anatomyImage.contents;
            end
            if ~isempty(configuration.pathStartpointsMask)
                mask_brain = Mask(configuration.pathStartpointsMask, configuration.isMaskFilled,...
                    configuration.xRange, configuration.yRange, configuration.zRange, 1);
                if configuration.do_resize
                    mask_brain = mask_brain.resize(configuration.sizeFactor);
                end
                if configuration.dilate > 0
                    mask_brain = mask_brain.dilate(configuration.dilate);
                end
            else
                mask_brain = Mask(configuration.pathMask, configuration.isMaskFilled,...
                    configuration.xRange, configuration.yRange, configuration.zRange);
                if configuration.do_resize
                    mask_brain = mask_brain.resize(configuration.sizeFactor);
                end
                if configuration.dilate > 0
                    mask_brain = mask_brain.dilate(configuration.dilate);
                end
                if ~isempty(obj.anatomyImage)
                    obj.anatomyImage(mask_brain.contents == 0) = 0;
                end
            end
            obj.indexROI = find(mask_brain.contents(obj.points) == 1); 
            scalarMask = repmat(mask_brain.contents,[1 1 1 obj.nData]);
            vectorMask = repmat(mask_brain.contents(:),[1 3 obj.nData]);
            obj.speed(scalarMask == 0) = 0;
            obj.diffusiveSpeed(scalarMask == 0) = 0;
            obj.augmentedSpeed(scalarMask == 0) = 0;
            obj.Pe(scalarMask == 0) = 0;
            obj.source(scalarMask == 0) = 0;
            obj.sFlux(scalarMask == 0) = 0;
            obj.velocity(vectorMask == 0) = 0;
            obj.diffusiveVelocity(vectorMask == 0) = 0;
            obj.augmentedVelocity(vectorMask == 0) = 0;
            obj.diffusionCoefficient(scalarMask == 0) = 0;
        end
        
        function visualize(obj)
            obj.visualizeSpeed;
            configuration = Configuration3D(obj.jsonFileName);
            if configuration.addSource
                obj.visualizeR;
            end
        end

        function visualizeSpeed(obj)
            configuration = Configuration3D(obj.jsonFileName);
            % Eul speed
            x = 1:configuration.trueSize(1);
            y = 1:configuration.trueSize(2);
            z = 1:configuration.trueSize(3);
            speed1 = obj.speed(obj.speed>0);
            maxset = maxk(speed1(:),round(0.05*length(speed1(:))));
            spdMax = maxset(end);
            for i = 1:obj.nData
            figure,
            tmp = obj.speed(:,:,:,i);
            hs=slice(y,x,z,tmp,y,x,z); 
            set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
            alpha('color'),alphamap(linspace(0,1,100))
            grid off, box off, axis image
            xticks(0:10:obj.trueSize(1)); yticks(0:5:obj.trueSize(2)); zticks(0:5:obj.trueSize(3));
            ax = gca; ax.FontSize = 10; 
            xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
            set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
            set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
            colormap(jet)
            clim([0,spdMax])
            view([242.1011   14.4475])
            set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
            grid on,
            xlim([0 obj.trueSize(1)]); ylim([0 obj.trueSize(2)]); zlim([0 obj.trueSize(3)])
            
            cb = colorbar;
            cb.Ticks = linspace(0, 1, 6);
            cb.TickLabels = num2cell(round((0:spdMax/5:spdMax)*100)/100);
            cb.FontSize = 8;
            cb.Label.String = 'speed (a.u.)';
            title(sprintf('Eulerian Speed map E%d -> E%d',configuration.timeInitial+(i-1)*configuration.timeJump,configuration.timeInitial+i*configuration.timeJump));
            %text(1,-3,15,'speed (a.u.)','Rotation',90,'FontSize',25);
            
            %saveas(gcf, sprintf('%s/%s/%s_EulAveSpeed_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Eul,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
            end
        end

        function visualizeR(obj)
            %% VISUALZER visualize eulerian influx/clearance rate maps
            configuration = Configuration3D(obj.jsonFileName);
            % Eul r
            x = 1:configuration.trueSize(1);
            y = 1:configuration.trueSize(2);
            z = 1:configuration.trueSize(3);
            src1 = obj.source(abs(obj.source)>0);
            maxset = maxk(src1(:),round(0.05*length(src1(:))));
            srcMax = maxset(end);
            for i = 1:obj.nData
            figure,
            tmp = obj.source(:,:,:,i);
            tmp(tmp==0) = NaN;
            hs=slice(y,x,z,tmp,y,x,z); 
            set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
            alpha('color'),alphamap([linspace(0,0.12,100)])
            grid off, box off, axis image
            xticks(0:10:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
            ax = gca; ax.FontSize = 10; 
            xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
            set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
            set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
            clim([-srcMax,srcMax])
            colormap(bluewhitered_drdt)
            view([242.1011   14.4475])
            set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
            grid on,
            xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
            
            cb = colorbar;
            cb.Ticks = linspace(-srcMax, srcMax, 5);
            cb.TickLabels = num2cell(round((-srcMax:srcMax/2:srcMax)*100)/100);
            cb.FontSize = 8;
            cb.Label.String = 'r (a.u.)';
            title(sprintf('Eulerian influx/clearance rate map E%d -> E%d',configuration.timeInitial+(i-1)*configuration.timeJump,configuration.timeInitial+i*configuration.timeJump));
            %text(1,-3,12,'relative source','Rotation',90,'FontSize',25);
            
            %saveas(gcf, sprintf('%s/%s/%s_EulAveR_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Eul,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
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
        
        function niftiMask = getNiftiMask(path,do_resize,domainSize,sizeFactor)
            % GETNIFTIMASK create mask in nifti format which will used to
            % be saved
            niftiMask = load_untouch_nii(path);
            if do_resize
                imageSize = round(domainSize .* sizeFactor);
                niftiMask.img = zeros(imageSize);
                niftiMask.hdr.dime.dim(2:4) = imageSize;
            else
                imageSize = domainSize;
                niftiMask.img = zeros(imageSize);
            end
        end
        
        function quantityNiftiForm = applyNiftiMask(mask,do_resize,sizeFactor,xRange,yRange,zRange,quantity)
            % APPLYNIFTIMASK save quantity in nifti format mask
            quantityNiftiForm = mask;
            if do_resize
                quantityNiftiForm.img(round(xRange(1) * sizeFactor):round(xRange(end) * sizeFactor),...
                    round(yRange(1) * sizeFactor):round(yRange(end) * sizeFactor),...
                    round(zRange(1) * sizeFactor):round(zRange(end) * sizeFactor)) = quantity;
            else
                quantityNiftiForm.img(xRange,yRange,zRange) = quantity;
            end           
        end
        
    end
end

