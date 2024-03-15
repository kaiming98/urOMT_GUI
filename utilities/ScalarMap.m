classdef ScalarMap < FlowFieldCharacteristic
    % Speed maps visualize scalar values of velocity(Peclet, density, ...) 
    % magnitude across the domain.
    % 
    % Inputs for instance construction: 
    %   tag  - specifies dataset configuration and constants
    %   pathlines - pathlines of the flow, which is used for calculating 
    %               scalar map
    %   analysisConstants - will be used in the following calculation and
    %                       visualization
    
    properties (Access = public)
        cluster
        density % dimension: L^-3M
        speed % dimension: T^-1L
        diffusiveSpeed 
        augmentedSpeed 
        timeStep % dimension: T
        flux % dimension: T^-1L^-2M
        drho_dt % drho/dt, dimension T^-1L^-3M
        Pe %Peclet number (ratio of advective to diffusive transport)
        source % r
        sFlux % source flux: \rho*r
        diffusionCoefficient
        featureName
        metricNameShort
    end

    properties (Access = protected)
        pathlinesClustered
        experimentMaskNifti
        addSource
    end
    
    methods

        function obj = ScalarMap(tag, pathlines, analysisConstants)

            arguments
                tag %{mustBeText}
                pathlines
                analysisConstants
            end
            
            % inherit functionality from superclass
            configuration = Configuration3D(tag);
            obj = obj@FlowFieldCharacteristic(configuration);
            obj.analysisConstants = analysisConstants;
            obj.addSource = configuration.addSource;
            obj.jsonFileName = tag;
            % transfer important properties from pathlines to scalarMap
            % instance.
            obj.experimentMask = pathlines.experimentMask;
            obj.cutoffThresholds = pathlines.cutoffThresholds;
            obj.calculateProperties(pathlines,configuration);
        end
        
        function printInfo(obj,startPoints)
            %PRINTINFO print useful information for scalar map analysis
            configuration = Configuration3D(obj.jsonFileName);
            switch obj.analysisConstants.metricName
                case 'AveragePointwiseEuclideanMetric'
                    obj.metricNameShort = 'AvgPwEuc';
                    obj.featureName = sprintf('ResampleFeature(nb_points=%d)',...
                                        obj.analysisConstants.nBundles);
                case 'CosineMetric'
                    obj.metricNameShort = 'Cos';
                    obj.featureName = 'VectorOfEndpointsFeature()';
            end
            obj.initializeOutDir(1,configuration,startPoints);
            title_str = sprintf(['============= initiating...' ...
                '\n\nLagrangian-Pathline (%s data, mask = %d, affSmooth =' ...
                ' %d, dilate = %d), \nanalysis type = %s\n \nflowType = %s,' ...
                ' isInterpolated = %d, nTimeInterval = %d(%s), %s, spErode = %d, spDilate = %d,' ...
                ' nEulerianSteps= %d, \ncutoffStr = %s, concThresh = %5.4f, ' ...
                ' speedThresh = %5.4f, imBoundFlag0 = %d, %s, slTol = %d, \nclusters, ' ...
                ' centroid = %d, nBundles = %d, metric = %s, qbthresh = %d, ' ...
                ' clusterTol = %d, clusterCutoff = %d, \ndiff = %s, tj = %d, nt= %d%s_%s_%s\n\n'], ...
                obj.tag, configuration.mask_number, configuration.smooth, configuration.dilate, ...
                'scalarmap', obj.analysisConstants.flowType, obj.analysisConstants.isInterpolated, ...
                obj.analysisConstants.nTimeInterval, obj.analysisConstants.distanceFlag, ...
                startPoints.maskString, startPoints.erodeMask, startPoints.dilateMask, ...
                obj.analysisConstants.nEulerianSteps, obj.analysisConstants.cutoffFlag,...
                obj.cutoffThresholds.density, obj.cutoffThresholds.speed, obj.analysisConstants.intensityRangeFlag,...
                startPoints.methodTitle, obj.streamlineLengthTol, obj.analysisConstants.centroidFlag,...
                obj.analysisConstants.nBundles, obj.metricNameShort, obj.analysisConstants.quickBundleThreshold,...
                obj.analysisConstants.clusterTol, obj.analysisConstants.clusterCutoff,...
                obj.sigmaStr, configuration.timeJump, configuration.nt,...
                startPoints.thresholdString, obj.paperFigStr, obj.dateStr);
            fprintf(title_str)
        end
                
        function exportToNifti(obj)
            %EXPORTTONIFTI Saves velocity and Peclet fields in NIfTI format.
            %   Neuroimaging Informatics Technology Initiative (NIfTI) is a
            %   type of file format commonly used in imaging informatics, 
            %   in contrast to DICOM which is the clinical standard.
            
            configuration = Configuration3D(obj.jsonFileName);
            % create output directory
            obj.createOutDir(obj.pathOutput,obj.tag);
            
            % create Nifti form mask
            obj.experimentMaskNifti = obj.getNiftiMask(obj.pathMask,configuration.do_resize,...
                configuration.rawSize,configuration.sizeFactor);

            speedNiftiForm = obj.applyNiftiMask(obj.experimentMaskNifti, configuration.do_resize,...
                configuration.sizeFactor, configuration.xRange, configuration.yRange, configuration.zRange, obj.speed);
            PeNiftiForm = obj.applyNiftiMask(obj.experimentMaskNifti, configuration.do_resize,...
                configuration.sizeFactor, configuration.xRange, configuration.yRange, configuration.zRange, obj.Pe);    
            if obj.addSource
                sourceNiftiForm = obj.applyNiftiMask(obj.experimentMaskNifti, configuration.do_resize,...
                    configuration.sizeFactor, configuration.xRange, configuration.yRange, configuration.zRange, obj.source); 
            end
%             lnMask = Mask(obj.pathMask,1);
%             obj.speed(lnMask.contents==0) = 0;
%             vtkwrite(sprintf('%s/%s/%s_LagSpeed_E%02d_%02d_%s_%s.vtk',...
%                 obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
%                 obj.timeFinal + obj.timeJump, obj.paperFigStr, obj.dateStr), 'structured_points',...
%                     'mask', obj.speed);
            % save to nifty (view in Amira later)
            save_untouch_nii(speedNiftiForm, sprintf('%s/%s/%s_LagSpeed_E%02d_%02d_%s_%s.nii',...
                obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
                obj.timeFinal + obj.timeJump, obj.paperFigStr, obj.dateStr));
            save_untouch_nii(PeNiftiForm, sprintf('%s/%s/%s_LagPe_E%02d_%02d_%s_%s.nii',...
                obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
                obj.timeFinal + obj.timeJump, obj.paperFigStr, obj.dateStr));
            if obj.addSource
                save_untouch_nii(sourceNiftiForm, sprintf('%s/%s/%s_LagSource_E%02d_%02d_%s_%s.nii',...
                    obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
                    obj.timeFinal + obj.timeJump, obj.paperFigStr, obj.dateStr));
                fprintf('Source Map in nifty format saved in %s/%s\n\n',...
                    obj.pathOutput, obj.outdir)      
            end
            
            fprintf('Speed and Peclet Map in nifty format saved in %s/%s\n\n',...
                obj.pathOutput, obj.outdir)      
        end    
        
        function visualize3D(obj)
            %VISUALIZE3D plot speed scalar map and peclet scalar map
            x = 1:obj.trueSize(1);
            y = 1:obj.trueSize(2);
            z = 1:obj.trueSize(3);
            
            figure,
            hsMask = slice(y, x, z, 0.03*max(obj.speed,[],'all')*obj.experimentMask.contents, x, y, z);
            set(hsMask, 'EdgeColor', 'none', 'FaceColor', '#ffcccb', 'FaceAlpha', 0.04);
            hold on
            hs = slice(y, x, z, obj.speed, x, y, z);
            set(hs, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.04);
            alpha('color'), alphamap(linspace(0, 1, 100))
            title(sprintf('Test: tag = %s, Speed Map', obj.tag), 'Fontsize', 20)
            box off
            axis image
            xlabel('x-axis')
            ylabel('y-axis')
            zlabel('z-axis')
            colormap(jet)
            set(gca, 'Color', [0.85, 0.85, 0.93]), set(gcf, 'unit',...
                'normalized', 'position', [0.1, 1, 0.4, 0.5], 'Color',...
                [0.85, 0.85, 0.93]), set(gcf, 'InvertHardcopy', 'off')
            colorbar, grid on,
            xlim([0 obj.trueSize(2)])
            ylim([0 obj.trueSize(1)])
            zlim([0 obj.trueSize(3)])
        
            saveas(gcf, sprintf('%s/%s/%s_LagSpeed3D_E%02d_%02d.png',...
                   obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
                   obj.timeFinal + obj.timeJump));
        
            figure,
            hsMask = slice(y, x, z, 0.03*max(obj.Pe,[],'all')*obj.experimentMask.contents, x, y, z);
            set(hsMask, 'EdgeColor', 'none', 'FaceColor','#ffcccb', 'FaceAlpha', 0.04);
            hold on
            hs = slice(y, x, z, obj.Pe, x, y, z);
            set(hs, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.04);
            alpha('color'), alphamap(linspace(0, 1, 100))
            title(sprintf('Test: tag = %s, Peclet Map', obj.tag), 'Fontsize', 20)
            box off
            axis image
            xlabel('x-axis')
            ylabel('y-axis')
            zlabel('z-axis')
            colormap(jet)
            set(gca, 'Color', [0.85, 0.85, 0.93]), set(gcf, 'unit',...
                'normalized', 'position', [0.1, 1, 0.4, 0.5], 'Color',...
                [0.85, 0.85, 0.93]), set(gcf, 'InvertHardcopy', 'off')
            colorbar
            xlim([0 obj.trueSize(2)])
            ylim([0 obj.trueSize(1)])
            zlim([0 obj.trueSize(3)])
            
            saveas(gcf, sprintf('%s/%s/%s_LagPe3D_E%02d_%02d.png',...
                   obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
                   obj.timeFinal + obj.timeJump)); 
            if obj.addSource
                % calculate alphamap
                maxSource = max(obj.source(:));
                minSource = min(obj.source(:));
                alpha_map = abs(linspace(minSource,maxSource,100))/(max(abs(maxSource),abs(minSource)));
                % plot result
                figure,
                hsMask = slice(y, x, z, 0.03*max(obj.source,[],'all')*obj.experimentMask.contents, x, y, z);
                set(hsMask, 'EdgeColor', 'none', 'FaceColor','#ffcccb', 'FaceAlpha', 0.04);
                hold on
                hs = slice(y, x, z, obj.source, x, y, z);
                set(hs, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.04);
                alpha('color'), alphamap(alpha_map)
                title(sprintf('Test: tag = %s, Source Map', obj.tag), 'Fontsize', 20)
                box off
                axis image
                xlabel('x-axis')
                ylabel('y-axis')
                zlabel('z-axis')
                colormap(jet)
                set(gca, 'Color', [0.85, 0.85, 0.93]), set(gcf, 'unit',...
                    'normalized', 'position', [0.1, 1, 0.4, 0.5], 'Color',...
                    [0.85, 0.85, 0.93]), set(gcf, 'InvertHardcopy', 'off')
                colorbar
                xlim([0 obj.trueSize(2)])
                ylim([0 obj.trueSize(1)])
                zlim([0 obj.trueSize(3)])

                saveas(gcf, sprintf('%s/%s/%s_LagSource3D_E%02d_%02d.png',...
                       obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
                       obj.timeFinal + obj.timeJump)); 
            end
        end
        
        function visualize(obj)
            %VISUALIZE2D Plots the scalar contour maps 
            % Contours are plotted on 3 equidistant slices spanning the 
            % Z directon.

            volume = {obj.speed, obj.Pe};
            dimension = size(obj.speed);
            variableNames = {'Velocity','Peclet'};

            for j = 1:length(volume)
                arraySlices = floor(linspace(1,dimension(3),4));
                figure(j); hold on
                phandles = contourslice(volume{j},[],[],arraySlices,100);
                view(3); axis tight
                set(phandles, 'LineWidth',2)
                mhandles = contourslice(obj.experimentMask.contents,[],[],arraySlices,1);
                set(mhandles,'EdgeColor','k','linewidth',3);
                colorbar;
                xlabel('x axis')
                ylabel('y axis')
                zlabel('z axis')
                title(sprintf('%s - %s',variableNames{j},obj.tag));
                set(findall(gcf,'-property','FontSize'),'FontSize',16)
                grid on
            end
        end

        function visualize2D(obj)
            
            %VISUALIZE2D plot speed scalar map and peclet scalar map slice
            % Plots largest XY slice by default
            % 'Slice' specifies 
            [sliceYZ,sliceXZ,sliceXY] = findLargestSlice(obj.speed);
            speedXY =  squeeze(obj.speed(:,:,sliceXY));
            speedXY  = speedXY';
            ind = find(speedXY>0);
            [x,y] = ind2sub(size(speedXY),ind);
            mask2D = squeeze(obj.experimentMask.contents(:,:,sliceXY));
            convolutionFilter = [0 -1/4 0;-1/4 1 -1/4;0 -1/4 0];
            maskConvolved = conv2(convolutionFilter,mask2D);
            maskBoundary = maskConvolved(2:end-1,2:end-1);
            indB =find(maskBoundary'>0);
            [xB,yB] = ind2sub(size(maskBoundary'),indB);
            
            % speed map 2D visualization
            figure,
            scatter(xB,yB,15,'ks','fill');
            hold on;
            scatter(x,y,15,speedXY(ind),'s','fill')
            colormap(jet)
            colorbar
            axis image
            title(sprintf('Test: tag = %s, Speed Map', obj.tag), 'Fontsize', 20)
            xlabel('x-axis')
            ylabel('y-axis')
            xlim([0 obj.trueSize(2)])
            ylim([0 obj.trueSize(1)])
            
            saveas(gcf, sprintf('%s/%s/%s_LagSpeed2D_E%02d_%02d.png',...
                   obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
                   obj.timeFinal + obj.timeJump));
            
            % Peclet map 2D visualization
            PeXY =  squeeze(obj.Pe(:,:,sliceXY));
            PeXY  = PeXY';
            ind = find(PeXY>0);
            [x,y] = ind2sub(size(PeXY),ind);
            figure,
            scatter(xB,yB,15,'ks','fill');
            hold on;
            scatter(x,y,15,PeXY(ind),'s','fill')
            colormap(jet)
            colorbar
            axis image
            title(sprintf('Test: tag = %s, Peclet Map', obj.tag), 'Fontsize', 20)
            xlabel('x-axis')
            ylabel('y-axis')
            xlim([0 obj.trueSize(2)])
            ylim([0 obj.trueSize(1)])
            
            saveas(gcf, sprintf('%s/%s/%s_LagPe2D_E%02d_%02d.png',...
                   obj.pathOutput, obj.outdir, obj.tag, obj.timeInitial,...
                   obj.timeFinal + obj.timeJump));
        end

        function calculateProperties(obj,pathlines,configuration)
            %CALCULATEPROPERTIES calculate scalar properties defined over 
            % the spatial domain based on the N largest clusters of pathlines.

            pathlines.centroid = double(pathlines.centroid);
            obj.pathlinesClustered = cell(1, pathlines.nClusters);     
            for ind = 1:pathlines.nClusters
                obj.pathlinesClustered{1, ind} = squeeze(pathlines.centroid(ind, :, :));
            end
        
            [~, clusterOrder] = sort(squeeze(pathlines.centroid(:, 1, 1))); %sl_centroid is nclus x npoints x 3
            
            %spatial grid size
            h1 = 1;
            h2 = 1;
            h3 = 1;
            
            % initialize temporary masks:
            obj.initializeTmpMap(configuration.trueSize);

            % make cluster masks
            % getting clustered pathline start points
            startPointTmp = zeros(configuration.trueSize);
        
            for iCluster = 1:pathlines.nClusters
                slines_tmp = pathlines.qualified([pathlines.index{clusterOrder(iCluster)}] + 1);
                rlines_tmp = pathlines.rho([pathlines.index{clusterOrder(iCluster)}] + 1);
                spdlines_tmp = pathlines.speed([pathlines.index{clusterOrder(iCluster)}] + 1);
                dspdlines_tmp = pathlines.diffusiveSpeed([pathlines.index{clusterOrder(iCluster)}] + 1);
                aspdlines_tmp = pathlines.augmentedSpeed([pathlines.index{clusterOrder(iCluster)}] + 1);
                tlines_tmp = pathlines.timeStep([pathlines.index{clusterOrder(iCluster)}] + 1);
                flines_tmp = pathlines.flux([pathlines.index{clusterOrder(iCluster)}] + 1);
                dlines_tmp = pathlines.drho_dt([pathlines.index{clusterOrder(iCluster)}] + 1);
                dCoefflines_tmp = pathlines.diffusionCoefficient([pathlines.index{clusterOrder(iCluster)}] + 1);
                if obj.addSource
                    srlines_tmp = pathlines.source([pathlines.index{clusterOrder(iCluster)}] + 1);
                    srFluxlines_tmp = pathlines.sFlux([pathlines.index{clusterOrder(iCluster)}] + 1);
                end
                
                %getting start points in world coords
                spx = cellfun(@(x) x(1, 1), pathlines.qualified([pathlines.index{clusterOrder(iCluster)}] + 1));
                spy = cellfun(@(x) x(1, 2), pathlines.qualified([pathlines.index{clusterOrder(iCluster)}] + 1));
                spz = cellfun(@(x) x(1, 3), pathlines.qualified([pathlines.index{clusterOrder(iCluster)}] + 1));
                %getting start points in matlab coords
                spm1 = round(spy ./ h1 + 0.5);
                spm2 = round(spx ./ h2 + 0.5);
                spm3 = round(spz ./ h3 + 0.5);
                indtmp = sub2ind(configuration.trueSize, spm1, spm2, spm3);
                startPointTmp(indtmp) = 1;
        
                n_slines = size(slines_tmp, 2);
        
                for ind_line = 1:n_slines
                    %convert back to MATLAB grid
                    sline = round((slines_tmp{ind_line} ./ [h1, h2, h3]) + 0.5);
                    [sline, ia, ~] = unique(sline, 'rows', 'stable');
                    rsl = rlines_tmp{ind_line}(ia);
                    ssl = spdlines_tmp{ind_line}(ia);
                    dssl = dspdlines_tmp{ind_line}(ia);
                    assl = aspdlines_tmp{ind_line}(ia);
                    tsl = tlines_tmp{ind_line}(ia);
                    fsl = flines_tmp{ind_line}(ia);
                    dsl = dlines_tmp{ind_line}(ia);
                    dCoeffsl = dCoefflines_tmp{ind_line}(ia);
        
                    subs_1 = sline(:, 2);
                    subs_2 = sline(:, 1);
                    subs_3 = sline(:, 3);
                    subs_1(subs_1 < 1) = 1;
                    subs_2(subs_2 < 1) = 1;
                    subs_3(subs_3 < 1) = 1;
                    subs_1(subs_1 > configuration.trueSize(1)) = configuration.trueSize(1);
                    subs_2(subs_2 > configuration.trueSize(2)) = configuration.trueSize(2);
                    subs_3(subs_3 > configuration.trueSize(3)) = configuration.trueSize(3);
                    inds = sub2ind(configuration.trueSize, subs_1, subs_2, subs_3);
                    
                    if obj.addSource
                        srl = srlines_tmp{ind_line}(ia);
                        obj.source(inds) = srl;
                        srFluxl = srFluxlines_tmp{ind_line}(ia);
                        obj.sFlux(inds) = srFluxl;
                    end
                    
                    obj.cluster(inds) = iCluster * 1;
                    obj.density(inds) = rsl;
                    obj.speed(inds) = ssl;
                    obj.diffusiveSpeed(inds) = dssl;
                    obj.augmentedSpeed(inds) = assl;
                    obj.timeStep(inds) = tsl;
                    obj.flux(inds) = fsl;
                    obj.drho_dt(inds) = dsl;
                    obj.diffusionCoefficient(inds) = dCoeffsl;
                end
        
                %centroid mask:
                if obj.analysisConstants.centroidFlag
                    obj.getCentroidMask(configuration.trueSize,iCluster);
                end
        
            end
            
            % calculate Peclet number
            obj.Pe = obj.speed ./ obj.diffusiveSpeed;
            obj.Pe(isnan(obj.Pe)) = 0;
            obj.Pe(obj.Pe == Inf) = 0;
            
            if obj.analysisConstants.maskFlag
                obj.applyMask(pathlines.experimentMask);
            end
            % only save within brain mask
            if ~isempty(configuration.pathStartpointsMask)
                mask_brain = Mask(configuration.pathStartpointsMask,configuration.isMaskFilled,...
                    configuration.xRange, configuration.yRange, configuration.zRange, 1);        
                obj.applyMask(mask_brain);
            end
        end

    end
    
    methods (Access = protected)
        
        function applyMask(obj,mask)
            %APPLYMASK make sure nothing outside of masked region
            obj.cluster(~mask.contents) = 0;
            obj.density(~mask.contents) = 0;
            obj.speed(~mask.contents) = 0;
            obj.diffusiveSpeed(~mask.contents) = 0;
            obj.augmentedSpeed(~mask.contents) = 0;
            obj.timeStep(~mask.contents) = 0;
            obj.flux(~mask.contents) = 0;
            obj.drho_dt(~mask.contents) = 0; 
            obj.Pe(~mask.contents) = 0;
            obj.diffusionCoefficient(~mask.contents) = 0;
            if obj.addSource
                obj.source(~mask.contents) = 0;
                obj.sFlux(~mask.contents) = 0;
            end
        end
        
        function getCentroidMask(obj,trueSize,iCluster)
            centroid_tmp = round(obj.pathlinesClustered{I(iCluster)});
            centroid_tmp = unique(centroid_tmp, 'rows', 'stable');
            subs_1 = centroid_tmp(:, 2);
            subs_2 = centroid_tmp(:, 1);
            subs_3 = centroid_tmp(:, 3);
            subs_1(subs_1 < 1) = 1;
            subs_2(subs_2 < 1) = 1;
            subs_3(subs_3 < 1) = 1;
            subs_1(subs_1 > trueSize(1)) = trueSize(1);
            subs_2(subs_2 > trueSize(2)) = trueSize(2);
            subs_3(subs_3 > trueSize(3)) = trueSize(3);
            inds = sub2ind(trueSize, subs_1, subs_2, subs_3);
            cent_tmp(inds) = iCluster * 1;
        end
        
        function initializeTmpMap(obj,trueSize)
            % INITIALIZETMPMAP create zero scalar maps with size trueSize
            obj.cluster = zeros(trueSize); 
            obj.density = zeros(trueSize);
            obj.speed = zeros(trueSize);
            obj.diffusiveSpeed = zeros(trueSize);
            obj.augmentedSpeed = zeros(trueSize);
            obj.timeStep = zeros(trueSize);
            obj.flux = zeros(trueSize);
            obj.drho_dt = zeros(trueSize);
            obj.diffusionCoefficient = zeros(trueSize);
            if obj.addSource
                obj.source = zeros(trueSize);
                obj.sFlux = zeros(trueSize);
            end
        end
        

        
    end
    
    methods (Static)
        
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
                niftiMask.hdr.dime.dim(2:4) = imageSize;
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

