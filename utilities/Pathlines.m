classdef Pathlines < FlowFieldCharacteristic
    % class of pathlines, containing pathlines, 
    % (rho, speed, diffusive speed, augmented speed, time step, flux, drho/dt,peclet) 
    % along the pathlines.   
    %
    % Inputs for instance construction: 
    %   tag - specifies dataset configuration and constants
    %   lagrangianSpecification - gives the properties information that
    %   pathlines calculation needs
    
    properties
        position % dimension: L
        rho %density; dimension: L^-3M
        speed % dimension: T^-1L
        diffusiveSpeed 
        augmentedSpeed 
        timeStep % dimension: T
        flux % dimension: T^-1L^-2M
        drho_dt %drho/dt; dimension T^-1L^-3M
        drho_dx %|\nabla rho|
        Pe
        source % r
        sFlux % source flux r*\rho
        diffusionCoefficient
        qualified %qualified pathlines
        PeInROI % peclet in anatomical region of interest
        speedInROI % speed in anatomical region of interest
        drho_dxInROI % |\nabla rho| in anatomical region of interest
        positionInROI % pathlines in anatomical region of interest
        sourceInROI % source in anatomical region of interest
        sFluxInROI % source flux in anatomical region of interest
        diffusionCoeffcientInROI
        startPointsInROI
        endPointsInROI
        displacementInROI
        directionInROI
        anatomyImage
        pathlineJump % how much to downsample pathlines when exporting to vtk format 
        centroid
        nClusters
        index
    end
    
    properties (Access = protected)
        addSource
        indexROI
        step %temporary time step
        nPoints
        length
        displacement
        normalizedDisplacement
        startPoints
        endPoints
        direction
        displacementLength
    end
    
    methods

        function obj = Pathlines(tag,lagrangianSpecification)

            arguments
                tag 
                lagrangianSpecification
            end
            
            configuration = Configuration3D(tag);
            obj = obj@FlowFieldCharacteristic(configuration);
            obj.jsonFileName = tag;
            % transfer important properties from lagrangianSpecification to
            % this pathlines instance.
            obj.experimentMask = lagrangianSpecification.experimentMask;
            obj.analysisConstants = lagrangianSpecification.analysisConstants;
            obj.cutoffThresholds = lagrangianSpecification.cutoffThresholds;
            obj.addSource = configuration.addSource;
            % set useful variables(starting points and voxel spatial size)
            sx = lagrangianSpecification.startPoints.x;
            sy = lagrangianSpecification.startPoints.y;
            sz = lagrangianSpecification.startPoints.z;
            h1 = 1; h2 = 1; h3 = 1;
            %initialize streamlines:
            nstartpoints = length(sx);
            obj.initialize(nstartpoints);
            %convert from matlab grid to cell-centered grid:
            s1 = (sy - 0.5) .* h1; %i/y-axis
            s2 = (sx - 0.5) .* h2; %j/x-axis
            s3 = (sz - 0.5) .* h3; %k/z-axis
            sp_123 = [s1, s2, s3];
            %current point i.e. list of current location in each
            %streamline that hasn't been terminated
            pcur = sp_123; 
            %keep track of the # of streamlines that have not yet been terminated
            npoints = length(pcur); 
            %% select qualified pathlines
            iPathlines = 0; %counter for added pathlines
            iPathlines = obj.update(lagrangianSpecification,npoints,obj.analysisConstants.nPathlines,iPathlines);
            obj.remove0(iPathlines);
            obj.doThresholding(obj.streamlineLengthTol);
            fprintf([' # of start points = %d\n # of effective pathlines ' ...
                'after pathline-number (pln) threshold = %d \n # of ' ...
                'effective pathlines after Euclidean dist (sl_tol)' ...
                'threshold = %d\n'], npoints, iPathlines, length(obj.position))
            obj.qualified = cellfun(@(x) x(:, [2, 1, 3]), obj.position, 'UniformOutput', false);
        end
        
        function doQuickBundle(obj,pathToPythonInterpreter)
            %DOQUICKBUNDLE Clustering algorithm for reducing number of
            %pathlines
            %
            % Large numbers of pathlines are hard to visualize, interact with, 
            % and interpret. The QuickBundles (QB) algorithm reduces complexity
            % by sequentially assigning each pathlines to its closest bundle. 
            % If for a given pathline its closest bundle is farther than 
            % a given threshold, a new bundle is created and the pathline is assigned to it.
            %
            % The python script called from this function outputs and saves
            % the following 2 matrices in matlab format: 
            %   'pl_centroid_array.mat' a summary of the clustered
            %   pathlines 
            %   'pli_array.mat' an array of the indices of each cluster of pathlines 
    
            obj.deleteOldInfo;
            pl_cur = obj.qualified;            
            save('pl_cur.mat', 'pl_cur');
            
            oldpath = getenv('PATH');
            setenv('PATH',pathToPythonInterpreter);
            !python run_dipyQB_pl.py
            
            setenv('PATH',oldpath);

        end
        
        function filterClusters(obj)
            %FILTERCLUSTERS filter the N(clusterCutoff) largest clusters.
            
            % add start points of pathlines considered for clustering analysis
            load(fullfile(cd, 'pli_array.mat'));
            load(fullfile(cd, 'pl_centroid_array.mat'));

            % only want clusters with more than tolerant # of pathlines:
            clusterLength = cellfun('size', pli_array, 2)';
            fprintf(' # of original clusters = %d\n', length(clusterLength))
        
            obj.index = pli_array(clusterLength > obj.analysisConstants.clusterTol);
            obj.centroid = pl_centroid_array(clusterLength > obj.analysisConstants.clusterTol, :, :);
            obj.nClusters = size(obj.index, 2); %# of clusters that made the cutoff
            clusterLength = cellfun('size', obj.index, 2)';
            fprintf(' # of clusters after cluster-length (clus_tol) threshold = %d\n',...
                    length(clusterLength))
            [~, rankClusterSize] = sort(clusterLength);
        
            %only keep largest N clusters where N = clusterCutoff:
            if obj.analysisConstants.clusterCutoff > 0       
                if obj.nClusters > obj.analysisConstants.clusterCutoff
                    ind_tmp = zeros(obj.nClusters, 1);
                    ind_tmp(rankClusterSize(end - obj.analysisConstants.clusterCutoff + 1:end)) = 1;
                    % to keep clusters in same order that they were returned 
                    % from quick bundle algorithm:
                    obj.index = obj.index(ind_tmp == 1);
                    obj.centroid = obj.centroid(ind_tmp == 1, :, :);
                    obj.nClusters = size(obj.index, 2); %# of clusters that made the cutoff
                    fprintf('nclus=%d\n', obj.nClusters)
                end  
            end
            
            fprintf(' # of clusters after max-cluster-number (clus_cutoff) threshold = %d\n',...
                    obj.nClusters)
        end
        
        function printInfo(obj,startPoints)
            %PRINTINFO print useful information for pathlines analysis
            configuration = Configuration3D(obj.jsonFileName);
            obj.initializeOutDir(2,configuration,startPoints);           

            title_str = sprintf(['============= initiating...'...
                '\n\nLagrangian-Pathline(%s data, mask = %d, affSmooth ='...
                '%d, dilate = %d), \nanalysis type = %s\n\nflowType = %s,'...
                'isInterpolated = %d, nTimeInterval = %d(%s), %s, spErode = %d, spDilate = %d,'...
                'nEulerianSteps= %d, \ncutoffStr = %s, concThresh = %5.4f,'...
                'speedThresh = %5.4f, imBoundFlag = %d, %s, slTol = %d, '...
                '\ndiff = %s, tj = %d, nt = %d%s_%s_%s\n\n'], ...
                obj.tag, configuration.mask_number, configuration.smooth,...
                configuration.dilate, 'vectors', obj.analysisConstants.flowType,...
                obj.analysisConstants.isInterpolated, obj.analysisConstants.nTimeInterval,...
                obj.analysisConstants.distanceFlag, startPoints.maskString, startPoints.erodeMask,...
                startPoints.dilateMask, obj.analysisConstants.nEulerianSteps,...
                obj.analysisConstants.cutoffFlag, obj.cutoffThresholds.density, obj.cutoffThresholds.speed,...
                obj.analysisConstants.intensityRangeFlag, startPoints.methodTitle,...
                obj.streamlineLengthTol, obj.sigmaStr, obj.timeJump, obj.nt,...
                startPoints.thresholdString, obj.paperFigStr, obj.dateStr);
            fprintf(title_str)

        end

        function checkOutputDir(obj,outputPath,tag)
            %CHECKOUTPUTDIR check if the directory for output is exist or
            %not
             if ~exist(sprintf('%s/%s',outputPath, obj.outdir), 'dir')
                fprintf('%s/%s\n Directory does not exist :(\n', outputPath, obj.outdir);
                return
             else
                fprintf('%s: %s Directory exists :)\n', tag, obj.outdir);
             end
        end
        
        function calculateExtraProperties(obj)
            %calculate extra properties of pathlines, including number of
            %points, length, displacement field...
            
            %set spatial grid size
            h1 = 1; h2 = 1; h3 = 1;
            
            fprintf('Total original %d pathlines\n', length(obj.position));

            %convert from cell-centered grid to matlab grid:
            obj.position = cellfun(@(x) [x(:, 1) ./ h1 + 0.5, x(:, 2) ./ h2 + 0.5,...
                x(:, 3) ./ h3 + 0.5], obj.position, 'UniformOutput', false);

            dispnor = cellfun(@(x) (x(end, :) - x(1, :)) / norm(x(end, :) - x(1, :)),...
                                    obj.position, 'UniformOutput', false);
            disp = cellfun(@(x) x(end, :) - x(1, :), obj.position, 'UniformOutput', false);
            directionAtStart = cellfun(@(x) (x(2, :) - x(1, :))/norm(x(2, :) - x(1, :)), obj.position, 'UniformOutput', false);
            startp = cellfun(@(x) x(1, :), obj.position, 'UniformOutput', false);
            endp = cellfun(@(x) x(end, :), obj.position, 'UniformOutput', false);

            obj.nPoints = cellfun(@(x) size(x, 1), obj.position); % number of points in each pathline
            obj.length = cellfun(@(x) sum(sqrt(sum(diff(x).^2, 2))), obj.position); % total length of path in each pathline
            obj.displacement = reshape([disp{:}]', 3, [])'; % displacement field
            obj.direction = reshape([directionAtStart{:}]',3,[])'; % direction of pathlines at startPoint
            obj.normalizedDisplacement = reshape([dispnor{:}]', 3, [])'; % normalized displacement field
            obj.startPoints = reshape([startp{:}]', 3, [])'; % start points
            obj.endPoints = reshape([endp{:}]', 3, [])'; % end points
            obj.displacementLength = sqrt(obj.displacement(:, 1).^2 +...
                obj.displacement(:, 2).^2 + obj.displacement(:, 3).^2);% displacement length in each pathline
        end
        
        function getPropertiesInMask(obj)
            % GETPROPERTIESINMASK mask sure pathlines and properties are 
            % within the anatomical region of interest.
            
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
                if ~isempty(obj.anatomyImage)
                    obj.anatomyImage(mask_brain.contents == 0) = 0;
                end
                % index of those starting within brain
                obj.indexROI = find(mask_brain.contents(sub2ind(obj.trueSize, obj.startPoints(:, 1),...
                    round(obj.startPoints(:, 2)), obj.startPoints(:, 3))) == 1); 


                obj.positionInROI = obj.position(obj.indexROI);
                obj.PeInROI = obj.Pe(obj.indexROI);
                obj.speedInROI = obj.speed(obj.indexROI);
                obj.drho_dxInROI = obj.drho_dx(obj.indexROI);
                obj.startPointsInROI = obj.startPoints(obj.indexROI,:);
                obj.endPointsInROI = obj.endPoints(obj.indexROI,:);
                obj.displacementInROI = obj.displacement(obj.indexROI,:);
                obj.directionInROI = obj.direction(obj.indexROI,:);
                obj.diffusionCoeffcientInROI = obj.diffusionCoefficient(obj.indexROI);
                if obj.addSource
                    obj.sourceInROI = obj.source(obj.indexROI);
                end
            else
                obj.positionInROI = obj.position;
                obj.PeInROI = obj.Pe;
                obj.speedInROI = obj.speed;
                obj.drho_dxInROI = obj.drho_dx;
                obj.directionInROI = obj.direction;
                obj.startPointsInROI = obj.startPoints;
                obj.endPointsInROI = obj.endPoints;
                obj.displacementInROI = obj.displacement;
                obj.diffusionCoeffcientInROI = obj.diffusionCoefficient;
                if obj.addSource
                    obj.sourceInROI = obj.source;
                end
            end
        end
        
        function exportToMat(obj)
            %create output directory
            obj.createOutDir(obj.pathOutput,obj.tag);
            obj.checkOutputDir(obj.pathOutput,obj.tag); 
            
            %configuration = Configuration3D(obj.jsonFileName);
            save(sprintf('%s/pathlines.mat',obj.pathOutput),'obj');
        end

        function exportToVtk(obj)
            %EXPORTTOVTK save the pathlines and propties in vtk format  

            %create output directory
            obj.createOutDir(obj.pathOutput,obj.tag);
            obj.checkOutputDir(obj.pathOutput,obj.tag); 
            
            configuration = Configuration3D(obj.jsonFileName);
            
            % write out to vtk
            vtkwrite_pathlines(sprintf('%s/%s/%s_pathlines_lentol_%.1f_jp_%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                obj.positionInROI(1:obj.analysisConstants.pathlineJump:end));
            vtkwrite_spdlines(sprintf('%s/%s/%s_Pelines_lentol_%.1f_jp_%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                obj.positionInROI(1:obj.analysisConstants.pathlineJump:end),...
                obj.PeInROI(1:obj.analysisConstants.pathlineJump:end));
            vtkwrite_spdlines(sprintf('%s/%s/%s_drhodxlines_lentol_%.1f_jp_%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                obj.positionInROI(1:obj.analysisConstants.pathlineJump:end),...
                obj.drho_dxInROI(1:obj.analysisConstants.pathlineJump:end));
            vtkwrite_spdlines(sprintf('%s/%s/%s_Spdlines_lentol_%.1f_jp_%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                obj.positionInROI(1:obj.analysisConstants.pathlineJump:end),...
                obj.speedInROI(1:obj.analysisConstants.pathlineJump:end));
            vtkwrite(sprintf('%s/%s/%s_disp_lentol_%.2f_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol, obj.outversion),...
                'structured_grid', obj.startPointsInROI(:, 1),...
                obj.startPointsInROI(:, 2), obj.startPointsInROI(:, 3), ...
                'vectors', 'vector_field', obj.displacementInROI(:, 1),...
                obj.displacementInROI(:, 2), obj.displacementInROI(:, 3));
            vtkwrite(sprintf('%s/%s/%s_directionOnStart_lentol_%.2f_jp=%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,obj.analysisConstants.pathlineJump,obj.outversion),...
                'structured_grid', obj.startPointsInROI(1:obj.analysisConstants.pathlineJump:end, 1),...
                obj.startPointsInROI(1:obj.analysisConstants.pathlineJump:end, 2), obj.startPointsInROI(1:obj.analysisConstants.pathlineJump:end, 3), ...
                'vectors', 'vector_field', obj.directionInROI(1:obj.analysisConstants.pathlineJump:end, 1),...
                obj.directionInROI(1:obj.analysisConstants.pathlineJump:end, 2), obj.directionInROI(1:obj.analysisConstants.pathlineJump:end, 3));
            if ~strcmp(configuration.diffusionCoefficientType,'constant')
                vtkwrite_spdlines(sprintf('%s/%s/%s_dCoefflines_lentol_%.1f_jp_%d_%s.vtk',...
                    obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                    obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                    obj.positionInROI(1:obj.analysisConstants.pathlineJump:end),...
                    obj.diffusionCoeffcientInROI(1:obj.analysisConstants.pathlineJump:end), 'precision', 6);                
            end
            if ~isempty(obj.anatomyImage)
                vtkwrite(sprintf('%s/%s/%s_anato_%s.vtk', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion), 'structured_points',...
                    'mask', obj.anatomyImage);
            else
                vtkwrite(sprintf('%s/%s/%s_mask_%s.vtk', obj.pathOutput,...
                    obj.outdir, obj.tag, obj.outversion), 'structured_points',...
                    'mask', obj.experimentMask.contents);
            end
            if obj.addSource
                vtkwrite_spdlines(sprintf('%s/%s/%s_Sourcelines_lentol_%.1f_jp_%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                obj.positionInROI(1:obj.analysisConstants.pathlineJump:end),...
                obj.sourceInROI(1:obj.analysisConstants.pathlineJump:end));
            end
            %{
            for i=1:length(obj.positionInROI)
                obj.positionInROI{i} = obj.positionInROI{i} + [configuration.xRange(1)-1, ... 
                    configuration.yRange(1)-1,configuration.zRange(1)-1];
                obj.startPointsInROI = obj.startPointsInROI + [configuration.xRange(1)-1, ... 
                    configuration.yRange(1)-1,configuration.zRange(1)-1];
            end
            vtkwrite_pathlines(sprintf('%s/%s/%s_pathlinesInOrigin_lentol_%.1f_jp_%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                obj.positionInROI(1:obj.analysisConstants.pathlineJump:end));
            
            vtkwrite_spdlines(sprintf('%s/%s/%s_PelinesInOrigin_lentol_%.1f_jp_%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                obj.positionInROI(1:obj.analysisConstants.pathlineJump:end),...
                obj.PeInROI(1:obj.analysisConstants.pathlineJump:end));
            
            vtkwrite_spdlines(sprintf('%s/%s/%s_SpdlinesInOrigin_lentol_%.1f_jp_%d_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol,...
                obj.analysisConstants.pathlineJump, obj.outversion), 'polydata', 'lines',...
                obj.positionInROI(1:obj.analysisConstants.pathlineJump:end),...
                obj.speedInROI(1:obj.analysisConstants.pathlineJump:end));
            vtkwrite(sprintf('%s/%s/%s_dispInOrigin_lentol_%.2f_%s.vtk',...
                obj.pathOutput, obj.outdir, obj.tag, obj.streamlineLengthTol, obj.outversion),...
                'structured_grid', obj.startPointsInROI(:, 1),...
                obj.startPointsInROI(:, 2), obj.startPointsInROI(:, 3), ...
                'vectors', 'vector_field', obj.displacementInROI(:, 1),...
                obj.displacementInROI(:, 2), obj.displacementInROI(:, 3));
            
            %}
            fprintf('Pathlines and Flux vectors in vtk format saved in %s/%s\n\n',...
                obj.pathOutput, obj.outdir)
        end
        function visualize(obj)
            % visualize pathlines in matlab
            obj.visualizePathlines;
            % visualize speed-lines in matlab
            obj.visualizeSpeedlines;
            % visualize peclet-lines in matlab
            obj.visualizePelines;
            % visualize velocity flux vectors in matlab
            obj.visualizeDisplacementVectors;
        end
        function visualizePathlines(obj)
            configuration = Configuration3D(obj.jsonFileName);
            %% pathlines
            figure,
            SL2 = obj.positionInROI;
            nSL = length(SL2);
            %colors = jet(round(max(PATH.displen)));
            if isempty(configuration.mask) 
                configuration.mask = Mask(configuration.pathMask,configuration.isMaskFilled,configuration.xRange,configuration.yRange,configuration.zRange);
                if configuration.do_resize
                    configuration.mask = configuration.mask.resize(configuration.sizeFactor);
                end

                if configuration.dilate>0
                    configuration.mask = configuration.mask.dilate(configuration.dilate,configuration.dilateDim);
                end
            end
            for ind = 1:10:nSL
                SL_tmp = SL2{ind};
                colors = jet(size(SL_tmp,1));
                hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
                %set(hlines,'EdgeColor',colors(round(PATH.displen(ind)),:));
                set(hlines,'FaceVertexCData',[colors;colors(end,:)],'EdgeColor','flat','FaceColor','none');
            end
            %
            hold on;
            [x, y, z] = meshgrid(1:configuration.trueSize(2), 1:configuration.trueSize(1), 1:configuration.trueSize(3));
            mskfv = isosurface(x,y,z,configuration.mask.contents,0.5);
            mskp = patch(mskfv);
            mskp.FaceColor = [.37,.37,.37];
            mskp.FaceAlpha= 0.1;
            mskp.EdgeColor = [.37,.37,.37];
            mskp.EdgeAlpha= 0;
            axis image
            colormap('jet');
            view([242.1011   14.4475])
            xticks(0:5:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
            ax = gca; ax.FontSize = 10; 
            xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
            set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
            set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
            set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
            grid on,
            xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
            
            cb = colorbar;
            set(cb,'YTick',[0,1])
            cb.TickLabels = {'start','end'};
            % text(1,-3,53,'end','Rotation',90,'FontSize',25);
            % text(1,-3,-1,'start','Rotation',90,'FontSize',25);
            title('pathlines')
            %saveas(gcf, sprintf('%s/%s/%s_LagPathlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
        end

        function visualizeSpeedlines(obj)
            configuration = Configuration3D(obj.jsonFileName);
            %% Speed-lines
            spdSet = cell2mat(obj.speedInROI');
            maxset = maxk(spdSet(:),round(0.05*length(spdSet(:))));
            spdmax = maxset(end);
            figure,
            Num = 100; % the larger, the more intervals in colors
            
            SL2 = obj.positionInROI;
            SL_spd = obj.speedInROI;
            nSL = length(SL2);
            colors = jet(Num);
            if isempty(configuration.mask) 
                configuration.mask = Mask(configuration.pathMask,configuration.isMaskFilled,configuration.xRange,configuration.yRange,configuration.zRange);
                if configuration.do_resize
                    configuration.mask = configuration.mask.resize(configuration.sizeFactor);
                end

                if configuration.dilate>0
                    configuration.mask = configuration.mask.dilate(configuration.dilate,configuration.dilateDim);
                end
            end
            for ind = 1:10:nSL
                SL_tmp = SL2{ind};
                SL_spd_tmp = SL_spd{ind};
                SL_spd_tmp(SL_spd_tmp>spdmax) = spdmax;
                SL_spd_rk = round(SL_spd_tmp/spdmax*Num);
                SL_spd_rk(SL_spd_rk<1) = 1;
                hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
                set(hlines,'FaceVertexCData',[colors(SL_spd_rk,:);colors(SL_spd_rk(end),:)],'EdgeColor','flat','FaceColor','none');
            end
            %
            hold on;
            [x, y, z] = meshgrid(1:configuration.trueSize(2), 1:configuration.trueSize(1), 1:configuration.trueSize(3));
            mskfv = isosurface(x,y,z,configuration.mask.contents,0.5);
            mskp = patch(mskfv);
            mskp.FaceColor = [.37,.37,.37];
            mskp.FaceAlpha= 0.1;
            mskp.EdgeColor = [.37,.37,.37];
            mskp.EdgeAlpha= 0;
            view([242.1011   14.4475])
            xticks(0:5:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
            ax = gca; ax.FontSize = 10; 
            xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
            set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
            set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
            set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
            grid on; axis image;
            xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
            colormap('jet'); grid on;
            
            cb = colorbar;
            cb.Ticks = linspace(0, 1, 6);
            cb.TickLabels = num2cell(round((0:spdmax/5:spdmax)*100)/100);
            cb.FontSize = 8;
            cb.Label.String = 'speed (a.u.)';
            %text(1,-3,15,'speed (a.u.)','Rotation',90,'FontSize',25);
            title('speed-lines')
            xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
            %saveas(gcf, sprintf('%s/%s/%s_LagSpdlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
        end

        function visualizePelines(obj)
            configuration = Configuration3D(obj.jsonFileName);
            %% pelines
            peSet = cell2mat(obj.PeInROI');
            maxset = maxk(peSet(:),round(0.05*length(peSet(:))));
            Pemax = maxset(end);
            figure,
            %Pemax = 1000;
            Num = 100; % the larger, the more intervals in colors
            
            SL2 = obj.positionInROI;
            SL_Pe = obj.PeInROI;
            nSL = length(SL2);
            colors = jet(Num);
            if isempty(configuration.mask) 
                configuration.mask = Mask(configuration.pathMask,configuration.isMaskFilled,configuration.xRange,configuration.yRange,configuration.zRange);
                if configuration.do_resize
                    configuration.mask = configuration.mask.resize(configuration.sizeFactor);
                end

                if configuration.dilate>0
                    configuration.mask = configuration.mask.dilate(configuration.dilate,configuration.dilateDim);
                end
            end
            for ind = 1:10:nSL
                SL_tmp = SL2{ind};
                SL_Pe_tmp = SL_Pe{ind};
                SL_Pe_tmp(SL_Pe_tmp>Pemax) = Pemax;
                SL_Pe_rk = round(SL_Pe_tmp/Pemax*Num);
                SL_Pe_rk(SL_Pe_rk<1) = 1;
                hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
                set(hlines,'FaceVertexCData',[colors(SL_Pe_rk,:);colors(SL_Pe_rk(end),:)],'EdgeColor','flat','FaceColor','none');
            end
            %
            hold on;
            [x, y, z] = meshgrid(1:configuration.trueSize(2), 1:configuration.trueSize(1), 1:configuration.trueSize(3));
            mskfv = isosurface(x,y,z,configuration.mask.contents,0.5);
            mskp = patch(mskfv);
            mskp.FaceColor = [.37,.37,.37];
            mskp.FaceAlpha= 0.1;
            mskp.EdgeColor = [.37,.37,.37];
            mskp.EdgeAlpha= 0;
            view([242.1011   14.4475])
            xticks(0:5:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
            ax = gca; ax.FontSize = 10; 
            xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
            set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
            set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
            set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
            grid on; axis image;
            xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
            colormap('jet'); grid on;
            
            cb = colorbar;
            cb.Ticks = linspace(0, 1, 6);
            cb.TickLabels = num2cell(round((0:Pemax/5:Pemax)));
            cb.FontSize = 8;
            cb.Label.String = '\it{Pe}';
            %text(1,-3,25,'\it{Pe}','Rotation',90,'FontSize',25);
            title('pelines')
            xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
            %saveas(gcf, sprintf('%s/%s/%s_LagPelines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
        end

        function visualizeDisplacementVectors(obj)
            configuration = Configuration3D(obj.jsonFileName);
            if isempty(configuration.mask) 
                configuration.mask = Mask(configuration.pathMask,configuration.isMaskFilled,configuration.xRange,configuration.yRange,configuration.zRange);
                if configuration.do_resize
                    configuration.mask = configuration.mask.resize(configuration.sizeFactor);
                end

                if configuration.dilate>0
                    configuration.mask = configuration.mask.dilate(configuration.dilate,configuration.dilateDim);
                end
            end
            %% velocity flux vector 
            figure,
            strid = 10;
            magnify = 1;
            [x, y, z] = meshgrid(1:configuration.trueSize(2), 1:configuration.trueSize(1), 1:configuration.trueSize(3));
            mskfv = isosurface(x,y,z,configuration.mask.contents,0.5);
            mskp = patch(mskfv);
            mskp.FaceColor = [.17,.17,.17];
            mskp.FaceAlpha= 0.031;
            mskp.EdgeColor = [.17,.17,.17];
            mskp.EdgeAlpha= 0;
            mskp.DisplayName = 'mask';
            
            grid on, axis image
            hold on,
            q = quiver3(obj.startPointsInROI(1:strid:end,2),obj.startPointsInROI(1:strid:end,1),obj.startPointsInROI(1:strid:end,3),obj.displacementInROI(1:strid:end,2)*magnify,obj.displacementInROI(1:strid:end,1)*magnify,obj.displacementInROI(1:strid:end,3)*magnify,...
                'color','r','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','off','DisplayName','flux vectors');
            %title(sprintf('Velocity Flux Vectors'),'FontSize',20, 'Interpreter', 'none'), 
            view([242.1011   14.4475])
            ax = gca; ax.FontSize = 10; 
            xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
            set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
            set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
            set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.4808 0.3259 0.2937 0.3625],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
            grid on; axis image;
            xticks(0:5:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
            %// Compute the magnitude of the vectors
            mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
                        reshape(q.WData, numel(q.UData), [])).^2, 2));
            
            %// Get the current colormap
            currentColormap = colormap(jet);
            [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
            cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
            cmap(:,:,4) = 255;
            cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
            set(q.Head, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:3,:,:), [], 4).');   
            set(q.Tail, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:2,:,:), [], 4).');
            
            view([242.1011   14.4475])
            set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
            %colorbar('FontSize',18), 
            grid on,
            xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
            
            cb = colorbar;
            cb.Ticks = linspace(0, 800, 6);
            cb.TickLabels = num2cell(0:160:800);
            cb.FontSize = 8;
            cb.Label.String = 'distance of movement (a.u.)';
            title('velocity flux vector')
            
            %saveas(gcf, sprintf('%s/%s/%s_LagFluxVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
        end


    end
    
 methods (Access = protected)
        
        function initialize(obj,nstartpoints)
            %initialize pathlines
            obj.position = cell(1, nstartpoints);
            obj.rho = cell(1, nstartpoints);
            obj.drho_dx = cell(1, nstartpoints);
            obj.augmentedSpeed = cell(1, nstartpoints);
            obj.speed = cell(1, nstartpoints);
            obj.diffusiveSpeed = cell(1, nstartpoints);
            obj.timeStep = cell(1, nstartpoints);
            obj.flux = cell(1, nstartpoints);
            obj.drho_dt = cell(1, nstartpoints);
            obj.Pe = cell(1, nstartpoints);
            obj.diffusionCoefficient = cell(1, nstartpoints);
            if obj.addSource
                obj.source = cell(1,nstartpoints);
                obj.sFlux = cell(1,nstartpoints);
            end
        end
        
        function iPathlines = update(obj,lagrangianSpecification,npoints,nPathlines,iPathlines)
            %update pathlines
            %npoints: number of starting points
            %nPathlines: threshold of unique points for pathlines
            %iPathlines: count number of pathlines with unique points>pln
            for iStartPoints=1:npoints
                pl_cur = squeeze(lagrangianSpecification.position(iStartPoints, :, :))';
                aind = any(~isnan(pl_cur), 2); %1 if row is not NaN, 0 if row is NaN
                pl_cur = pl_cur(aind, :);
                %check that unique Pathlines has more than 1 point (remove not-move coordinates)
                [pl_cur, ia, ~] = unique(pl_cur, 'rows', 'stable');

                if size(pl_cur, 1) > nPathlines
                    iPathlines = iPathlines + 1;
                    obj.position{iPathlines} = pl_cur;

                    plr_cur = lagrangianSpecification.rho(iStartPoints, aind)';
                    pldrdx_cur = lagrangianSpecification.drho_dx(iStartPoints, aind)';
                    pls_cur = lagrangianSpecification.speed(iStartPoints, aind)';
                    plds_cur = lagrangianSpecification.diffusiveSpeed(iStartPoints, aind)';
                    plas_cur = lagrangianSpecification.augmentedSpeed(iStartPoints, aind)';
                    plt_cur = lagrangianSpecification.timeStep(iStartPoints, aind)';
                    plflx_cur = lagrangianSpecification.flux(iStartPoints, aind)';
                    pldr_cur = lagrangianSpecification.drho_dt(iStartPoints, aind)';
                    pldc_cur = lagrangianSpecification.diffusionCoefficient(iStartPoints, aind)';
                    if obj.addSource
                        plsr_cur = lagrangianSpecification.source(iStartPoints, aind)';
                        obj.source{iPathlines} = plsr_cur(ia);
                        plsrFlux_cur = lagrangianSpecification.sFlux(iStartPoints, aind)';
                        obj.sFlux{iPathlines} = plsrFlux_cur(ia);
                    end

                    obj.rho{iPathlines} = plr_cur(ia);
                    obj.drho_dx{iPathlines} = pldrdx_cur(ia);
                    obj.speed{iPathlines} = pls_cur(ia);
                    obj.diffusiveSpeed{iPathlines} = plds_cur(ia);
                    obj.augmentedSpeed{iPathlines} = plas_cur(ia);
                    obj.timeStep{iPathlines} = plt_cur(ia);
                    obj.flux{iPathlines} = plflx_cur(ia);
                    obj.drho_dt{iPathlines} = pldr_cur(ia);
                    obj.Pe{iPathlines} = pls_cur(ia) ./ (plds_cur(ia) + eps);
                    obj.diffusionCoefficient{iPathlines} = pldc_cur(ia);
                end
            end
        end
        
        function remove0(obj,iPathlines)
            % REMOVE0 remove extra created cells
            obj.position = obj.position(1:iPathlines);
            obj.rho = obj.rho(1:iPathlines);
            obj.drho_dx = obj.drho_dx(1:iPathlines);
            obj.speed = obj.speed(1:iPathlines);
            obj.diffusiveSpeed = obj.diffusiveSpeed(1:iPathlines);
            obj.augmentedSpeed = obj.augmentedSpeed(1:iPathlines);
            obj.timeStep = obj.timeStep(1:iPathlines);
            obj.flux = obj.flux(1:iPathlines);
            obj.drho_dt = obj.drho_dt(1:iPathlines);
            obj.Pe = obj.Pe(1:iPathlines);
            obj.diffusionCoefficient = obj.diffusionCoefficient(1:iPathlines);
            if obj.addSource
                obj.source = obj.source(1:iPathlines);
                obj.sFlux = obj.sFlux(1:iPathlines);
            end
        end
        
        function doThresholding(obj,tolerance)
            % DOTHRESHOLDING remove pathlines with euclidean length less
            % than tolerance
            
            %get cell array with euclidean length of sl between first and last point
            eucdist = cellfun(@(x) sqrt(sum((x(end, :) - x(1, :)).^2)), obj.position); 
            % only keep streamlines whose Euclidean length between first 
            % and last points is larger than the threshold
            obj.position = obj.position(eucdist > tolerance); 
            obj.rho = obj.rho(eucdist > tolerance);
            obj.drho_dx = obj.drho_dx(eucdist>tolerance);
            obj.speed = obj.speed(eucdist > tolerance);
            obj.diffusiveSpeed = obj.diffusiveSpeed(eucdist > tolerance);
            obj.augmentedSpeed = obj.augmentedSpeed(eucdist > tolerance);
            obj.timeStep = obj.timeStep(eucdist > tolerance);
            obj.flux = obj.flux(eucdist > tolerance);
            obj.drho_dt = obj.drho_dt(eucdist > tolerance);
            obj.Pe = obj.Pe(eucdist > tolerance);
            obj.diffusionCoefficient = obj.diffusionCoefficient(eucdist > tolerance);
            if obj.addSource
                obj.source = obj.source(eucdist > tolerance);
                obj.sFlux = obj.sFlux(eucdist > tolerance);
            end
        end
        
    end
    
        
    methods (Static)
        
        function deleteOldInfo
            % make sure don't use pathlines information from previous timestep
            File1 = fullfile(cd, 'pl_cur.mat');
            File2 = fullfile(cd, 'pli_array.mat');
            File3 = fullfile(cd, 'pl_centroid_array.mat');

            if exist(File1, 'file')
                delete(File1);
            end

            if exist(File2, 'file')
                delete(File2);
            end

            if exist(File3, 'file')
                delete(File3);
            end
        end
        
    end

end

