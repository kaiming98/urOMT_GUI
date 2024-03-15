classdef FlowFieldCharacteristic < handle
    % A superclass for calculating flow field characterstics from the
    % results of the optimized velocity field found using the
    % regularized optimal mass transport (rOMT) algorithm. 
    %
    % Inputs for instance construction: 
    %   configuration  - instance of the Configuration class where the data is
    %                    from


    properties (Access=public)
        tag %dataset identifier (e.g., 'C294')
        jsonFileName
        analysisConstants
        experimentMask Mask
        cutoffThresholds
    end

    properties (Access=protected)
        streamlineLengthTol double = 2;%2;glymphatic 
        %output strings
        sigmaStr
        paperFigStr
        dateStr
        outdir
        outdirLong
        outversion
        %configuration constants
        pathMask
        pathOutput
        trueSize
        nt
        timeInitial
        timeFinal
        timeJump
    end
    
    methods

        function obj = FlowFieldCharacteristic(configuration) 
            % transfer important constants from configuration to the instance 
            obj.tag = configuration.tag;
            obj.pathMask = configuration.pathMask;
            obj.pathOutput = configuration.pathOutput;   
            obj.nt = configuration.nt;
            obj.timeInitial = configuration.timeInitial;
            obj.timeFinal = configuration.timeFinal;
            obj.timeJump = configuration.timeJump;
            obj.trueSize = configuration.trueSize;           
        end
               
        function initializeOutDir(obj,analysisNum,configuration,startPoints)
            %INITIALIZEOUTDIR set output directory related strings
            obj.sigmaStr = num2str(configuration.sigma);
            obj.paperFigStr = sprintf('set0%02d', analysisNum);
            formatOut = 'mmddyy';
            obj.dateStr = datestr(now, formatOut);
            obj.outversion = sprintf('%s_%s', obj.paperFigStr, obj.dateStr);
            obj.outdir = sprintf('LPPA_%s', obj.outversion);
            switch analysisNum
                case 1
                    obj.outdir = 'lagrangian/speed';
                case 2
                    obj.outdir = 'lagrangian/pathlines';
                case 3
                    obj.outdir = 'eulerian';
            end
            
            switch analysisNum
                case 1
                    obj.outdirLong = sprintf(['LPPA_%s_mask%d_type_%s_%s_'...
                        'img%d_mdt%d%s_%s_erode%d_dilate%d_nEul%d_cutoffStr_'...
                        '%s_concThresh%5.4f_spdThresh%5.4f_imBoundFlag_%d_%s_'...
                        'slTol%d_clusters_centroid%d_nbp%d_metric%s_'...
                        'qbthresh%d_clusTol%d_clusCutoff%d_diff%s_tj%d_'...
                        'nt%d%s_%s_%s'], obj.tag, configuration.mask_number,...
                        'speedmap', obj.analysisConstants.flowType, obj.analysisConstants.isInterpolated,...
                        obj.analysisConstants.nTimeInterval, obj.analysisConstants.distanceFlag,...
                        startPoints.maskString, startPoints.erodeMask, startPoints.dilateMask,...
                        obj.analysisConstants.nEulerianSteps, obj.analysisConstants.cutoffFlag,...
                        obj.cutoffThresholds.density, obj.cutoffThresholds.speed, obj.analysisConstants.intensityRangeFlag,...
                        startPoints.method, obj.streamlineLengthTol, obj.analysisConstants.centroidFlag,...
                        obj.analysisConstants.nBundles, obj.metricNameShort,...
                        obj.analysisConstants.quickBundleThreshold, obj.analysisConstants.clusterTol,...
                        obj.analysisConstants.clusterCutoff, obj.sigmaStr, obj.timeJump,...
                        obj.nt, startPoints.thresholdString, obj.paperFigStr, obj.dateStr);
                case 2
                    obj.outdirLong = sprintf(['LPPA_%s_mask%d_type_%s_%s_'...
                        'img%d_mdt%d%s_%s_erode%d_dilate%d_nEul%d_cutoffStr_'...
                        '%s_concThresh%5.4f_spdThresh%5.4f_imBoundFlag_%d_%s_slTol%d_'...
                        'diff%s_tj%d_nt%d%s_%s_%s'], obj.tag,...
                        configuration.mask_number, 'vectors', obj.analysisConstants.flowType,...
                        obj.analysisConstants.isInterpolated, obj.analysisConstants.nTimeInterval,...
                        obj.analysisConstants.distanceFlag,startPoints.maskString, startPoints.erodeMask,...
                        startPoints.dilateMask, obj.analysisConstants.nEulerianSteps,...
                        obj.analysisConstants.cutoffFlag, obj.cutoffThresholds.density,...
                        obj.cutoffThresholds.speed, obj.analysisConstants.intensityRangeFlag,...
                        startPoints.method, obj.streamlineLengthTol, obj.sigmaStr,...
                        obj.timeJump, obj.nt, startPoints.thresholdString, obj.paperFigStr, obj.dateStr);
            end
        end   
        
        function createOutDir(obj,out_dir,tag)
            out_dir = sprintf('%s/%s', out_dir, obj.outdir);
            if ~exist(out_dir, 'dir')
                mkdir(out_dir)
            end

            fid = fopen(sprintf('%s/%s_record_%s_%s.txt', out_dir, tag,...
                obj.paperFigStr, obj.dateStr), 'a+');
            fprintf(fid, '%s/%s_record_%s_%s\noutdir-long: %s', out_dir,...
                tag, obj.paperFigStr, obj.dateStr, obj.outdirLong);  
        end
        
        
    end

end

