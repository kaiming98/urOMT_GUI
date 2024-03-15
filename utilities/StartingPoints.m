classdef StartingPoints < handle
    %STARTINGPOINTS Subset of points used in calculating lagrangian
    %Specification of flow field
    %   Detailed explanation goes here

    properties
        x %x coordinates of StartingPoints
        y %y coordinates of StartingPoints
        z %z coordinates of StartingPoints
        method
        methodTitle
        mask
        pathMask
        maskFlag = 0
        maskName = "olf_nomain"
        erodeMask = 0
        dilateMask = 0
        thresholdFlag = 1
        threshold = 12
        thresholdString = ''
        maskString = 'altSPmsk0'
    end
    
    properties (Access = protected)
        indexThresholded
    end

    methods
        function obj = StartingPoints(pathMask,maskName)
            arguments
                pathMask %{mustBeText}
                maskName %{mustBeText}
            end
            obj.pathMask = pathMask;
            obj.maskName = maskName;
            if isempty(obj.pathMask)
                obj.maskFlag = 0;
            end

            if obj.maskFlag
                obj.maskString = sprintf('altSPmsk_%s', obj.maskName);
            end

            if obj.erodeMask > 0 && obj.dilateMask > 0
                warning(['starting points mask cannot be both eroded and '...
                    'dilated. starting points mask will be eroded (no dilation) by default']);
                obj.dilateMask = 0;
            end             

        end

        function loadMask(obj,experimentMask,experiment)
            
            if obj.maskFlag
                obj.mask = Mask(obj.pathMask, experiment.isMaskFilled,...
                    experiment.xRange, experiment.yRange, experiment.zRange);
                if experiment.do_resize
                    obj.mask = obj.mask.resize(experiment.sizeFactor);
                end
                obj.mask.matrixForm(~experimentMask.matrixForm) = 0;
            else
                obj.mask = experimentMask;
            end
            
            if obj.erodeMask > 0
                obj.mask = obj.mask.erode(obj.erodeMask);
            elseif obj.dilateMask > 0
                %first find if original mask was dilated for rOMT run
                dilateFactor = experiment.dilate;
                if dilateFactor > 0
                    maskDIL = experimentMask.dilate(dilateFactor);
                else
                    maskDIL = experimentMask;
                end
                %now dilate spmask and remove any points outside of mask used for rOMT
                obj.mask = obj.mask.dilate(obj.dilateMask);
                obj.mask.contents(~maskDIL.contents) = 0;
            end
        end
        
        function doThresholding(obj,tag)
            % starting points thresholding:Only using the starts points 
            % from max(dpsnrv)>startpointsThresh
            configuration = Configuration3D(tag);
            if obj.thresholdFlag
                if ~isempty(configuration.pathMaxDensity)
                    maxDensity = nii2mat(configuration.pathMaxDensity, ...
                        configuration.xRange,configuration.yRange,configuration.zRange);
                else
                    configuration = Configuration3D(tag);
                    configuration.loadVolume;

                    maxDensityArray = zeros(1,length(configuration.Data));
                    for j = 1:length(configuration.Data)
                        maxDensityArray(j) = max(max(max(configuration.Data(j).volume)));
                    end
                    [~,iMaxDensity] = max(maxDensityArray);
                    maxDensity = configuration.Data(iMaxDensity).volume;
                end

                obj.indexThresholded = find((obj.mask.contents > 0) & (maxDensity > obj.threshold));
                obj.thresholdString = sprintf('_dpsnrv_min_%d', obj.threshold);
            else
                obj.indexThresholded = find(obj.mask.contents > 0);
            end
                      
        end
        
    end
end