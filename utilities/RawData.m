classdef RawData < Configuration
    %RAWDATA Data without any preprocessing and functions that analyze raw
    %data
    %   Detailed explanation goes here
    
    properties
        height
        width
        depth
        contents % intensities of imaging data.
        rawMask
        niftyForm % used for saving as nifty file
    end
    
    methods
        function obj = RawData(tag)
            %RAWDATA Construct an instance of this class
            %   Construct an instance of RawDara based on Configuration and
            %   load raw data without doing any processing.
            
            addpath('./Sensitivities', './analyzeFlows', genpath('./utilities'))

            obj = obj@Configuration(tag);
            obj.pathOutput = sprintf('./output/%s/RawData',obj.tag);
            obj.nFiles = obj.data_end - obj.data_start + 1;
            obj.load;
            if ~isempty(obj.pathMask)
                X = load_untouch_nii(obj.pathMask);
                obj.rawMask = double(X.img);
            else
                obj.rawMask = ones(obj.height,obj.width,obj.depth);
            end         
        end
        
        function load(obj)
            %LOAD load the data without doing any processing            
            if strcmp(obj.extension,'.gz') % packed 
                X = load_untouch_nii(obj.pathData);
                obj.niftyForm = X;
                obj.niftyForm.hdr.dime.dim(1) = 3;
                obj.niftyForm.hdr.dime.dim(5) = 1;
                obj.niftyForm.hdr.dime.pixdim(5) = 1;
                obj.contents = double(X.img);
                [obj.height,obj.width,obj.depth,~] = size(obj.contents);               
            elseif strcmp(obj.extension,'.nii') % unpacked
                X = load_untouch_nii(sprintf('%s%02d%s', obj.pathData, ...
                        obj.data_start, '.nii'));
                [obj.height,obj.width,obj.depth] = size(X.img);
                obj.niftyForm = X;
                obj.contents = zeros(obj.height,obj.width,obj.depth,obj.nFiles);
                obj.contents(:,:,:,1) = double(X.img);
                for kk = 2:obj.nFiles
                    X = load_untouch_nii(sprintf('%s%02d%s', obj.pathData, ...
                        obj.data_start + kk - 1, '.nii'));   
                    obj.contents(:,:,:,kk) = double(X.img);
                end
            else
                error("Data must be packed ('*nii.gz') or unpacked ('*.nii') NIfTI files")
            end
        end
        
        function [peakValue, timeToPeak] = peak(obj,mask)
            %PEAK Calculate maximum value over time for each voxel in ROI
            %and the time that maximum value occurs.
            if nargin < 2
                mask = obj.rawMask;
            end
            peakValue = zeros(obj.height,obj.width,obj.depth);
            timeToPeak = zeros(obj.height,obj.width,obj.depth);
            ind = find(mask == 1);
            [x, y, z] = ind2sub([obj.height,obj.width,obj.depth], ind);
            xmin = min(x);
            xmax = max(x);
            ymin = min(y);
            ymax = max(y);
            zmin = min(z);
            zmax = max(z);
            [peakValue(xmin:xmax,ymin:ymax,zmin:zmax), timeToPeak(xmin:xmax,ymin:ymax,zmin:zmax)] = ...
                max(obj.contents(xmin:xmax,ymin:ymax,zmin:zmax,:),[],4);
            peakValue(xmin:xmax,ymin:ymax,zmin:zmax) = peakValue(xmin:xmax,ymin:ymax,zmin:zmax).* ...
                                                        mask(xmin:xmax,ymin:ymax,zmin:zmax);
            timeToPeak(xmin:xmax,ymin:ymax,zmin:zmax) = timeToPeak(xmin:xmax,ymin:ymax,zmin:zmax).* ...
                                                        mask(xmin:xmax,ymin:ymax,zmin:zmax);
        end
        
        function timeValue = timeToHalfPeak(obj,mask)
            if nargin < 2
                mask = obj.rawMask;
            end
            [peakValue, timeToPeak] = obj.peak(mask);
            timeValue = zeros(obj.height,obj.width,obj.depth);
            
            % redundant calculation, may be saved in class mask later.
            ind = find(mask == 1);
            [x, y, z] = ind2sub([obj.height,obj.width,obj.depth], ind);
            xmin = min(x);
            xmax = max(x);
            ymin = min(y);
            ymax = max(y);
            zmin = min(z);
            zmax = max(z);
            
            for i = xmin:xmax
                for j = ymin:ymax
                    for k = zmin:zmax
                        if mask(i,j,k)
                            [~,timeValue(i,j,k)] = min(abs(obj.contents(i,j,k,1:timeToPeak(i,j,k))-0.5*peakValue(i,j,k)));
                        end
                    end
                end
            end
        end
        
        function exportTimeToHalfPeak(obj)
            timeValue = obj.timeToHalfPeak;
            if ~exist(obj.pathOutput, 'dir')
                mkdir(obj.pathOutput)
            end
            niftyOutput = obj.niftyForm;
            niftyOutput.img = timeValue;
            save_untouch_nii(niftyOutput, sprintf('%s/timeToHalfPeak.nii', obj.pathOutput));
        end
    end
end

