classdef Mask
    %MASK Subsets of points used in computing. 
    
    properties
        contents %the 3D array of binary coordinates
        niftiForm
        path
        xRange
        yRange
        zRange
        threshold = 0
        resizeFactor = 0
        dilateFactor = 0
        erodeFactor = 0
    end
    
    methods
        
        function obj = Mask(pathMask,isMaskFilled,xRange,yRange,zRange,threshold)
            if nargin == 6
                obj.threshold = threshold;
            end
            if nargin < 5 
                [xRange, yRange, zRange] = find_range(pathMask);
            end
            if isempty(isMaskFilled)
                isMaskFilled = 1;
            end
            obj.path = pathMask;
            obj.xRange = xRange;
            obj.yRange = yRange;
            obj.zRange = zRange;
            if isempty(obj.path)
                obj.path = '';
            end
            if isfile(obj.path)
                A = nii2mat(obj.path, obj.xRange, obj.yRange, obj.zRange);
                obj.contents = zeros(size(A));
                obj.contents(A > obj.threshold) = 1;
            else
                obj.contents = ones(length(obj.xRange), length(obj.yRange), length(obj.zRange));
            end
            if ~isMaskFilled
                obj = obj.fill;
            end
        end
        
        function dilatedMask = dilate(obj,d)
            %DILATE dilates a mask by a factor d
            [xr, yr, zr] = meshgrid(-d:d, -d:d, -d:d);
            strel = (xr / d).^2 + (yr / d).^2 + (zr / d).^2 <= 1;
            dilatedMask = obj;
            dilatedMask.contents = imdilate(obj.contents, strel);
            dilatedMask.dilateFactor = d;
        end
        
        function erodedMask = erode(obj,e)
            [xer, yer, zer] = meshgrid(-e:e, -e:e, -1:1);
            strel_e = (xer / d).^2 + (yer / obj.d).^2 + (zer).^2 <= 1;
            erodedMask = obj;
            erodedMask.contents = imerode(obj.contents, strel_e);
            erodedMask.erodeFactor = e;
        end
        
        function resizedMask = resize(obj,s)
            %RESIZE resize a mask by a factor s
            resizedMask = obj;
            resizedMask.contents = resizeMatrix(obj.contents, round(s .* size(obj.contents)), 'linear');
            resizedMask.contents(resizedMask.contents ~= 1) = 0;
            resizedMask.resizeFactor = s;
        end
        
        function filledMask = fill(obj)
            %FILL fill the hollow mask
            filledMask = obj;            
            for i = 1:size(obj.contents,3)
                filledMask.contents(:,:,i) = imfill(obj.contents(:,:,i));
            end
        end
        
    end
end

