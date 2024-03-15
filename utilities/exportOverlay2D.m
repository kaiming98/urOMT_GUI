function imOut = exportOverlay2D(im, overlaidProperty, tag, direction, level)
%EXPORTOVERLAID2D display and export images that overlay the properties on
%the data in 2D.
%   Inputs:
%       im - 2D data image
%       overlaidProperty - property we want used (3d).
%       tag - configuration tag
%       direction - the direction we fixed coordinates, for example, 'z' 
%       means we are looked at xy-plane.
%       level - slice level. Along the direction, show the level
%       coordinate.
%   colormap: jet(blue to red)
obj = Configuration3D(tag);
overlaidProperty = squeeze(overlaidProperty);
im = double(im);
[height,width,~] = size(im);
%load mask
if isempty(obj.mask) 
    if ~isempty(obj.pathMask)
        mask = Mask(obj.pathMask,obj.isMaskFilled,obj.xRange,obj.yRange,obj.zRange);
        if obj.do_resize
            mask = mask.resize(obj.sizeFactor);
        end

        if obj.dilate>0
            mask = mask.dilate(obj.dilate);
        end
    else
        error("missing path to mask")
    end
end
switch direction
    case 'z'
        propertyMap = zeros(obj.rawSize(2),obj.rawSize(1),3);
        largeMask = zeros(obj.rawSize(2),obj.rawSize(1));
        overlaidProperty = squeeze(overlaidProperty(:,:,level));
        A = min(overlaidProperty,[],'all');
        B = max(overlaidProperty,[],'all');
        for i = obj.yRange
            for j = obj.xRange
                if mask.contents(j-obj.xRange(1)+1,i-obj.yRange(1)+1,level)>0
                    largeMask(i,j) = 1;
                    propertyMap(i,j,1) = 255/(B-A)*(overlaidProperty(j-obj.xRange(1)+1,i-obj.yRange(1)+1)-A);
                    propertyMap(i,j,2) = 0;
                    propertyMap(i,j,3) = 255/(A-B)*(overlaidProperty(j-obj.xRange(1)+1,i-obj.yRange(1)+1)-B);
                end
            end
        end
    case 'x'
        propertyMap = zeros(obj.rawSize(3),obj.rawSize(2),3);
        largeMask = zeros(obj.rawSize(3),obj.rawSize(2));
        overlaidProperty = squeeze(overlaidProperty(level,:,:));
        A = min(overlaidProperty,[],'all');
        B = max(overlaidProperty,[],'all');
        for i = obj.zRange
            for j = obj.yRange
                if mask.contents(level,j-obj.yRange(1)+1,i-obj.zRange(1)+1)>0
                    largeMask(i,j) = 1;
                    propertyMap(i,j,1) = 255/(B-A)*(overlaidProperty(j-obj.yRange(1)+1,i-obj.zRange(1)+1)-A);
                    propertyMap(i,j,2) = 0;
                    propertyMap(i,j,3) = 255/(A-B)*(overlaidProperty(j-obj.yRange(1)+1,i-obj.zRange(1)+1)-B);
                end
            end
        end
    case 'y'
        propertyMap = zeros(obj.rawSize(3),obj.rawSize(1),3);
        largeMask = zeros(obj.rawSize(3),obj.rawSize(1));
        overlaidProperty = squeeze(overlaidProperty(:,level,:));
        A = min(overlaidProperty,[],'all');
        B = max(overlaidProperty,[],'all');
        for i = obj.zRange
            for j = obj.xRange
                if mask.contents(j-obj.xRange(1)+1,level,i-obj.zRange(1)+1)>0
                    largeMask(i,j) = 1;
                    propertyMap(i,j,1) = 255/(B-A)*(overlaidProperty(j-obj.xRange(1)+1,i-obj.zRange(1)+1)-A);
                    propertyMap(i,j,2) = 0;
                    propertyMap(i,j,3) = 255/(A-B)*(overlaidProperty(j-obj.xRange(1)+1,i-obj.zRange(1)+1)-B);
                end
            end
        end
end
largeMask = repmat(largeMask,[1,1,3]);
largeMask = imresize(largeMask,[height,width]);
propertyMap = imresize(propertyMap,[height,width]);
imOut = im.*(largeMask==0)+propertyMap.*(largeMask>0);
imOut = uint8(imOut);
imshow(imOut);
end

