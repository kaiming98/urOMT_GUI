function [sliceYZ,sliceXZ,sliceXY] = findLargestSlice(data)
    %FINDLARGESTSLICE find the regular slice that nonzero region is enclosed by
    %the largest rectangle.
    
    % find sliceYZ, the x coordinate of the largest yz-slice(parallel to yz-plane)
    area = 0;
    for x=1:size(data,1)
        slice = squeeze(data(x,:,:));
        ind = find(slice>0);
        [y,z] = ind2sub(size(slice),ind);
        if isempty(y)
            areaTmp=0;
        else
            areaTmp = (max(y)-min(y))*(max(z)-min(z));
        end
        if areaTmp>area
            area = areaTmp;
            sliceYZ = x;
        end
    end
    
    %find sliceXZ
    area = 0;
    for y=1:size(data,2)
        slice = squeeze(data(:,y,:));
        ind = find(slice>0);
        [x,z] = ind2sub(size(slice),ind);
        if isempty(x)
            areaTmp=0;
        else
            areaTmp = (max(x)-min(x))*(max(z)-min(z));
        end
        if areaTmp>area
            area = areaTmp;
            sliceXZ = y;
        end
    end
    
    % find sliceXY
    area = 0;
    for z=1:size(data,3)
        slice = squeeze(data(:,:,z));
        ind = find(slice>0);
        [x,y] = ind2sub(size(slice),ind);
        if isempty(x)
            areaTmp=0;
        else
            areaTmp = (max(x)-min(x))*(max(y)-min(y));
        end
        if areaTmp>area
            area = areaTmp;
            sliceXY = z;
        end
    end
end

