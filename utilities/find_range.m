function [x, y, z, msk_size] = find_range(mask_path,dilate)
    %find the rectangle around the dilated mask
    %x = min(x(msk==1))-2 : max(x(msk==1))+2 
    %msk_size: mask size
    if nargin == 1
        dilate = 0;
    end
    mask = load_untouch_nii(mask_path);
    msk = mask.img;
    msk = bwmorph3(msk, 'fill');
    if dilate > 0
        [xr, yr, zr] = meshgrid(-dilate:dilate, -dilate:dilate, -dilate:dilate);
        strel = (xr / dilate).^2 + (yr / dilate).^2 + (zr / dilate).^2 <= 1;
        msk = imdilate(msk, strel);
    end
    ind = find(msk == 1);
    msk_size = size(msk);
    [x, y, z] = ind2sub(msk_size, ind);
    x1 = min(x);
    x2 = max(x);
    y1 = min(y);
    y2 = max(y);
    z1 = min(z);
    z2 = max(z);
    x1 = max([1, x1 - 2]);
    y1 = max([1, y1 - 2]);
    z1 = max([1, z1 - 2]);
    x2 = min([msk_size(1), x2 + 2]);
    y2 = min([msk_size(2), y2 + 2]);
    z2 = min([msk_size(3), z2 + 2]);
    x = x1:x2;
    y = y1:y2;
    z = z1:z2;
end
