function im_out = nii2mat(data_path,x1_range,x2_range,x3_range,size_factor,up_thresh)
% this function extracts the img from the NIfTI file and returns the desired
% range and desired rescaled image as a double with no processing or 
% interpolation
%
% Input:
% >>> x1_range     = range of rows
% >>> x2_range     = range of cols
% >>> x3_range     = range of depth/layers
% >>> data_path    = path to .nii file
% >>> size_factor  = factor by which img should be resized. default, set to 1 if
%                    should not resize
% Output:
% <<< im_out       = image matrix

arguments
    data_path string
    x1_range = 1:256
    x2_range = 1:256
    x3_range = 1:256
    size_factor = 1
    up_thresh = 30000

end

D=load_untouch_nii(sprintf('%s',data_path));
im=double(D.img(x1_range,x2_range,x3_range,:));
im(im<0)=0;
im(im>up_thresh) = up_thresh;

if size_factor~=1 
    s_in = size(im);
    s_out = round(s_in*size_factor);
    im = resizeMatrix(im,s_out,'cubic');
end
    
im_out=im;
