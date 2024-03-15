function  nii2vtk(pathData)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
data  = load_untouch_nii(pathData);
vtkwrite(sprintf('%s.vtk', pathData(1:end-4)), 'structured_points',...
    'mask', data.img);
end

