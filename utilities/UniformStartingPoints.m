classdef UniformStartingPoints < StartingPoints
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        percentage (1,1) double = 40;
        total % total number of points used
    end

    methods

        function update(obj,newPercentage)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                newPercentage (1,1) double    
            end
            obj.percentage = newPercentage;
            obj.method = sprintf('spDISTunif_spPerc%d', obj.percentage);
            obj.methodTitle = sprintf('spDISTunif-spPerc%d', obj.percentage);     
        end
        
        function select(obj,tag)
            %find indices of all voxels inside the ROI-SP
            startPointMaskvol = sum(obj.mask.contents(:)); %volume of mask used to select start points
            obj.total = round(obj.percentage * startPointMaskvol / 100);

            if ~exist(sprintf('%s_%s_startingPoints_Percentage%d_number%d.mat',...
                    tag, obj.maskName, obj.percentage, obj.total), 'file')
                [indexSampled, ~] = datasample(obj.indexThresholded, obj.total, 'Replace', false);
                save(sprintf('%s_%s_startingPoints_Percentage%d_number%d.mat',...
                    tag, obj.maskName,obj.percentage, obj.total), 'indexSampled');
            else
                load(sprintf('%s_%s_startingPoints_Percentage%d_number%d.mat',...
                    tag, obj.maskName, obj.percentage, obj.total));
            end

            [obj.y, obj.x, obj.z] = ind2sub(Configuration3D(tag).trueSize, sort(indexSampled, 'ascend'));
        end
    end
    
end