classdef OrderedStartingPoints < StartingPoints
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        jump (1,1) double = 5;
    end

    methods

        function updateJump(obj,newJump)
            arguments
                obj
                newJump (1,1) double    
            end
            obj.jump = newJump;
            obj.method = sprintf('spDISTordered_fs%d', obj.jump);
            obj.methodTitle = sprintf('spDISTordered-fs%d', obj.jump);
                      
        end
        
        function select(obj,tag)
            %find indices of all voxels inside the ROI-SP
            [yIndex, xIndex, zIndex] = ind2sub(Configuration3D(tag).trueSize, obj.indexThresholded); 
            obj.y = yIndex(1:obj.jump:end);
            obj.x = xIndex(1:obj.jump:end);
            obj.z = zIndex(1:obj.jump:end);
        end
    end
end