function wrapper(pathData,pathMask,pathOutputFolder,addSource,timeInitial,timeFinal,alpha,isDataPacked)
%WRAPPER wrapper for running a rOMT/urOMT demo
%   inputs:
%   pathMask - path to mask file.
%   pathData - path to data file.
%   pathOutputFolder - path to the folder that save the outputs.
%   addSource - 1 if using urOMT; 0 if using rOMT.
%   alpha - parameter controlling the balance between the flow and the
%   source.
%   timeInitial - the first frame number.
%   timeFinal - the final frame number.
%   isDataPacked - 1 if the input is a 4D volume; 0 if the inputs is a
%   series of 3D images.

arguments
    pathData;
    pathMask;
    pathOutputFolder;
    addSource;
    timeInitial;
    timeFinal;
    alpha double = 10000;
    isDataPacked double = 1;
end

addpath('./Sensitivities', './analyzeFlows', genpath('./utilities'))
%% write parameters into a json file, json file name is jsonFileName+'.json'
obj.pathMask = pathMask;
obj.pathData = pathData;
obj.alpha = alpha;
obj.addSource = addSource;
obj.timeInitial = timeInitial;
obj.timeFinal = timeFinal;
obj.pathOutputFolder = pathOutputFolder;
obj.isDataPacked = isDataPacked;
encodedJson = jsonencode(obj,PrettyPrint=true);
jsonFileName = 'templete';
fid=fopen(sprintf('%s.json',jsonFileName),'w');
fprintf(fid, encodedJson); 
fclose(fid);

% %% run OMT 
% runRegularizedOMT(jsonFileName);
% process and visualize results
runPostprocessing(jsonFileName);

%% visualize lagrangian results
% visualize pathlines
visualizePathlines(jsonFileName);
% visualize speed-lines
visualizeSpeedlines(jsonFileName);
% visualize pe-lines
visualizePelines(jsonFileName);
% visualize velocity flux vectors
visualizeDisplacementVectors(jsonFileName);

%% visualize eulerian results
% visualize eulerian speed
visualizeEulSpeed(jsonFileName);
% visualize eulerian source
visualizeEulSource(jsonFileName);
end