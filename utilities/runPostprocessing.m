function [pathlines, eulerian] = runPostprocessing(tag)
%RUNPOSTPROCESSING Lagrangian representation and visualization of rOMT results
%   Detailed explanation goes here

    % Setup 
    addpath('./Sensitivities', './analyzeFlows', genpath('./utilities'))
    %{
    % load scalar map analysis constants
    analysisConstants = loadJsonFile('scalarMap.json');
    
    % Calculate Lagrangian specification
    lagrangianSpecification = LagrangianSpecification(tag,analysisConstants);
    
    % Create new pathlines instance which contains density (speed,flux,...) along pathlines
    pathlines = Pathlines(tag,lagrangianSpecification);
    
    % Run QuickBundle algorithm in external python process
    pathlines.doQuickBundle(pathToPython);
    
    % Filter the N largest clusters.
    pathlines.filterClusters;
    
    % Create scalarMap instance
    scalarMap = ScalarMap(tag, pathlines, analysisConstants);
    
    % print useful information
    scalarMap.printInfo(lagrangianSpecification.startPoints);
    
    % Save results in NIfTI format 
    scalarMap.exportToNifti;
    
    % Plot velocity and Peclet fields within an anatomical region of interest
    %scalarMap.visualize3D;
    %}
    %% Here we save pathlines without bundling in vtk format for visualizing in a another software. 

    % load pathlines analysis constants
    analysisConstants = loadJsonFile('pathlines.json');
    
    % calculate Lagrangian specification
    lagrangianSpecification = LagrangianSpecification(tag,analysisConstants);
    
    % Create new pathlines instance at higher density (i.e., smaller jump value)
    pathlines = Pathlines(tag,lagrangianSpecification);
    
    % print useful information
    pathlines.printInfo(lagrangianSpecification.startPoints);
    
    % calculate extra properties, including number of points, length, displacement field...
    pathlines.calculateExtraProperties;
    
    % make sure all properties are within the anatomical region of interest
    pathlines.getPropertiesInMask;
    
    % save results in mat file
    pathlines.exportToMat;

    % save pathlines results into vtk format
    pathlines.exportToVtk;

    % visualize pathlines/spdlines/pelines/velocity flux vectors in matlab
    % pathlines.visualize;

    %% Here we save the eulerian visualization into vtk format
    eulerian = EulerianSpecification(tag);
    
    % save eulerian results into nifti format
    eulerian.exportToNifti;

    % save results in mat file
    eulerian.exportToMat;

    % save eulerian results into vtk format
    eulerian.exportToVtk;

    % visualize eulerian speed/source in matlab
    % eulerian.visualize;

end











