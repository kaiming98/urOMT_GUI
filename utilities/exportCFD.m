function exportCFD(tag,varargin)
%EXPORTCFD Summary of this function goes here
%   Detailed explanation goes here
    
    p = inputParser;
    addRequired(p,'tag',@ischar)
    defaultFormat = 'VTK';
    addParameter(p,'format',defaultFormat);
    parse(p,tag,varargin{:});


    addpath('./Sensitivities', './analyzeFlows', genpath('./utilities'))

    configuration = Configuration3D(tag);
    configuration.loadVolume;
    
    ti = configuration.timeInitial;
    tj = configuration.timeJump;
    tf = configuration.timeFinal;
    nt = configuration.nt;
    
    xrange = 1:configuration.trueSize(2);
    yrange = 1:configuration.trueSize(1);
    zrange = 1:configuration.trueSize(3);
    [x, y, z] = meshgrid(xrange,yrange,zrange);
    iTime = 0;
    
    for t = ti:tj:tf %iterate over images
    
        U = importdata(sprintf('%s/%s_v_%d_%d_t_%d.mat', ...
            configuration.pathOutput, configuration.tag, ...
            t, t + tj, (t - ti) / tj + 1));
        
        rho = importdata(sprintf('%s/%s_d_%d_%d_t_%d.mat', ...
            configuration.pathOutput, configuration.tag, ...
            t, t + tj, (t - ti) / tj + 1));
    
        U = reshape(U, configuration.dimension * prod(configuration.n), configuration.nt);
        rho = reshape(rho, configuration.n');

        eps = 0.001;
        [w2, w1, w3] = gradient(log(rho + 2 * eps));
        du = configuration.sigma * [w1(:), w2(:), w3(:)];
        diffusiveSpeed = reshape(sqrt(sum(du.^2, 2)), configuration.trueSize);
    
            for t2 = 1:nt %iterate over numerical timesteps
                iTime = iTime + 1;
                u = reshape(U(:, t2), [], 3);
                v1 = reshape(u(:, 1), configuration.trueSize);
                v2 = reshape(u(:, 2), configuration.trueSize);
                v3 = reshape(u(:, 3), configuration.trueSize);
                
                velocityMagnitude = reshape(sqrt(sum(u.^2, 2)), configuration.trueSize);
                
                Peclet = velocityMagnitude./diffusiveSpeed;
                Peclet(isinf(Peclet)) = 0;

                switch upper(p.Results.format)
                    case 'PLT'
                        % create tecplot data structure
                        tdata = [];
                        tdata.title = sprintf('%s_%d_%d_t_%d.mat', configuration.tag, ...
                            t, t + tj, (t - ti) / tj + 1);
                        tdata.Nvar = 9;
                        tdata.vformat = 1*ones(1,tdata.Nvar); %float precision for each variable
                        tdata.varnames = {'x','y','z','Vx','Vy','Vz','density','Peclet','mask'};
                        tdata.cubes(1).x = x;
                        tdata.cubes(1).y = y;
                        tdata.cubes(1).z = z;
                        tdata.cubes(1).v(1,:,:,:) = v1;
                        tdata.cubes(1).v(2,:,:,:) = v2;
                        tdata.cubes(1).v(3,:,:,:) = v3;
                        tdata.cubes(1).v(4,:,:,:) = rho;
                        tdata.cubes(1).v(5,:,:,:) = Peclet;
                        tdata.cubes(1).v(6,:,:,:) = configuration.mask.contents;
                        tdata.cubes.solutiontime = t + (t2-1)/10;
                        filename = sprintf('tecplot/%s_%.1f.plt',tag,tdata.cubes.solutiontime);
                        mat2tecplot(tdata, filename);

                    case 'VTK'
                        filename = sprintf('vtk/%s_%03d.vtk',tag,iTime);
                        vtkwrite(filename,'structured_grid',x,y,z, ...
                            'vectors','vector_field',v1,v2,v3, ...
                            'scalars','density',rho, ...
                            'scalars','Peclet',Peclet, ...
                            'scalars','mask',configuration.mask.contents); 
                
                end

          
            end

    end

end