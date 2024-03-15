classdef Configuration3D < Configuration
    % Subclass for handling three-dimensional data 
    
    properties
        Zc
    end
    
    methods
        
        function rho = advectionDiffusion(obj,rho0, u)

            if isempty(obj.B)
                obj.initializeGrid;
            end

            rho = zeros(prod(obj.n), obj.nt + 1);
            rho(:, 1) = rho0;

            u = reshape(u, 3 * prod(obj.n), obj.nt);
    
            for i = 2:obj.nt + 1
                
                % separate spatial velocity components per time point
                U1 = reshape(u(1:prod(obj.n), i - 1), obj.n');
                U2 = reshape(u(prod(obj.n) + 1:2 * prod(obj.n), i - 1), obj.n');
                U3 = reshape(u(2 * prod(obj.n) + 1:end, i - 1), obj.n');
                
                % sends the rho and x+dx, y+dy, z+dz (spatial translation
                % based on velocities)
                S = dTrilinears3d(rho(:, i - 1),...
                    obj.Xc + obj.dt * U1,...
                    obj.Yc + obj.dt * U2,...
                    obj.Zc + obj.dt * U3, ...
                    1, 1, 1, obj.boundaryConditions);
     
                switch obj.diffusionCoefficientType
                    case 'constant'
                        % advection step
                        % applies the S operator computer before which is the
                        % advection step in the advection diffusion equation. 
                        rho(:, i) = S * rho(:, i - 1);
                        % diffusion step
                        [rho(:, i), pcgflag] = pcg(obj.B, rho(:, i));
                        if pcgflag ~= 0
                            warning('MATLAB:pcgExitFlag', 'Warning: advecDiff.m >>> while finding rho(:,%d), pcg exit flag = %d', i, pcgflag)
                        end
                    case 'anisotropic'
                        diffusionResult = obj.calculateAnisotropicDiffusion(rho(:, i - 1));
                        rho(:, i) = S * rho(:, i - 1) + ...
                            obj.dt*diffusionResult;
                    case 'autoAnisotropic'
                        diffusionResult = obj.calculateAnisotropicDiffusion(rho(:, i - 1));
                        rho(:, i) = S * rho(:, i - 1) + ...
                            obj.dt*diffusionResult;   
                    case 'brightness'
                        diffusionResult = obj.calculateBrightnessDiffusion(rho(:, i - 1));
                        rho(:, i) = S * rho(:, i - 1) + ...
                            obj.dt*diffusionResult; 
                end
    
            end
    
            rho = rho(:, 2:end);
 
        end
        
        function rho = advectionDiffusionSource(obj, rho0, u, r)
            if isempty(obj.B)
                obj.initializeGrid;
            end
            
            rho = zeros(prod(obj.n), obj.nt + 1);
            rho(:, 1) = rho0;

            u = reshape(u, 3 * prod(obj.n), obj.nt);
            r = reshape(r, prod(obj.n), obj.nt);
            for i = 2:obj.nt + 1
                
                rhoSource = (1 + obj.dt * r(:, i - 1)) .* rho(:, i - 1);

                U1 = reshape(u(1:prod(obj.n), i - 1), obj.n');
                U2 = reshape(u(prod(obj.n) + 1:2 * prod(obj.n), i - 1), obj.n');
                U3 = reshape(u(2 * prod(obj.n) + 1:end, i - 1), obj.n');

                % sends the rho and x+dx, y+dy, z+dz (spatial translation
                % based on velocities)
                S = dTrilinears3d(rhoSource,...
                    obj.Xc + obj.dt * U1,...
                    obj.Yc + obj.dt * U2,...
                    obj.Zc + obj.dt * U3, ...
                    1, 1, 1, obj.boundaryConditions);
                
                switch obj.diffusionCoefficientType
                    case 'constant'
                        rho(:, i) = S * rhoSource;
                        [rho(:, i), pcgflag] = pcg(obj.B, rho(:, i));
                        if pcgflag ~= 0
                            warning('MATLAB:pcgExitFlag', 'Warning: advecDiff.m >>> while finding rho(:,%d), pcg exit flag = %d', i, pcgflag)
                        end
                    case 'anisotropic'
                        rho(:, i) = S * rhoSource + ...
                            obj.dt*obj.calculateAnisotropicDiffusion(rho(:, i - 1));
                    case 'autoAnisotropic'
                        diffusionResult = obj.calculateAnisotropicDiffusion(rho(:, i - 1));
                        rho(:, i) = S * rhoSource(:, i - 1) + ...
                            obj.dt*diffusionResult;   
                    case 'brightness'
                        diffusionResult = obj.calculateBrightnessDiffusion(rho(:, i - 1));
                        rho(:, i) = S * rhoSource(:, i - 1) + ...
                            obj.dt*diffusionResult; 
                end

            end

            rho = rho(:, 2:end);
        end
                
    end

    methods (Access = protected)

        function getCellCenteredData(obj)

            cellDimension1 = ones(obj.trueSize(1),1);
            cellDimension2 = ones(obj.trueSize(2),1);
            cellDimension3 = ones(obj.trueSize(3),1);
                    
            [obj.Xc, obj.Yc, obj.Zc] = getCellCenteredGrid(cellDimension1, cellDimension2, cellDimension3);    

            boundaryConditionType = {'ccn' 'ccn' 'ccn'}; % Cell Centered Nuemann
            obj.Grad = getCellCenteredGradMatrix(boundaryConditionType, cellDimension1, cellDimension2, cellDimension3);
        end

        function M = calculateM(obj, RHO0, U, R)
        % calculates intermediary matrix M used in getVelocity()
        % 
        % M.S is a matrix, that when multiplied against a density
        % distribution, returns a new, advected density disitrubution. 
        %
        % M.Tx, M.Ty, M.Tz are spatial derivatives (should be renamed
        % because T refers to density, which is the variable used in the
        % dTrilinears3d function but not in the rest of this code. 
        % perhaps rename to getAdvectionIntermediaryMatrices(). 

            for k = 1:obj.nt
                U1 = reshape(U(1:prod(obj.n), k), obj.n');
                U2 = reshape(U(prod(obj.n) + 1:2 * prod(obj.n), k), obj.n');
                U3 = reshape(U(2 * prod(obj.n) + 1:end, k), obj.n');
                if ~obj.addSource
                    [M.S{k}, M.Tx{k}, M.Ty{k}, M.Tz{k}] = dTrilinears3d(RHO0(:, k), ...
                        obj.Xc + obj.dt * U1, obj.Yc + obj.dt * U2, obj.Zc + obj.dt * U3,...
                        1, 1, 1, obj.boundaryConditions);
                else
                    M.R{k} = 1 + obj.dt * R(:, k);
                    M.Rho{k} = RHO0(:, k);
                    [M.S{k}, M.Tx{k}, M.Ty{k}, M.Tz{k}] = dTrilinears3d(M.R{k}.*RHO0(:, k), ...
                        obj.Xc + obj.dt * U1, obj.Yc + obj.dt * U2, obj.Zc + obj.dt * U3,...
                        1, 1, 1, obj.boundaryConditions);                    
                end
            end            
        end
        
        function diffusionResult = calculateAnisotropicDiffusion(obj,rho)
            % CALCULATEANISOTROPICDIFFUSION calculate the diffusion term
            % \nabla\cdot(\sigma(abs(\nabla rho))\nabla rho).
            % Here we only consider cell-centered neumann boundary
            % conditions(ccn).
            Rho = reshape(rho,obj.trueSize);
            % boundary condition
            bc = {'ccn','ccn','ccn'};
            % spatial grid size, h = [dx,dy].
            h = {ones(1,obj.trueSize(2)-1), ones(1,obj.trueSize(1)-1), ones(1,obj.trueSize(3)-1)}; 
            differenceXMinus = obj.doXMinusDifference(Rho,bc{1},h{1});
            differenceXPlus = obj.doXPlusDifference(Rho,bc{1},h{1});
            differenceYMinus = obj.doYMinusDifference(Rho,bc{2},h{2});
            differenceYPlus = obj.doYPlusDifference(Rho,bc{2},h{2});
            differenceZMinus = obj.doZMinusDifference(Rho,bc{3},h{3});
            differenceZPlus = obj.doZPlusDifference(Rho,bc{3},h{3});
            
            diffusionResult = obj.doXMinusDifference(obj.sigmaFunctional(differenceXPlus.^2 + ...
                obj.minmod(differenceYMinus,differenceYPlus).^2 + ...
                obj.minmod(differenceZMinus,differenceZPlus).^2).*differenceXPlus,bc{1},h{1}) + ...
                obj.doYMinusDifference(obj.sigmaFunctional(differenceYPlus.^2 + ...
                obj.minmod(differenceXMinus,differenceXPlus).^2 + ...
                obj.minmod(differenceZMinus,differenceZPlus).^2).*differenceYPlus,bc{2},h{2}) + ...
                obj.doZMinusDifference(obj.sigmaFunctional(differenceZPlus.^2 + ...
                obj.minmod(differenceXMinus,differenceXPlus).^2 + ...
                obj.minmod(differenceYMinus,differenceYPlus).^2).*differenceZPlus,bc{3},h{3});
            diffusionResult = diffusionResult(:);
        end
        
        function diffusionJacobian = calculateDiffusionJacobian(obj,rho)
            %CALCULATEDIFFUSIONJACOBIAN calculate the jacobian matrix of
            %anisotropic diffuion term D w.r.t the density rho.
            % rho(i,j,k) i index for y axis; j index for x axis; k index 
            % for z axis.
            Rho = reshape(rho,obj.trueSize);
            % boundary condition
            bc = {'ccn','ccn','ccn'};
            % spatial grid size, h = [dx,dy,dz].
            h = {ones(1,obj.trueSize(2)-1), ones(1,obj.trueSize(1)-1), ones(1,obj.trueSize(3)-1)}; 
            differenceXMinus = obj.doXMinusDifference(Rho,bc{1},h{1});
            differenceXPlus = obj.doXPlusDifference(Rho,bc{1},h{1});
            differenceYMinus = obj.doYMinusDifference(Rho,bc{2},h{2});
            differenceYPlus = obj.doYPlusDifference(Rho,bc{2},h{2});
            differenceZMinus = obj.doZMinusDifference(Rho,bc{3},h{3});
            differenceZPlus = obj.doZPlusDifference(Rho,bc{3},h{3});

            % gradient rho at different place
            gradientXPlus = differenceXPlus.^2 + obj.minmod(differenceYMinus,differenceYPlus).^2 + ...
                obj.minmod(differenceZMinus,differenceZPlus).^2;
            gradientXMinus = zeros(size(gradientXPlus)); 
            gradientXMinus(:,2:end,:) = gradientXPlus(:,1:end-1,:);
            gradientXMinus(:,1,:) = obj.minmod(differenceYMinus(:,1,:),differenceYPlus(:,1,:)).^2 + ...
                obj.minmod(differenceZMinus(:,1,:),differenceZPlus(:,1,:)).^2; 
            
            gradientYPlus = differenceYPlus.^2 + obj.minmod(differenceXMinus,differenceXPlus).^2 + ...
                obj.minmod(differenceZMinus,differenceZPlus).^2;
            gradientYMinus = zeros(size(gradientYPlus)); 
            gradientYMinus(2:end,:,:) = gradientYPlus(1:end-1,:,:);
            gradientYMinus(1,:,:) = obj.minmod(differenceXMinus(1,:,:),differenceXPlus(1,:,:)).^2 + ...
                obj.minmod(differenceZMinus(1,:,:),differenceZPlus(1,:,:)).^2;
            
            gradientZPlus = differenceZPlus.^2 + obj.minmod(differenceXMinus,differenceXPlus).^2 + ...
                obj.minmod(differenceYMinus,differenceYPlus).^2;
            gradientZMinus = zeros(size(gradientZPlus)); 
            gradientZMinus(:,:,2:end) = gradientZPlus(:,:,1:end-1);
            gradientZMinus(:,:,1) = obj.minmod(differenceXMinus(:,:,1),differenceXPlus(:,:,1)).^2 + ...
                obj.minmod(differenceYMinus(:,:,1),differenceYPlus(:,:,1)).^2;
            % diffusion coefficient at different place
            sigmaXPlus = obj.sigmaFunctional(gradientXPlus);
            sigmaXMinus = obj.sigmaFunctional(gradientXMinus);            
            sigmaYPlus = obj.sigmaFunctional(gradientYPlus);
            sigmaYMinus = obj.sigmaFunctional(gradientYMinus);            
            sigmaZPlus = obj.sigmaFunctional(gradientZPlus);
            sigmaZMinus = obj.sigmaFunctional(gradientZMinus);
            % diffusion coefficient derivative at different place
            dsigmaXPlus = obj.sigmaFunctionalDerivative(gradientXPlus);
            dsigmaXMinus = obj.sigmaFunctionalDerivative(gradientXMinus);
            dsigmaYPlus = obj.sigmaFunctionalDerivative(gradientYPlus);
            dsigmaYMinus = obj.sigmaFunctionalDerivative(gradientYMinus);
            dsigmaZPlus = obj.sigmaFunctionalDerivative(gradientZPlus);
            dsigmaZMinus = obj.sigmaFunctionalDerivative(gradientZMinus);
            RhoGhost = obj.addGhostPoints(Rho);

            % \partial D_{i,j,k}/\partial \rho_{i,j,k}
            fijkRhoijk = -sigmaXPlus - sigmaXMinus + differenceXPlus.*dsigmaXPlus.* ...
                (-2*differenceXPlus + obj.minmod2Derivative(RhoGhost(3:end,2:end-1,2:end-1),Rho,RhoGhost(1:end-2,2:end-1,2:end-1),2) + ...
                obj.minmod2Derivative(RhoGhost(2:end-1,2:end-1,3:end),Rho,RhoGhost(2:end-1,2:end-1,1:end-2),2)) - ...
                2*differenceXMinus.^2.*dsigmaXMinus - sigmaYPlus - sigmaYMinus + ...
                differenceYPlus.*dsigmaYPlus.*( -2*differenceYPlus + obj.minmod2Derivative(RhoGhost(2:end-1,3:end,2:end-1),Rho,RhoGhost(2:end-1,1:end-2,2:end-1),2) + ...
                obj.minmod2Derivative(RhoGhost(2:end-1,2:end-1,3:end),Rho,RhoGhost(2:end-1,2:end-1,1:end-2),2)) - ...
                2*differenceYMinus.^2.*dsigmaYMinus - sigmaZPlus - sigmaZMinus + ...
                differenceZPlus.*dsigmaZPlus.*( -2*differenceZPlus + obj.minmod2Derivative(RhoGhost(2:end-1,3:end,2:end-1),Rho,RhoGhost(2:end-1,1:end-2,2:end-1),2) + ...
                obj.minmod2Derivative(RhoGhost(3:end,2:end-1,2:end-1),Rho,RhoGhost(1:end-2,2:end-1,2:end-1),2)) - ...
                2*differenceZMinus.^2.*dsigmaZMinus;
            % \partial D_{i,j,k}/\partial \rho_{i,j+1,k}
            fijkRhoijPlusk = sigmaXPlus + 2*(differenceXPlus.^2).*dsigmaXPlus + ...
                (differenceYPlus.*dsigmaYPlus + differenceZPlus.*dsigmaZPlus).* ...
                obj.minmod2Derivative(RhoGhost(2:end-1,3:end,2:end-1),Rho,RhoGhost(2:end-1,1:end-2,2:end-1),1);
            % \partial D_{i,j,k}/\partial \rho_{i+1,j,k}
            fijkRhoiPlusjk = sigmaYPlus + 2*(differenceYPlus.^2).*dsigmaYPlus + ...
                (differenceXPlus.*dsigmaXPlus + differenceZPlus.*dsigmaZPlus).* ...
                obj.minmod2Derivative(RhoGhost(3:end,2:end-1,2:end-1),Rho,RhoGhost(1:end-2,2:end-1,2:end-1),1);
            % \partial D_{i,j,k}/\partial \rho_{i,j,k+1}
            fijkRhoijkPlus = sigmaZPlus + 2*(differenceZPlus.^2).*dsigmaZPlus + ...
                (differenceXPlus.*dsigmaXPlus + differenceYPlus.*dsigmaYPlus).* ...
                obj.minmod2Derivative(RhoGhost(2:end-1,2:end-1,3:end),Rho,RhoGhost(2:end-1,2:end-1,1:end-2),1);
            % \partial D_{i,j,k}/\partial \rho_{i,j-1,k}
            fijkRhoijMinusk = sigmaXMinus - differenceXMinus.*dsigmaXMinus.*(-2*differenceXMinus + ...
                obj.minmod2Derivative(RhoGhost(3:end,1:end-2,2:end-1),RhoGhost(2:end-1,1:end-2,2:end-1),RhoGhost(1:end-2,1:end-2,2:end-1),2) + ...
                obj.minmod2Derivative(RhoGhost(2:end-1,1:end-2,3:end),RhoGhost(2:end-1,1:end-2,2:end-1),RhoGhost(2:end-1,1:end-2,1:end-2),2)) + ...
                (differenceYPlus.*dsigmaYPlus + differenceZPlus.*dsigmaZPlus).* ...
                obj.minmod2Derivative(RhoGhost(2:end-1,3:end,2:end-1),Rho,RhoGhost(2:end-1,1:end-2,2:end-1),3);
            % \partial D_{i,j,k}/\partial \rho_{i-1,j,k}
            fijkRhoiMinusjk = sigmaYMinus - differenceYMinus.*dsigmaYMinus.*(-2*differenceYMinus + ...
                obj.minmod2Derivative(RhoGhost(1:end-2,3:end,2:end-1),RhoGhost(1:end-2,2:end-1,2:end-1),RhoGhost(1:end-2,1:end-2,2:end-1),2) + ...
                obj.minmod2Derivative(RhoGhost(1:end-2,2:end-1,3:end),RhoGhost(1:end-2,2:end-1,2:end-1),RhoGhost(1:end-2,2:end-1,1:end-2),2)) + ...
                (differenceXPlus.*dsigmaXPlus + differenceZPlus.*dsigmaZPlus).* ...
                obj.minmod2Derivative(RhoGhost(3:end,2:end-1,2:end-1),Rho,RhoGhost(1:end-2,2:end-1,2:end-1),3);
            % \partial D_{i,j,k}/\partial \rho_{i,j,k-1}
            fijkRhoijkMinus = sigmaZMinus - differenceZMinus.*dsigmaZMinus.*(-2*differenceZMinus + ...
                obj.minmod2Derivative(RhoGhost(3:end,2:end-1,1:end-2),RhoGhost(2:end-1,2:end-1,1:end-2),RhoGhost(1:end-2,2:end-1,1:end-2),2) + ...
                obj.minmod2Derivative(RhoGhost(2:end-1,3:end,1:end-2),RhoGhost(2:end-1,2:end-1,1:end-2),RhoGhost(2:end-1,1:end-2,1:end-2),2)) + ...
                (differenceXPlus.*dsigmaXPlus + differenceYPlus.*dsigmaYPlus).* ...
                obj.minmod2Derivative(RhoGhost(2:end-1,2:end-1,3:end),Rho,RhoGhost(2:end-1,2:end-1,1:end-2),3);
            % \partial D_{i,j,k}/\partial \rho_{i-1,j+1,k}
            fijkRhoiMinusjPlusk = -differenceYMinus.*obj.minmod2Derivative(RhoGhost(1:end-2,3:end,2:end-1),RhoGhost(1:end-2,2:end-1,2:end-1),RhoGhost(1:end-2,1:end-2,2:end-1),1).*dsigmaYMinus;
            % \partial D_{i,j,k}/\partial \rho_{i+1,j-1,k}
            fijkRhoiPlusjMinusk = -differenceXMinus.*obj.minmod2Derivative(RhoGhost(3:end,1:end-2,2:end-1),RhoGhost(2:end-1,1:end-2,2:end-1),RhoGhost(1:end-2,1:end-2,2:end-1),1).*dsigmaXMinus;
            % \partial D_{i,j,k}/\partial \rho_{i-1,j,k+1}
            fijkRhoiMinusjkPlus = -differenceYMinus.*obj.minmod2Derivative(RhoGhost(1:end-2,2:end-1,3:end),RhoGhost(1:end-2,2:end-1,2:end-1),RhoGhost(1:end-2,2:end-1,1:end-2),1).*dsigmaYMinus;
            % \partial D_{i,j,k}/\partial \rho_{i+1,j,k-1}
            fijkRhoiPlusjkMinus = -differenceZMinus.*obj.minmod2Derivative(RhoGhost(3:end,2:end-1,1:end-2),RhoGhost(2:end-1,2:end-1,1:end-2),RhoGhost(1:end-2,2:end-1,1:end-2),1).*dsigmaZMinus;
            % \partial D_{i,j,k}/\partial \rho_{i,j+1,k-1}
            fijkRhoijPluskMinus = -differenceZMinus.*obj.minmod2Derivative(RhoGhost(2:end-1,3:end,1:end-2),RhoGhost(2:end-1,2:end-1,1:end-2),RhoGhost(2:end-1,1:end-2,1:end-2),1).*dsigmaZMinus;
            % \partial D_{i,j,k}/\partial \rho_{i,j-1,k+1}
            fijkRhoijMinuskPlus = -differenceXMinus.*obj.minmod2Derivative(RhoGhost(2:end-1,1:end-2,3:end),RhoGhost(2:end-1,1:end-2,2:end-1),RhoGhost(2:end-1,1:end-2,1:end-2),1).*dsigmaXMinus;
            % \partial D_{i,j,k}/\partial \rho_{i-1,j-1,k}
            fijkRhoiMinusjMinusk = -differenceYMinus.*obj.minmod2Derivative(RhoGhost(1:end-2,3:end,2:end-1),RhoGhost(1:end-2,2:end-1,2:end-1),RhoGhost(1:end-2,1:end-2,2:end-1),3).*dsigmaYMinus - ...
                differenceXMinus.*obj.minmod2Derivative(RhoGhost(3:end,1:end-2,2:end-1),RhoGhost(2:end-1,1:end-2,2:end-1),RhoGhost(1:end-2,1:end-2,2:end-1),3).*dsigmaXMinus;
            % \partial D_{i,j,k}/\partial \rho_{i,j-1,k-1}
            fijkRhoijMinuskMinus = -differenceZMinus.*obj.minmod2Derivative(RhoGhost(2:end-1,3:end,1:end-2),RhoGhost(2:end-1,2:end-1,1:end-2),RhoGhost(2:end-1,1:end-2,1:end-2),3).*dsigmaZMinus - ...
                differenceXMinus.*obj.minmod2Derivative(RhoGhost(2:end-1,1:end-2,3:end),RhoGhost(2:end-1,1:end-2,2:end-1),RhoGhost(2:end-1,1:end-2,1:end-2),3).*dsigmaXMinus;
            % \partial D_{i,j,k}/\partial \rho_{i-1,j,k-1}
            fijkRhoiMinusjkMinus = -differenceYMinus.*obj.minmod2Derivative(RhoGhost(1:end-2,2:end-1,3:end),RhoGhost(1:end-2,2:end-1,2:end-1),RhoGhost(1:end-2,2:end-1,1:end-2),3).*dsigmaYMinus - ...
                differenceZMinus.*obj.minmod2Derivative(RhoGhost(3:end,2:end-1,1:end-2),RhoGhost(2:end-1,2:end-1,1:end-2),RhoGhost(1:end-2,2:end-1,1:end-2),3).*dsigmaZMinus;
            
            %Put information into the Jacobian sparse matrix jacobian(D)
            nRho = prod(obj.trueSize);
            ii1 = 1:nRho; jj1 = 1:nRho; uu1 = fijkRhoijk(:);
            A1 = sparse(ii1,jj1,uu1,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),1:obj.trueSize(1)-1,1:obj.trueSize(3));
            ind2 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii2 = ii1(ind2); jj2 = ii2+1; uu2 = fijkRhoiPlusjk(:); uu2 = uu2(ind2);
            A2 = sparse(ii2,jj2,uu2,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),2:obj.trueSize(1),1:obj.trueSize(3));
            ind3 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii3 = ii1(ind3); jj3 = ii3-1; uu3 = fijkRhoiMinusjk(:); uu3 = uu3(ind3);
            A3 = sparse(ii3,jj3,uu3,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2)-1,1:obj.trueSize(1),1:obj.trueSize(3));
            ind4 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii4 = ii1(ind4); jj4 = ii4+obj.trueSize(1); uu4 = fijkRhoijPlusk(:); uu4 = uu4(ind4);
            A4 = sparse(ii4,jj4,uu4,nRho,nRho);
            
            [X,Y,Z] = meshgrid(2:obj.trueSize(2),1:obj.trueSize(1),1:obj.trueSize(3));
            ind5 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii5 = ii1(ind5); jj5 = ii5-obj.trueSize(1); uu5 = fijkRhoijMinusk(:); uu5 = uu5(ind5);
            A5 = sparse(ii5,jj5,uu5,nRho,nRho);
            
            [X,Y,Z] = meshgrid(2:obj.trueSize(2),1:obj.trueSize(1)-1,1:obj.trueSize(3));
            ind6 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii6 = ii1(ind6); jj6 = ii6-obj.trueSize(1)+1; uu6 = fijkRhoiPlusjMinusk(:); uu6 = uu6(ind6);
            A6 = sparse(ii6,jj6,uu6,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2)-1,2:obj.trueSize(1),1:obj.trueSize(3));
            ind7 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii7 = ii1(ind7); jj7 = ii7+obj.trueSize(1)-1; uu7 = fijkRhoiMinusjPlusk(:); uu7 = uu7(ind7);
            A7 = sparse(ii7,jj7,uu7,nRho,nRho);
            
            [X,Y,Z] = meshgrid(2:obj.trueSize(2),2:obj.trueSize(1),1:obj.trueSize(3));
            ind8 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii8 = ii1(ind8); jj8 = ii8-obj.trueSize(1)-1; uu8 = fijkRhoiMinusjMinusk(:); uu8 = uu8(ind8);
            A8 = sparse(ii8,jj8,uu8,nRho,nRho);
            
            XYSize = obj.trueSize(1)*obj.trueSize(2);
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),1:obj.trueSize(1),1:obj.trueSize(3)-1);
            ind9 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii9 = ii1(ind9); jj9 = ii9+XYSize; uu9 = fijkRhoijkPlus(:); uu9 = uu9(ind9);
            A9 = sparse(ii9,jj9,uu9,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),1:obj.trueSize(1),2:obj.trueSize(3));
            ind10 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii10 = ii1(ind10); jj10 = ii10-XYSize; uu10 = fijkRhoijkMinus(:); uu10 = uu10(ind10);
            A10 = sparse(ii10,jj10,uu10,nRho,nRho);
            
            [X,Y,Z] = meshgrid(2:obj.trueSize(2),1:obj.trueSize(1),2:obj.trueSize(3));
            ind11 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii11 = ii1(ind11); jj11 = ii11-XYSize-obj.trueSize(1); uu11 = fijkRhoijMinuskMinus(:); uu11 = uu11(ind11);
            A11 = sparse(ii11,jj11,uu11,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),2:obj.trueSize(1),2:obj.trueSize(3));
            ind12 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii12 = ii1(ind12); jj12 = ii12-XYSize-1; uu12 = fijkRhoiMinusjkMinus(:); uu12 = uu12(ind12);
            A12 = sparse(ii12,jj12,uu12,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),1:obj.trueSize(1)-1,2:obj.trueSize(3));
            ind13 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii13 = ii1(ind13); jj13 = ii13-XYSize+1; uu13 = fijkRhoiPlusjkMinus(:); uu13 = uu13(ind13);
            A13 = sparse(ii13,jj13,uu13,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2)-1,1:obj.trueSize(1),2:obj.trueSize(3));
            ind14 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii14 = ii1(ind14); jj14 = ii14-XYSize+obj.trueSize(1); uu14 = fijkRhoijPluskMinus(:); uu14 = uu14(ind14);
            A14 = sparse(ii14,jj14,uu14,nRho,nRho);
            
            [X,Y,Z] = meshgrid(2:obj.trueSize(2),1:obj.trueSize(1),1:obj.trueSize(3)-1);
            ind15 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii15 = ii1(ind15); jj15 = ii15+XYSize-obj.trueSize(1); uu15 = fijkRhoijMinuskPlus(:); uu15 = uu15(ind15);
            A15 = sparse(ii15,jj15,uu15,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),2:obj.trueSize(1),1:obj.trueSize(3)-1);
            ind16 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii16 = ii1(ind16); jj16 = ii16+XYSize-1; uu16 = fijkRhoiMinusjkPlus(:); uu16 = uu16(ind16);
            A16 = sparse(ii16,jj16,uu16,nRho,nRho);
            
            diffusionJacobian = A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + ...
                                A10 + A11 + A12 + A13 + A14 + A15 + A16;            
        end
        
        function diffusionResult = calculateBrightnessDiffusion(obj,rho)
            % CALCULATEBRIGHTNESSDIFFUSION calculate the brightness
            % diffusion term, i.e. \nabla\cdot(\rho\nabla\rho)
            % Here we only consider cell-centered neumann boundary
            % conditions(ccn).
            Rho = reshape(rho,obj.trueSize);
            % boundary condition
            bc = {'ccn','ccn','ccn'}; 
            h = {ones(1,obj.trueSize(2)-1), ones(1,obj.trueSize(1)-1), ones(1,obj.trueSize(3)-1)}; 
            differenceXPlus = obj.doXPlusDifference(Rho,bc{1},h{1});
            differenceYPlus = obj.doYPlusDifference(Rho,bc{2},h{2});
            differenceZPlus = obj.doZPlusDifference(Rho,bc{3},h{3});  
            RhoGhost = obj.addGhostPoints(Rho);
            RhoXPlus = RhoGhost(2:end-1,3:end,2:end-1);
            RhoYPlus = RhoGhost(3:end,2:end-1,2:end-1);
            RhoZPlus = RhoGhost(2:end-1,2:end-1,3:end);
            diffusionResult = obj.doXMinusDifference(obj.sigmaFunctional((Rho+RhoXPlus)/2).*differenceXPlus,bc{1},h{1}) + ...
                obj.doYMinusDifference(obj.sigmaFunctional((Rho+RhoYPlus)/2).*differenceYPlus,bc{2},h{2}) + ...
                obj.doZMinusDifference(obj.sigmaFunctional((Rho+RhoZPlus)/2).*differenceZPlus,bc{3},h{3});
            diffusionResult = obj.sigma*diffusionResult(:);
        end
        
        function diffusionJacobian = calculateBDiffusionJacobian(obj,rho)
            % CALCULATEBDIFFUSIONJACOBIAN calculate jacobian of brightness
            % diffusion w.r.t density rho.
            
            Rho = reshape(rho,obj.trueSize);
            % boundary condition
            bc = {'ccn','ccn','ccn'}; 
            h = {ones(1,obj.trueSize(2)-1), ones(1,obj.trueSize(1)-1), ones(1,obj.trueSize(3)-1)}; 
            differenceXPlus = obj.doXPlusDifference(Rho,bc{1},h{1});
            differenceYPlus = obj.doYPlusDifference(Rho,bc{2},h{2});
            differenceZPlus = obj.doZPlusDifference(Rho,bc{3},h{3}); 
            differenceXMinus = obj.doXMinusDifference(Rho,bc{1},h{1});
            differenceYMinus = obj.doYMinusDifference(Rho,bc{2},h{2});
            differenceZMinus = obj.doZMinusDifference(Rho,bc{3},h{3}); 
            RhoGhost = obj.addGhostPoints(Rho);
            meanRhoXPlus = (RhoGhost(2:end-1,3:end,2:end-1)+Rho)/2;
            meanRhoYPlus = (RhoGhost(3:end,2:end-1,2:end-1)+Rho)/2;
            meanRhoZPlus = (RhoGhost(2:end-1,2:end-1,3:end)+Rho)/2;
            meanRhoXMinus = (RhoGhost(2:end-1,1:end-2,2:end-1)+Rho)/2;
            meanRhoYMinus = (RhoGhost(1:end-2,2:end-1,2:end-1)+Rho)/2;
            meanRhoZMinus = (RhoGhost(2:end-1,2:end-1,1:end-2)+Rho)/2;
            
            
            DijkRhoijk = obj.sigmaFunctionalDerivative(meanRhoXPlus).*differenceXPlus/2 - ...
                obj.sigmaFunctional(meanRhoXPlus) - obj.sigmaFunctional(meanRhoXMinus) - ...
                obj.sigmaFunctionalDerivative(meanRhoXMinus).*differenceXMinus/2 + ...
                obj.sigmaFunctionalDerivative(meanRhoYPlus).*differenceYPlus/2 - ...
                obj.sigmaFunctional(meanRhoYPlus) - obj.sigmaFunctional(meanRhoYMinus) - ...
                obj.sigmaFunctionalDerivative(meanRhoYMinus).*differenceYMinus/2 + ...
                obj.sigmaFunctionalDerivative(meanRhoZPlus).*differenceZPlus/2 - ...
                obj.sigmaFunctional(meanRhoZPlus) - obj.sigmaFunctional(meanRhoZMinus) - ...
                obj.sigmaFunctionalDerivative(meanRhoZMinus).*differenceZMinus/2;
            
            DijkRhoiPlusjk = obj.sigmaFunctionalDerivative(meanRhoYPlus).*differenceYPlus/2 + ...
                obj.sigmaFunctional(meanRhoYPlus);           
            DijkRhoiMinusjk = -obj.sigmaFunctionalDerivative(meanRhoYMinus).*differenceYMinus/2 + ...
                obj.sigmaFunctional(meanRhoYMinus);            
            DijkRhoijPlusk = obj.sigmaFunctionalDerivative(meanRhoXPlus).*differenceXPlus/2 + ...
                obj.sigmaFunctional(meanRhoXPlus);            
            DijkRhoijMinusk = -obj.sigmaFunctionalDerivative(meanRhoXMinus).*differenceXMinus/2 + ...
                obj.sigmaFunctional(meanRhoXMinus);
            DijkRhoijkPlus = obj.sigmaFunctionalDerivative(meanRhoZPlus).*differenceZPlus/2 + ...
                obj.sigmaFunctional(meanRhoZPlus);
            DijkRhoijkMinus = -obj.sigmaFunctionalDerivative(meanRhoZMinus).*differenceZMinus/2 + ...
                obj.sigmaFunctional(meanRhoZMinus);
            
            
            nRho = prod(obj.trueSize);
            ii1 = 1:nRho; jj1 = 1:nRho; uu1 = DijkRhoijk(:);
            A1 = sparse(ii1,jj1,uu1,nRho,nRho);
            
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),1:obj.trueSize(1)-1,1:obj.trueSize(3));
            ind2 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii2 = ii1(ind2); jj2 = ii2+1; uu2 = DijkRhoiPlusjk(:); uu2 = uu2(ind2);
            A2 = sparse(ii2,jj2,uu2,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),2:obj.trueSize(1),1:obj.trueSize(3));
            ind3 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii3 = ii1(ind3); jj3 = ii3-1; uu3 = DijkRhoiMinusjk(:); uu3 = uu3(ind3);
            A3 = sparse(ii3,jj3,uu3,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2)-1,1:obj.trueSize(1),1:obj.trueSize(3));
            ind4 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii4 = ii1(ind4); jj4 = ii4+obj.trueSize(1); uu4 = DijkRhoijPlusk(:); uu4 = uu4(ind4);
            A4 = sparse(ii4,jj4,uu4,nRho,nRho);
            
            [X,Y,Z] = meshgrid(2:obj.trueSize(2),1:obj.trueSize(1),1:obj.trueSize(3));
            ind5 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii5 = ii1(ind5); jj5 = ii5-obj.trueSize(1); uu5 = DijkRhoijMinusk(:); uu5 = uu5(ind5);
            A5 = sparse(ii5,jj5,uu5,nRho,nRho);
            
            XYSize = obj.trueSize(1)*obj.trueSize(2);
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),1:obj.trueSize(1),1:obj.trueSize(3)-1);
            ind6 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii6 = ii1(ind6); jj6 = ii6+XYSize; uu6 = DijkRhoijkPlus(:); uu6 = uu6(ind6);
            A6 = sparse(ii6,jj6,uu6,nRho,nRho);
            
            [X,Y,Z] = meshgrid(1:obj.trueSize(2),1:obj.trueSize(1),2:obj.trueSize(3));
            ind7 = sub2ind(obj.trueSize,Y(:),X(:),Z(:));           
            ii7 = ii1(ind7); jj7 = ii7-XYSize; uu7 = DijkRhoijkMinus(:); uu7 = uu7(ind7);
            A7 = sparse(ii7,jj7,uu7,nRho,nRho);
            
            diffusionJacobian = A1 + A2 + A3 + A4 + A5 + A6 + A7;  
        end
       
    end

    methods (Static)
        
        function checkDimension(dimension)

            arguments
                dimension (1,1) double
            end

            if ~(dimension == 3)
                error("Dimension of data should be 3")
            end

        end  
        
        function matrixOut = doXPlusDifference(matrixIn,BoundaryConditionType,h)
            [h1,h2,h3] = size(matrixIn);
            h = reshape(h,[1,h2-1]);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(:,2:end-1,:);
                matrixNPlus = matrixIn(:,3:end,:);
                matrixOut = (matrixNPlus-matrixN).*repmat(1./h(2:end),[h1,1,h3]);
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(:,1:end-1,:);
                matrixNPlus = matrixIn(:,2:end,:);
                matrixOut = (matrixNPlus-matrixN).*repmat(1./h,[h1,1,h3]);
                matrixOut(:,h2,:) = 0;
            end
        end
        
        function matrixOut = doXMinusDifference(matrixIn,BoundaryConditionType,h)
            [h1,h2,h3] = size(matrixIn);
            h = reshape(h,[1,h2-1]);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(:,2:end-1,:);
                matrixNMinus = matrixIn(:,1:end-2,:);
                matrixOut = (matrixN - matrixNMinus).*repmat(1./h(1:end-1),[h1,1,h3]);
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(:,2:end,:);
                matrixNMinus = matrixIn(:,1:end-1,:);
                matrixOut(:,2:h2,:) = (matrixN - matrixNMinus).*repmat(1./h,[h1,1,h3]);
                matrixOut(:,1,:) = 0;
            end
        end
        
        function matrixOut = doYPlusDifference(matrixIn,BoundaryConditionType,h)
            [h1,h2,h3] = size(matrixIn);
            h = reshape(h,[h1-1,1]);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(2:end-1,:,:);
                matrixNPlus = matrixIn(3:end,:,:);
                matrixOut = (matrixNPlus-matrixN).*repmat(1./h(2:end),[1,h2,h3]);
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(1:end-1,:,:);
                matrixNPlus = matrixIn(2:end,:,:);
                matrixOut = (matrixNPlus-matrixN).*repmat(1./h,[1,h2,h3]);
                matrixOut(h1,:,:) = 0;
            end
        end
        
        function matrixOut = doYMinusDifference(matrixIn,BoundaryConditionType,h)
            [h1,h2,h3] = size(matrixIn);
            h = reshape(h,[h1-1,1]);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(2:end-1,:,:);
                matrixNMinus = matrixIn(1:end-2,:,:);
                matrixOut = (matrixN - matrixNMinus).*repmat(1./h(1:end-1),[1,h2,h3]);
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(2:end,:,:);
                matrixNMinus = matrixIn(1:end-1,:,:);
                matrixOut(2:h1,:,:) = (matrixN - matrixNMinus).*repmat(1./h,[1,h2,h3]);
                matrixOut(1,:,:) = 0;
            end
        end
        
        function matrixOut = doZPlusDifference(matrixIn,BoundaryConditionType,h)
            [h1,h2,h3] = size(matrixIn);
            h = reshape(h,[1,h3-1]);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(:,:,2:end-1);
                matrixNPlus = matrixIn(:,:,3:end);
                matrixOut = (matrixNPlus-matrixN).*permute(repmat(1./h(2:end),[h1,1,h2]),[1,3,2]);
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(:,:,1:end-1);
                matrixNPlus = matrixIn(:,:,2:end);
                matrixOut = (matrixNPlus-matrixN).*permute(repmat(1./h,[h1,1,h2]),[1,3,2]);
                matrixOut(:,:,h3) = 0;
            end
        end
        
        function matrixOut = doZMinusDifference(matrixIn,BoundaryConditionType,h)
            [h1,h2,h3] = size(matrixIn);
            h = reshape(h,[1,h3-1]);
            if strcmp(BoundaryConditionType,'ccd')
                matrixN = matrixIn(:,:,2:end-1);
                matrixNMinus = matrixIn(:,:,1:end-2);
                matrixOut = (matrixN - matrixNMinus).*permute(repmat(1./h(1:end-1),[h1,1,h2]),[1,3,2]);
            end
            if strcmp(BoundaryConditionType,'ccn')
                matrixN = matrixIn(:,:,2:end);
                matrixNMinus = matrixIn(:,:,1:end-1);
                matrixOut(:,:,2:h3) = (matrixN - matrixNMinus).*permute(repmat(1./h,[h1,1,h2]),[1,3,2]);
                matrixOut(:,:,1) = 0;
            end
        end
        
        function rhoGhost = addGhostPoints(rho)
            rhoGhost = zeros(size(rho)+2);
            rhoGhost(2:end-1,2:end-1,2:end-1) = rho;
            rhoGhost(1,:,:) = rhoGhost(2,:,:);
            rhoGhost(end,:,:) = rhoGhost(end-1,:,:);
            rhoGhost(:,1,:) = rhoGhost(:,2,:);
            rhoGhost(:,end,:) = rhoGhost(:,end-1,:);
            rhoGhost(:,:,1) = rhoGhost(:,:,2);
            rhoGhost(:,:,end) = rhoGhost(:,:,end-1);
        end
        
    end
end
