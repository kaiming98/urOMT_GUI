function [rho] = advecDiff(rho0, u, obj)
    %% rho(:,1...T) = transport (rho0 - initial mass density,
    %                                 u - velocity vector for nt time steps,
    %                                     size(prod(n),nt)

    % Advection: Semi- Lagrangian  rho^(i,ad) = S(U^(i)*rho^(i-1);

    % Dispersion: Implicitly: (I - dt*Div*D*Grad)*rho^(i) = rho^(i,ad)

    % rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
    %             u(:,1)                   u(:,2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if obj.dimension == 2
        %% 2D version
        n = obj.n;
        rho = zeros(prod(n), obj.nt + 1);
        u = reshape(u, 2 * prod(n), obj.nt);

        if ~isempty(obj.B)
            obj = obj.initializeCellCenteredData;
        else
            B = obj.B;
        end

        rho(:, 1) = rho0;

        for i = 2:obj.nt + 1
            U1 = reshape(u(1:prod(n), i - 1), n');
            U2 = reshape(u(prod(n) + 1:end, i - 1), n');

            S = dTrilinears2d(rho(:, i - 1), obj.Xc + obj.dt * U1, obj.Yc + obj.dt * U2, ...
                obj.h1(1), obj.h2(1), obj.bc);

            rho(:, i) = S * rho(:, i - 1);          

            % diffusion step
            [rho(:, i), pcgflag] = pcg(B, rho(:, i));

            if pcgflag ~= 0
                warning('MATLAB:pcgExitFlag', 'Warning: advecDiff.m >>> while finding rho(:,%d), pcg exit flag = %d', i, pcgflag)
            end
        end

        rho = rho(:, 2:end);

    elseif obj.dimension == 3
        %% 3D version
        %%

        n = obj.n;
        rho = zeros(prod(n), obj.nt + 1);
        u = reshape(u, 3 * prod(n), obj.nt);


        if isempty(obj.B)
            obj = obj.initializeCellCenteredData;
        end
        
        B = obj.B;
        rho(:, 1) = rho0;

        for i = 2:obj.nt + 1
            U1 = reshape(u(1:prod(n), i - 1), n');
            U2 = reshape(u(prod(n) + 1:2 * prod(n), i - 1), n');
            U3 = reshape(u(2 * prod(n) + 1:end, i - 1), n');


            S = dTrilinears3d(rho(:, i - 1), obj.Xc + obj.dt * U1, obj.Yc + obj.dt * U2, obj.Zc + obj.dt * U3, ...
                obj.h1(1), obj.h2(1), obj.h3(1), obj.boundaryConditions);


            % advection step
            rho(:, i) = S * rho(:, i - 1);

            % diffusion step
            [rho(:, i), pcgflag] = pcg(B, rho(:, i));

            if pcgflag ~= 0
                warning('MATLAB:pcgExitFlag', 'Warning: advecDiff.m >>> while finding rho(:,%d), pcg exit flag = %d', i, pcgflag)
            end

        end

        rho = rho(:, 2:end);

    else
        warning('In advecDiff.m: dimension of data should be either 2 or 3')
    end

end
