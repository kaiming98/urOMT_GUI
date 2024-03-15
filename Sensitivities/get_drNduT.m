function [drNduT] = get_drNduT(M,nt,dt,par,y)
%% Sensitivity of rho(:,nt) w.r.t 'u' transpose times a vector x
%  y... vector length(prod(n),1)

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)
%               S1                      S2

%% 
n                  = par.n;
sensTx             = zeros(par.dimension*prod(n),nt);
sens               = y;
switch par.diffusionCoefficientType
    case 'constant'
        for i = nt:-1:1
            % Sensitivity:          
            [sensI,pcgflag1]   = pcg(par.B',sens);
            if pcgflag1 ~= 0
                warning('MATLAB:pcgExitFlag','Warning: get_drNduT4.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
            end
            if par.dimension==2
                y2 = [M.Tx{i},M.Ty{i}]'*(dt*sensI);
            elseif par.dimension==3
                y2 = [M.Tx{i},M.Ty{i},M.Tz{i}]'*(dt*sensI);
            end
            sensTx(:,i)              = y2;

            if  i>1
                if ~par.addSource
                    sens = M.S{i}'*sensI;
                else
                    sens = M.R{i}.*(M.S{i}'*sensI);
                end
            end
        end
    case {'anisotropic', 'autoAnisotropic', 'brightness'}
        for i = nt:-1:1
            % Sensitivity:    
            sensI = sens;
            if par.dimension==2
                y2 = [M.Tx{i},M.Ty{i}]'*(dt*sensI);
            elseif par.dimension==3
                y2 = [M.Tx{i},M.Ty{i},M.Tz{i}]'*(dt*sensI);
            end
            sensTx(:,i)              = y2;

            if  i>1 
                if ~par.addSource
                    sens = (M.S{i}+dt*M.dD{i})'*sensI;
                else
                    sens = dt*M.dD{i}'*sensI + M.R{i}.*(M.S{i}'*sensI);
                end
            end
        end
end

drNduT = sensTx(:);
end
