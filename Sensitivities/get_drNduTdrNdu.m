function [drNduTdrNdu] = get_drNduTdrNdu(M,nt,dt,par,x)
%% Created by Xinan Chen on 04/19/21
%% This function combines get_drNduT.m and get_drNdu.m to compute the second term in Hessian
%% Sensitivity of rho(:,end) transpose w.r.t 'u' * rho(:,end) w.r.t 'u' * full vector

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)
%%
n                  = par.n;
X                  = reshape(x,prod(n)*par.dimension,nt);
sensx              = zeros(prod(n),nt);

switch par.diffusionCoefficientType
    case 'constant'
        % first part:  rho(:,end) w.r.t 'u' * full vector
        for i = 1:nt
            % Sensitivity:
            if  i>1
                for j = 1:i-1
                    if ~par.addSource
                        [sensx(:,j),pcgflag2] = pcg(par.B,M.S{i}*sensx(:,j));
                    else
                        [sensx(:,j),pcgflag2] = pcg(par.B,M.S{i}*(M.R{i}.*sensx(:,j)));
                    end
                    if pcgflag2 ~= 0
                        warning('MATLAB:pcgExitFlag','Warning: get_drNduTdrNdu.m >>> while finding drho(:,n)/du_%d, pcg exit pcgflag2 = %d',j,pcgflag2)
                    end
                end
            end
            if par.dimension==2
                y1 = [M.Tx{i},M.Ty{i}]*(dt*X(:,i));
            elseif par.dimension==3
                y1 = [M.Tx{i},M.Ty{i},M.Tz{i}]*(dt*X(:,i));    
            end
            [sensx(:,i),pcgflag3] = pcg(par.B,y1);
            if pcgflag3 ~= 0
                warning('MATLAB:pcgExitFlag','Warning: get_drNduTdrNdu3.m >>> while finding drho(:,%d)/du, pcg exit pcgflag3 = %d',i,pcgflag3)
            end
        end

        sens = sum(sensx,2);


        % second part: rho(:,end) transpose w.r.t 'u' * full vector

        sensTx             = zeros(par.dimension*prod(n),nt);


        for i = nt:-1:1
            % Sensitivity:          
            [sensI,pcgflag1]   = pcg(par.B',sens);
            if pcgflag1 ~= 0
                warning('MATLAB:pcgExitFlag','Warning: get_drNduTdrNdu.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
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
        % first part:  rho(:,end) w.r.t 'u' * full vector
        for i = 1:nt
            % Sensitivity:
            if  i>1
                if ~par.addSource
                    for j = 1:i-1
                        sensx(:,j) = (M.S{i}+dt*M.dD{i})*sensx(:,j);
                    end
                else
                    for j = 1:i-1
                        sensx(:,j) = dt*M.dD{i}*sensx(:,j) + M.S{i}*(M.R{i}.*sensx(:,j));
                    end
                end
            end
            if par.dimension==2
                y1 = [M.Tx{i},M.Ty{i}]*(dt*X(:,i));
            elseif par.dimension==3
                y1 = [M.Tx{i},M.Ty{i},M.Tz{i}]*(dt*X(:,i));    
            end
            sensx(:,i) = y1;
        end

        sens = sum(sensx,2);


        % second part: rho(:,end) transpose w.r.t 'u' * full vector

        sensTx             = zeros(par.dimension*prod(n),nt);


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

drNduTdrNdu = sensTx(:);
end

