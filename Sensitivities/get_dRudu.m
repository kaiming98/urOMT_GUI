function dRdu=get_dRudu(u,nt,par)
% -----------------------------
% Created by:       Rena Elkin
% Last modified on: 04/13/2017
% -----------------------------
%This function returns the derivative of the regularization term
%R(u)=||grad(u)||^2

G=par.Grad'*par.Grad;%G=-par.Grad'*par.Grad;
U=vec2mat(u(:),par.dimension*nt);
dRdu=G*U;
dRdu=dRdu(:)';
end

