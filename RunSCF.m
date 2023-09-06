function [Vs, Res, VV] = RunSCF(G, V0, sigma, maxit, mtol)
% scfsub runs level-shifted SCF for NEPv G(V)*V = V*D.
%
% INPUT:
%	G: a function handle for G(V);
%	V0: initial vectors
%	sigma: level-shift; if sigma = 0, then runs plain SCF
%	maxit, mtol: max iteration number and relative residual tolerance
%
% OUTPUT:
%	Vs: approximate solution
%	Res: a vector containing relative residuals of each iteration
%	VV:	VV{i} is the i-th iteration 

n = size(V0,1); 
k = size(V0,2);
Gv = G(V0);
D0 = V0'*(Gv*V0);
Res = norm(Gv*V0-V0*D0,1);
VV = [];
VV{1} = V0;
for i = 1:maxit
    [Vi,Di] = eig(Gv + sigma*V0*V0');
    [~,idx] = sort(real(diag(Di)), 'descend'); 
    V0 = Vi(:,idx(1:k));
    Gv = G(V0);
    D0 = Di(idx(1:k),idx(1:k)) - sigma*eye(k);
    Res =[ Res, norm(Gv*V0-V0*D0,1)/norm(Gv,1)];
    VV{i} = V0;
    if(Res(end)< mtol*n)
        break
    end
end
if i == maxit 
	disp('SCF not converged');
end

Vs = V0;  

end % END OF RunSCF
