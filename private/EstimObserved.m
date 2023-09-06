function factor = EstimObserved(Res, nn, emode)
% function EstimObserved estimates the observed convergence rate for the
% residual sequence Res, using linear interpolation with the last nn
% samples (default nn=10).
% mode: 	0 for standard case;  1 for possible `later' convergence;

if nargin < 2, nn = 10; end 
if nargin < 3, emode = 0; end % 

if emode == 0
   	ne = find(Res< max(Res(end)*2,1.0E-11), 1); % index reaching tol~1.0E-11
else
   	ne = find(Res< Res(end)*10, 1); 			% index reaching tol~mtol*10
end

ns = find(Res< Res(end)*1.0E8, 1);  % index reaching tol*1.0E8
ns = max([ns, ne-nn]);		% use at most nn samples
dd = polyfit([ns:ne],log(Res(ns:ne)),1);
factor = exp(dd(1));

return 
% END OF EstimObserved


