function [phi, psi, Hphi, Hpsi, gradHphi, gradHpsi] = BuildSumTrRatio(A, B, D, alpha)
% function BuildSumTrRatio generates coefficient functions and the
% derivatives of NEPv for the sum of trace-ratio problem, i.e., (7.1): 
%
% 	max (1-alpha) * [tr(X'AX)/tr(X'BX)] + alpha *[tr(X'D)/sqrt(tr(X'BX))]
% 	s.t. X'X = I.
%
% Output function handles:
%
%	phi(X), psi(X), Hphi(X), Hpsi(X), gradHphi(X,E), gradHpsi(X,E)
%

ra = @(X) trace(X'*A*X);
rb = @(X) trace(X'*B*X);
srb = @(X) sqrt(rb(X));
tdx = @(X) trace(X'*D);

phi = @(X) (1-alpha)*ra(X)/rb(X);
psi = @(X) alpha/srb(X);

Hphi = @(X) 2/rb(X)*( (1-alpha)*A - phi(X)*B);
Hpsi = @(X) -psi(X)/rb(X)*B;

gradHphi = @(X,Y) -2*trace(Y'*B*X)/rb(X) * Hphi(X) - 2*trace(Y'*Hphi(X)*X)/rb(X)*B;
gradHpsi = @(X,Y) -3*trace(Y'*Hpsi(X)*X)/rb(X) *B;

end % END OF BuildSumTrRatio
