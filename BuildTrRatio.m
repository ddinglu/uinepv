function [phi, psi, Hphi, Hpsi, gradHphi, gradHpsi] = BuildTraceRatio(A, B, D, theta)
% function BuildTrRatio generates coefficient functions and the derivatives 
% of NEPv for the trace-ratio problem, i.e., (7.2):
%
% 	max [tr(X'AX+X'D)]/[tr(X'BX)^theta]
% 	s.t. X'X = I.
%
% Output function handles:
%
%	phi(X), psi(X), Hphi(X), Hpsi(X), gradHphi(X,E), gradHpsi(X,E)
%
rb = @(X) trace(X'*B*X);
ra = @(X) trace(X'*A*X);
tdx = @(X) trace(X'*D);

if theta == 0 
	% psi is a constant 1 
	phi = @(X) ra(X);
	psi = @(X) 1;

	Hphi = @(X) 2*A;
	Hpsi = @(X) 0*B;

	gradHphi = @(X,Y) 0*B;
	gradHpsi = @(X,Y) 0*B;

else
	srb = @(X)(rb(X))^theta;

	phi = @(X) ra(X)/srb(X);
	psi = @(X) 1/srb(X);

	Hphi = @(X) 2*psi(X)*(A - theta*ra(X)/rb(X)*B);
	Hpsi = @(X) -2*theta*psi(X)/rb(X)*B;


	Hphi1 = @(X) 2*psi(X)*(A - ra(X)/rb(X)*B); % Hphi with theta = 1
	gradHphi = @(X,Y) -2*theta*trace(Y'*B*X)/rb(X) * Hphi1(X) - 2*theta*trace(Y'*Hphi(X)*X)/rb(X)*B;
	gradHpsi = @(X,Y) -2*(theta+1)*trace(Y'*Hpsi(X)*X)/rb(X) *B;
end

end % END OF BuildTrRatio




