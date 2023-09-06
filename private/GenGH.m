function [H, G] = GenGH(D, phi, psi, Hphi, Hpsi, gradHphi, gradHpsi)
% Function GenGH generates coefficient matrix functions:
%	H(X) 
%	G(X) is the aligned H([X])

tdx = @(X) trace(X'*D);
H = @(X) Hphi(X) + tdx(X)*Hpsi(X)+psi(X)*(D*X'+X*D');
alignx = @(X) X*GetQ(X'*D);
G = @(X) H( alignx(X) );

return;
% END OF GenGH
