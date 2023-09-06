function [DG, DGX] = GenGradG(D, D1, P, Xs, psi, Hpsi, gradHphi, gradHpsi)
% function GenGradG constructs the directional derivative of G DG(X,E) 
% using the long assembling formula from the paper.
% Recall: D=D1*P' is a factorization of D1. 
% 
% The second output DGX is a function handle DGX(E) = DG(Xs, E), i.e.,g
% DG evaluated at X=Xs.

% part I: DG
M1 = @(X) GetA(X'*D1); % Tested: M1 = sqrt( (X'D1) * (X'D1)' )
gradM1 = @(X,E) sylvester( M1(X), M1(X), D1'*X*E'*D1 + D1'*E*X'*D1);
tra = @(X) trace(M1(X)); % Tested
Dtra = @(X,E) trace(gradM1(X,E));
Qx = @(X) GetQ(X'*D); 	
Q1x = @(X) GetQ(X'*D1); 
gradQx = @(X,E) ( E'*D1/M1(X) - Q1x(X)*gradM1(X,E)/M1(X) ) * P';

%
DG = @(X,E) gradHphi(X,E) + Dtra(X,E) * Hpsi(X) + tra(X)*gradHpsi(X,E) ...
   	+ trace(E'*Hpsi(X)*X) * (X*Qx(X)*D' + D*Qx(X)'*X' ) ...
	+ psi(X)*(D*gradQx(X,E)'*X' + X*gradQx(X,E)*D'+ E*Qx(X)*D' + D*Qx(X)'*E');


% part II: DGX
[Qxs, Mxs] = GetPolar(Xs'*D);    
[Q1s, M1s] = GetPolar(Xs'*D1);   
gradM1s = @(E) sylvester(M1s, M1s, D1'*Xs*E'*D1 + D1'*E*Xs'*D1);
tras = trace(M1s); 
Dtras = @(E) trace(gradM1s(E));
gradQxs = @(E) ( E'*D1/M1s - Q1s * gradM1s(E)/M1s ) * P';

DGX = @(E) gradHphi(Xs,E) + Dtras(E) * Hpsi(Xs) + tras*gradHpsi(Xs,E) ...
   	+ trace(E'*Hpsi(Xs)*Xs) * (Xs*Qxs*D' + D*Qxs'*Xs' ) ...
	+ psi(Xs)*(D*gradQxs(E)'*Xs' + Xs*gradQxs(E)*D'+ E*Qxs*D' + D*Qxs'*E');

return;
end


% auxiliary functions for polar factors Q and A

function [Q] = GetQ(M)
[U,S,V] = svd(M,0);
Q = U*V';
end

function [A] = GetA(M)
[U,S,V] = svd(M,0);
r = rank(S); % can be improved
A = V(:,1:r)*S(1:r,1:r)*V(:,1:r)';
end


% -----------------------------------------------------------------------
% This part is an older version where Xs is fixed
% -----------------------------------------------------------------------
%
%[Qx, Mx] = GetPolar(Xs'*D);    
%[Q1, M1] = GetPolar(Xs'*D1);   % polar for the condensed if D neq D1. Note M1 = Mx
%
%gradM1 = @(E) sylvester(M1, M1, D1'*Xs*E'*D1 + D1'*E*Xs'*D1);
%tra = trace(M1); 
%Dtra = @(E) trace(gradM1(E));
%
%gradQx = @(E) ( E'*D1/M1 - Q1 * gradM1(E)/M1 ) * P';
%
%%
%DG = @(E) gradHphi(Xs,E) + Dtra(E) * Hpsi(Xs) + tra*gradHpsi(Xs,E) ...
%   	+ trace(E'*Hpsi(Xs)*Xs) * (Xs*Qx*D' + D*Qx'*Xs' ) ...
%	+ psi(Xs)*(D*gradQx(E)'*Xs' + Xs*gradQx(E)*D'+ E*Qx*D' + D*Qx'*E');
%
%return;

