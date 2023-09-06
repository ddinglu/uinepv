function [sps, lb] = EstimSprdLS(G, GradGX, Xs, sigmas)
% function EstimSprdLS estimates the spectral radius of the linear operator 
% corresponding the level-shifted NEPv with sigma at Xs 
% INPUT:
%   G(X)        - coefficient matrix of the aligned NEPv;
%   GradGX(X,E) - directional derivative of G at X in the direction of E;
%   Xs          - target solution of the NEPv.
%   sigmas		- a vector of level-shifts to evaluate spectral radius 
%
% OUTPUT:
%	sps			- spectral radius at sigmas
%	lb			- estimated lower bound for convergent level-shift

n = size(Xs,1);
k = size(Xs,2);

% compute eigenvalue decomposition
Gx = G(Xs);  % smallest to largest eigenvalues
[Vx, Dx] = eig(Gx); 
[~,idx] = sort(real(diag(Dx)), 'descend');
Vx = Vx(:,idx); Dx = Dx(idx,idx);

% construct linear operators
lam = diag(Dx);
lam1 = lam(1:k); 
lam2 = lam(k+1:end);
Ss = @(s) -1./(lam2 - lam1'-s);
Vs = Vx(:,1:k); Vst = Vx(:,k+1:end);
Op = @(X,s) Ss(s).*( Vst'*GradGX(Vs, Vst*X)*Vs ) + s*Ss(s).*X; 

% compute spectral radius
sps = zeros(size(sigmas)); 
for ii = 1:length(sigmas)
    s = sigmas(ii);
    Afun = @(x) reshape(Op(reshape(x, n-k, k), s),[],1);
    lammax = eigs(Afun, (n-k)*k, 1);
    sps(ii) = abs(lammax);
end


% estimate theoretical lower bound for a convergent level-shift
% - estimate spectral range of the linear operator Qs 
Opqs = @(X) -Vst'*GradGX(Vs, Vst*X)*Vs - diag(lam2)*X + X*diag(lam1);
Afun = @(x) reshape(Opqs(reshape(x, n-k, k)),[],1);
mux = eigs(Afun, (n-k)*k, 2, 'bothendsreal');
mu_max = max(mux);
mu_min = min(mux);
deltas = -lam(k+1) + lam(k);
lb = mu_max/2-deltas; % theoretical lower bound

% MAY BE ADDED IN THE FUTURE: the estimate optimal level-shift
%spans = -lam(end) + lam(1);
%rhobnd = @(s) max(abs([mu_max./(s+deltas)-1;  mu_min./(s+spans)-1]));
%fs = @(s) mu_max./(s+deltas) + mu_min./(s+spans)-2;
%s_opt = fzero(fs, [0,3000]); % bnd optimal 

return 
% END OF EstimSprd
