function [observe_hist, sprd_hist, res_hist] = SolveNEPv_mod(A, B, D, D1, P, mode, param, sigma, mtol, reusev0, V0, maxit)
% function SolveNEPv_mod runs SCF and returns the observed and estimated
% convergence factor of SCF.
% 
% In SolveNEPv_mod: repeated initial guess is tried until good
% estimation of observed rates are obtained.
%
% Input: 
%	A, B, D, D1, P -- coefficient matrices
%	mode	-- mode = 1 for NEPv (7.1), and 2 for NEPv (7.2)
%	para	-- vector of parameters alpha for NEPv (7.1) and theta for NEPv (7.2)
%	sigma	-- level shifts
%	tol 	-- convergence tolerance
%	V0 		-- initial V0
% 	reusev0 -- use fixed V0 if reusev0==true; may be slow
%
% Output: 
%	observe_hist 	-- observed convergence rate
%	sprd_hist 		-- spectral radius
%	res_hist 		-- residual history of SCF


n = size(D,1);
k = size(D,2);

if nargin<12
	maxit = 100000;
end

if nargin < 11
	[V0,E0] = eig(A,B);
	[~,idx1] = sort(real(diag(E0)), 'descend'); 
	V0 = orth( V0(:,idx1(1:k)) );
end
Vint = V0;

if nargin<10
	reusev0 = false;
end

observe_hist = [];     	% observed rates 
sprd_hist = [];     	% spectral radius

for jj = 1:length(param)

	pj = param(jj);

    %   1. construct the problem
	if mode == 1 % sum trace ratio problem
		[phi, psi, Hphi, Hpsi, gradHphi, gradHpsi] = BuildSumTrRatio(A, B, D, pj);
	else % trace ratio problem
		[phi, psi, Hphi, Hpsi, gradHphi, gradHpsi] = BuildTrRatio(A, B, D, pj);
	end

    %   2. construct the coefficient matrices H(x), G(X) 
	
	[H, G] = GenGH(D, phi, psi, Hphi, Hpsi, gradHphi, gradHpsi);

    %   3. plain SCF iteration: try 20 random V0 until good estimation of observed rates obtained.
	if reusev0, V0 = Vint; end
	V00 = V0;
	cdiff = inf;
	for ii = 1:20 
		epsln = 1/n;
	   	V0 = orth(V00+randn(n,k)*epsln);
		[V0, Res, VV] = RunSCF(G, V0, -1*sigma, maxit, mtol/n);

	    %  4. observed convergence rate
		Conv= EstimObserved(Res,10,1);  	% observed convergence rate by plain SCF 

		if ii == 1 % Get reference spectral radius (no need to recalculate)
			%  5. spectral radius 
			% -- construct gradient of G 
			[GradGX, GradGX0]  = GenGradG(D, D1, P, V0, psi, Hpsi, gradHphi, gradHpsi);

			% -- estimate corresponding spectral radius
			%Convest_sp = EstimSprd(G, GradGX, V0); % less efficient 
			Convest_sp = EstimSprds(G, GradGX0, V0); % more efficient
		end
	
		if abs(Convest_sp-Conv)<1.0E-3
			mc = Conv;
			mres = Res;
			break;

		elseif abs(Convest_sp-Conv)<cdiff
			cdiff = abs(Convest_sp-Conv);
			mc = Conv;
			mres = Res;
		end
	end

 	res_hist{jj} = mres;
	observe_hist = [observe_hist,mc];
    sprd_hist = [sprd_hist,Convest_sp];
	disp(['Observed = ', num2str(mc), ';  Spectral radius = ', num2str(Convest_sp)]);
end
% END OF SolveNEPv_mod
