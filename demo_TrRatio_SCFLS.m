% NEPv 7.2, level-shifted SCF 
% 3-by-3 example, k = 1 and k = 2
% phi(X) = tr(X^TAX)/tr(X^TBX)^theta
% psi(X) = 1/tr(X^TBX)^theta


clear; close all;
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultTextInterpreter','latex');
rng(0);

ex = 2; % ex = 1, single vector D; 2, two vectors D

% ------------------------------------------------------------------------------- 
% 1. Generate testing matrices and parameters
% ------------------------------------------------------------------------------- 

if ex == 1
	A =[ 
	  -3.242  -0.450   1.807
	  -0.450  -1.630   0.790
	   1.807   0.790   0.226];
	B =[
	   0.592   1.873   0.175
	   1.873   6.332   0.617
	   0.175   0.617   0.488];
    D =[ 
        -9.122
        0.421
        3.134];
    n = 3; k = 1;
	D1 = D; 
	P = eye(k);

	theta = 0;
	sigma = 0;

	lsigma = -10; rsigma = 10; % range of level-shift to plot spectral radius
	SIGMAPLT = [-9.82, -9.4, -8.9, 0, 10]; % selected shift for convergence plot

elseif ex == 2

	A =[ 
	   1.145  -0.095   0.514
      -0.095   0.838   1.022
       0.514   1.022  -1.223];
	B =[
	   0.582  -0.037   0.025
      -0.037   0.183   0.043
       0.025   0.043   0.239];
	D =[ 
	   0.760   0.258
   	   0.011   0.774
       0.180   0.520];
	n = 3; k=2;
	D1 = D;
	P = eye(k);

	theta = 3.0;
	sigma = 100;

	lsigma = 0; rsigma = 400; % range of level-shift to plot spectral radius
	SIGMAPLT = [2, 18, 120, 240, 360]; % selected shift for convergence plot
end

maxit = 1000; mtol = 1.0E-13;


% ------------------------------------------------------------------------------- 
% 2. Generate NEPv and solve with SCF
% ------------------------------------------------------------------------------- 
[phi, psi, Hphi, Hpsi, gradHphi, gradHpsi] = BuildTrRatio(A, B, D, theta);
[H, G] = GenGH(D, phi, psi, Hphi, Hpsi, gradHphi, gradHpsi);

[Vab,Eab] = eig(A,B);
ee = diag(Eab); [~,idx0] = sort(real(ee), 'descend'); 
V00 = Vab(:,idx0(1:k));
  
[Vs, Res, VV] = RunSCF(G, V00, sigma, maxit, mtol);


% ------------------------------------------------------------------------------- 
% 3.  Plot spectral radius rho(L_sigma) as a function of level-shift sigma
% ------------------------------------------------------------------------------- 
[GradGX, GradGX0]  = GenGradG(D, D, eye(k), Vs, psi, Hpsi, gradHphi, gradHpsi); % gradients
SIGMAS = linspace(lsigma, rsigma, 3000);

[sps, s_theoretic] = EstimSprdLS(G, GradGX, Vs, SIGMAS);

figure(2);
plot(SIGMAS, sps, '-', 'linewidth', 2); hold on;
xlabel('$\sigma$')
ylabel('spectral radius')
plot([lsigma-2, rsigma], [1,1], '--k')
plot([s_theoretic, s_theoretic], [min(sps)-0.1, max(sps)+0.1], '--k')
xlim([lsigma-0.05*abs(rsigma-lsigma), rsigma])


% ------------------------------------------------------------------------------- 
% 4.  Plot convergence of SCF at selected level-shifts 
% ------------------------------------------------------------------------------- 
if ex == 1 % selected level-shift
	SIGMA = sort([SIGMAPLT, s_theoretic]); % k=1
	xfst1 =[1,1,1,1,1,1]; % manually set label offset for better display
	xfst2 =[0.5, 0.3, 0.5, 0.5, 0, -1]; 
	yfst2 =[0.0, 0.04, -0.01, -0.01, -0.1, -0.1];

elseif ex == 2
	SIGMA = sort([SIGMAPLT, s_theoretic]); % k=1
	xfst1 =[0,0,0,0,0,0]; % manually set label offset for better display
	xfst2 =[10,10,10,-0,-0,-0]; 
	yfst2 =[0,0,0,-.05,-.05,-.05];
end

% 4.1 Estimate spectral radius over SIGMA
[spLs, ~] = EstimSprdLS(G, GradGX, Vs, SIGMA); % 

% 4.2 Run SCF and plot observed convergence rate
sigmaLs = []; convLs = [];
 
V0 = orth(Vs+randn(n,k)*1.0E-2); % initial vectors
for ii = 1:length(SIGMA)
	if spLs(ii) < 1 % only show converged case
    	s = SIGMA(ii);
        [Vsls, Resls, VVls] = RunSCF(G, V0, s, maxit, mtol);
        convLs = EstimObserved(Resls, 5, 1);

		figure(1); %convergence history
        semilogy(Resls, '-', 'linewidth', 2); hold on;
		if ex == 1 
			xvar = min([50, max(find(Resls>mtol*10))]);
		elseif ex == 2
			xvar = min([80, max(find(Resls>mtol*10))]);
		end
		yvar = Resls(xvar)*2;
        %strr = ['$\sigma = $', num2str(s,'%4.2f')];
       	strr = ['$\sigma_',num2str(ii), '$'];
		text(xvar+xfst1(ii), yvar, strr, 'fontsize',15); % text label for sigma
		
		figure(2);
   		plot(s, convLs, 'or', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'r');
		text(s+xfst2(ii), convLs+yfst2(ii), strr, 'Fontsize', 15);
    end
end

if ex == 1 
	figure(1);
	xlim([0,60]); 
	ylim([1.0E-13, 1.0E-1])
	figure(2)
	ylim([min(sps) - 0.1, max(sps) + 0.1])
elseif ex == 2
	figure(1);
	xlim([0,100]); 
	ylim([1.0E-13, 1.0E-0])
	figure(2)
	ylim([0.2, 1.2])
end

xlabel('SCF iteration $i$'); ylabel('$\mbox{NRes}(X_i)$')

% ------------------------------------------------------------------------------- 
% 5. MISC: display optimal level-shift
% ------------------------------------------------------------------------------- 
[spls_opt, idx] = min(sps);
sigma_opt = SIGMAS(idx);
[~, idx] = find(sps<1, 1,'first');
sigma_cut = SIGMAS(idx);

disp('sigma_opt, rho(opt)')
[sigma_opt, spls_opt] 
disp('sigma_cut, sigma_lb')
[sigma_cut, s_theoretic]

%
return;
% END
