% NEPv for sum of trace ratio optimization (7.1): 
%
% 	max (1-alpha) * [tr(X'AX)/tr(X'BX)] + alpha *[tr(X'D)/sqrt(tr(X'BX))]
% 	s.t. X'X = I.
%
% Test level-shifted SCF on 3-by-3 examples with k = 1 and k = 2.
%

clear; close all;
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultTextInterpreter','latex');
rng(1);

ex = 1; % ex = 1, single vector D; 2, two vectors D

% ------------------------------------------------------------------------
% 1. Generate testing matrices and parameters
% ------------------------------------------------------------------------
A =[ 
  -3.242  -0.450   1.807
  -0.450  -1.630   0.790
   1.807   0.790   0.226];
B =[
   0.592   1.873   0.175
   1.873   6.332   0.617
   0.175   0.617   0.488];

if ex == 1
	D =[ 
	  -9.122
	   0.421
	   3.134];

	n = 3; k = 1;
	rrb = 800;
	alpha = .6;

elseif ex == 2
	D =[ 
	  -1.430   2.768
	  -0.120  -0.630
	   1.098   2.229];

	n = 3; k = 2;
	rrb = 100;
	alpha = .5;
end

maxit = 1000; mtol = 1.0E-13;
sigma = 100; % can try other level-shifts


% ------------------------------------------------------------------------
% 2. Generate NEPv and solve with SCF
% ------------------------------------------------------------------------
[phi, psi, Hphi, Hpsi, gradHphi, gradHpsi] = BuildSumTrRatio(A, B, D, alpha);
[H, G] = GenGH(D, phi, psi, Hphi, Hpsi, gradHphi, gradHpsi);
   
[V0,D0] = eig(A,B); % initialize V0 by largest eigenvectors of A - lam*B
[~,idx] = sort(real(diag(D0)), 'descend'); 
V0 = orth( V0(:,idx(1:k)) );
[Vs, Res, VV] = RunSCF(G, V0, sigma, maxit, mtol);


% ------------------------------------------------------------------------
% 3.  Plot spectral radius rho(L_sigma) as a function of level-shift sigma
% ------------------------------------------------------------------------
SIGMA2 = linspace(0, rrb, 3000);

[GradGX, GradGX0]  = GenGradG(D, D, eye(k), Vs, psi, Hpsi, gradHphi, gradHpsi); % gradients

[spLs2, s_theoretic] = EstimSprdLS(G, GradGX, Vs, SIGMA2); % estimate spectral radius

figure(2);
plot(SIGMA2, spLs2, '-', 'linewidth', 2); hold on;
xlabel('$\sigma$')
ylabel('spectral radius')
plot([-5, SIGMA2(end)], [1,1], '--k')
plot([s_theoretic, s_theoretic], [0, max(spLs2)*1.05], '--k')
xlim([-5, SIGMA2(end)])
ylim([min(spLs2)*0.7, max(spLs2)*1.05])


% ------------------------------------------------------------------------
% 4.  Plot convergence of SCF at selected level-shifts 
% ------------------------------------------------------------------------
if ex == 1 % selected level-shift
	SIGMA = sort([10.5, 41.5, 200, 450, 800, s_theoretic]); % k=1
	xfst1 =[0,-17,0,0,0,0]; % manually set label offset for better display
	xfst2 =[15,6,15,-15,-15,-50]; 
	yfst2 =[-.05,-.02,-.05,-.1,-.1,-.1];

elseif ex == 2
	SIGMA = sort([1.9, 4.5, 20, 50, 100, s_theoretic]); % k=2
	xfst1 =[0,0,0,0,0,0]; % manually set label offset for better display
	xfst2 =[2,2,2,-2,-2,-7]; 
	yfst2 =[0,0,0,-.1,-.1,-.1];
end

% 4.1 Estimate spectral radius over SIGMA
[spLs, ~] = EstimSprdLS(G, GradGX, Vs, SIGMA); % 

% 4.2 Run SCF and plot observed convergence rate
sigmaLs = []; convLs = [];
V0 = orth(Vs+randn(n,k)*1.0E-1); % initial vectors
 
for ii = 1:length(SIGMA)
	if spLs(ii) < 1 % only show converged case
    	s = SIGMA(ii);
        [Vsls, Resls, VVls] = RunSCF(G, V0, s, maxit, mtol);
        convLs = EstimObserved(Resls, 5, 1);

		figure(1); %convergence history
        semilogy(Resls, '-', 'linewidth', 2); hold on;
        xvar = min([180, max(find(Resls>mtol*10))]);
		yvar = Resls(xvar)*2;
        %strr = ['$\sigma = $', num2str(s,'%4.2f')];
       	strr = ['$\sigma_',num2str(ii), '$'];
		text(xvar+xfst1(ii), yvar, strr, 'fontsize',15); % text label for sigma
		
		figure(2);
   		plot(s, convLs, 'or', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'r');
		text(s+xfst2(ii), convLs+yfst2(ii), strr, 'Fontsize', 15);
    end
end
figure(1);
xlim([0,250]);
ylim([1.0E-11, 1.0E2]);
ylim([1.0E-13, 1.0E0]);
xlabel('SCF iteration $i$'); ylabel('$\mbox{NRes}(X_i)$')

% ------------------------------------------------------------------------
% 5. MISC: display optimal level-shift
% ------------------------------------------------------------------------
[spls_opt, idx] = min(spLs2);
sigma_opt = SIGMA2(idx);
[~, idx] = find(spLs2<1, 1,'first');
sigma_cut = SIGMA2(idx);

disp('sigma_opt, rho(opt)')
[sigma_opt, spls_opt] 
disp('sigma_cut, sigma_lb')
[sigma_cut, s_theoretic]

%
return;
% END
