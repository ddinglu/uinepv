% Sum of trace ratio optimization (7.1): 
%
% 	max (1-alpha) * [tr(X'AX)/tr(X'BX)] + alpha *[tr(X'D)/sqrt(tr(X'BX))]
% 	s.t. X'X = I.
%
% Test D with full and reduced column.
%

close all; clear all;
set(0,'defaultTextInterpreter','latex');
rng(0); 

% -------------------------------------------------------------------------
% 0. Uncomment to run `full' / 'deficient' rank test
% -------------------------------------------------------------------------

ex = 'full'; exlabel = 'k'; % full rank test
%ex = 'deficient'; exlabel = 'r_D'; % deficient rank test

% -------------------------------------------------------------------------
% 1. Testing matrices and parameters
% -------------------------------------------------------------------------
maxit = 100000; mtol = 1.0E-13;

% -- set size of random testing problem: k < n 
% n = 200, 20 samples, timing estimate, 10m; (as used in the paper)
% n = 400, 20 samples, timing estimate, 20m
% n = 800, 20 samples, timing estimate, 1h

n = 60; k = 40; 

% -- set number of sampled alpha to evaluate spectral radius and observed rate
nsample = 200; 
alphan = flip(linspace(0,1,nsample+1)); 	% for spectral radius plot
sample_alpha = flip(linspace(0,1,10+1)); 	% for observed rate 

% -- set coefficient matrices and parameters
A = diag(2*ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
B = diag([1:n]);
D00 = randn(n,k+10); 
corder = colororder; 
ind = 1;

RKS = 10:10:k; % testing rank ell = 10, 20, ..., k

% ------------------------------------------------------------------------
% 2. Main loop over the rank rk in RKS
% ------------------------------------------------------------------------
for rk = RKS 
	if strcmp(ex, 'full') % full rank D
    	D1 = D00(:,1:rk);
		P = eye(rk); 
		D = D1;
		nk = rk; % size D
	else
    	D1 = D00(:,1:rk);
		nk = k+10; % size D
		P = eye(nk,rk); %P = orth(randn(nk,rk)); 
    	D = D1*P'; 
	end

	% ---------------------------------------------------------------------
    % 2.1 Solve NEPv by SCF and check the rate of convergence
	% ---------------------------------------------------------------------
	disp('--- ') 
	disp('Check observed rate:') 
	[sample_conv, sample_sprd, res_hist] = SolveNEPv_mod(A, B, D, D1, P, 1, sample_alpha, 0, mtol);
	disp('--- ') 
	disp('Generate spectral radius:')
	[conv_hist, sprd_hist] = SolveNEPv(A, B, D, D1, P, 1, alphan, 0, mtol);

    SPRD_rk{ind} = sprd_hist; % spectral radius
    Conv_rk{ind} = conv_hist; % observed rate of convergence 

    SPSPRD_rk{ind} = sample_sprd;
    SPConv_rk{ind} = sample_conv;

	% ---------------------------------------------------------------------
    % 2.2 Plot curves of spectral radius and observed rate of convergence
	% ---------------------------------------------------------------------
    figure(1);
    h1{ind} = plot(alphan, sprd_hist, '-', 'LineWidth', 2, 'DisplayName',['$',exlabel,'=$ ',num2str(rk)], 'color', corder(ind,:)); hold on;

    figure(2); 
    h2{ind} = plot(alphan, sprd_hist, '-', 'LineWidth', 2, 'DisplayName',['$',exlabel,'=$ ',num2str(rk)], 'color', corder(ind,:)); hold on;
    h3{ind} = plot(sample_alpha, sample_conv, 'or','MarkerSize',6,'LineWidth',1, 'color', corder(ind,:)); hold on;
	ind = ind + 1;

end % loop over rank

figure(1); % spectral radius 
xlabel('$\alpha$');
ylabel('spectral radius');
ylim([0,1.0])
legend("show",'Location','southeast','Interpreter','latex','FontSize',12);

figure(2); % rate of convergence
xlabel('$\alpha$');
ylabel('convergence rate');
ylim([0.6,1.0])

% ------------------------------------------------------------------------
% 3. MISC: Save computation result for large simulation
% ------------------------------------------------------------------------
if n >= 200
	save(['nepv7_1_',ex,'_n',num2str(n),'.mat']);
end

return; 
