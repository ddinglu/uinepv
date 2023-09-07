% NEPv for trace-ratio optimization (7.2):
%
% 	max [tr(X'AX+X'D)]/[tr(X'BX)^theta]
% 	s.t. X'X = I.
%
% Test SCF on 3-by-3 examples with k = 1 and k = 2.
%

close all; clear all;
set(0,'defaultTextInterpreter','latex');
rng(0); 

ex = 1; % ex = 1, single vector D; 2, two vectors D

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

    nsample = 201;      % set number of sampled theta
    thetan = linspace(-0.5, 1.5, nsample);
	sigma = 0;

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

	nntt = 19;
	nsample = 201;
	thetan = linspace(0, 6, nsample);
	sigma = -40;
end

corder = colororder;
maxit = 100000; mtol = 1.0E-13;


% ------------------------------------------------------------------------------- 
% 2. Solve NEPv by SCF
% ------------------------------------------------------------------------------- 
V00 = orth(rand(n,k)+randn(n,k)); % common random initial vector for SCF
[conv_hist, sprd_hist] = SolveNEPv(A, B, D, D1, P, 2, thetan, sigma, mtol, 0, V00);


% ------------------------------------------------------------------------------- 
% 3. Estimate rates of convergence at sampled thetas where convergence happens
% ------------------------------------------------------------------------------- 
if ex == 1
    sample_theta = linspace(-0.5, 1.5, 21);
elseif ex == 2
	Ia = find(sprd_hist > 1, 1, 'first'); % find range for non-convergent theta
	Ib = find(sprd_hist > 1, 1, 'last'); 
	sa = thetan(Ia); sb = thetan(Ib);
	sample_theta = linspace(0, 6, 25);
	idxa = find(sample_theta<sa);
	idxb = find(sample_theta>sb);
	sample_theta=sample_theta([idxa,idxb]);
end
[sample_conv, sample_sprd, res_hist] = SolveNEPv(A, B, D, D1, P, 2, sample_theta, 0, mtol, 1, V00);


% ------------------------------------------------------------------------------- 
% 4. Plot convergence history and curve of spectral radius
% ------------------------------------------------------------------------------- 
% 4.1. Plot spectral radius
figure(1); 
h1 = plot(thetan, sprd_hist, '-k', 'LineWidth', 2); hold on;
h2 = plot(sample_theta, sample_conv, 'or','MarkerSize',8,'LineWidth',2); hold on;
if ex == 1
	plot([0,0], [0,1], '--k'); 
	plot([1,1], [0,1], '--k');
elseif ex == 2
	plot(thetan, ones(size(thetan)),':k')
	plot([1,1], [-1,2], '--k');
end
xlabel('$\theta$');
ylabel('convergence rate');
Lgd = legend([h1,h2], 'spectral radius', 'observed rate','Location','northeast');
Lgd.AutoUpdate = 'off'; % do not update legend 


% 4.2. Plot convergence history of SCF at manually selected theta
if ex == 1
	idxxx = [1,5,6,7,16];
   	idxmk = [1,2,3,4,5]; % theta index
   	txtx = [37,40,40,27,7];
   	txty = [0.5e-11, 0.1E-7,0.5E-3,0.5E-11,0.5E-11];
	ysftmk = [-0.05,-.05,0.05,0.05,0.05]; % marker shift for sampled theta
	xsftmk = [0.02,0.01,0.02,0.02,0.02]; % marker shift for sampled theta

elseif ex == 2
    idxxx = [1,5,9,11,16];
   	idxmk = [1,2,3,5,6]; % theta index
    txtx = [10,18,40,40,40];
    txty = [1.0e-12, 1.0E-12,0.5E-8,1.0E-2,1.0E-5];

	xsftmk = [0.02,0.02,0.02,-0.1,-0.3]; % marker shift for sampled theta
	ysftmk = [-0.05,-0.1,-0.1,-0.1,-0.1]; 

	div_theta = 3.0; % show divergged theta
	[~, idx_nc] = min(abs(thetan-div_theta));
	div_theta = thetan(idx_nc);
	sprd = sprd_hist(idx_nc);
	[~, ~, Res_hist] = SolveNEPv(A, B, D, D1, P, 2, div_theta, 0, mtol, 1, V00);
	Res = Res_hist{1};

	figure(2)
	semilogy(Res, '-', 'linewidth', 2); hold on;
	text(40,2.0E0,'$\theta_4$', 'Fontsize',15);
        
	figure(1)
   	plot(div_theta, sprd, 'xr','MarkerSize', 10, 'LineWidth',2);
	text(div_theta-0.1, sprd+0.1,'$\theta_4$','Fontsize',15);
end



for i = 1:length(idxxx)
   theta = sample_theta(idxxx(i));
   convsp = sample_conv(idxxx(i));
   Res = res_hist{idxxx(i)};
   figure(1) % mark solid circle in spectral radius plot
   plot(theta, convsp, 'or','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','r');
   text(theta+xsftmk(i), convsp+ysftmk(i),['$\theta_', num2str(idxmk(i)),'$'],'Fontsize',15);
   figure(2) % draw convergence history
   semilogy(Res, '-', 'linewidth', 2); hold on;
   text(txtx(i),txty(i),['$\theta_', num2str(idxmk(i)),'$'],'Fontsize',15);
end

xlabel('SCF iteration $i$');
ylabel('$\mbox{NRes}(X_i)$');

if ex==1
	figure(2)
	xlim([0,50])
	ylim([1.0E-13,1.0E0]) % for k=1
	figure(1)
	ylim([0,1.0]);
else
	figure(2)
	xlim([0,50])
    ylim([1.0E-13,2.0E1]) % for k=2
	figure(1)
	ylim([0,1.3]);
end

% show sampled convergence rates
disp('sampled theta, sprd, conv_est')
if ex == 1
	thetannn=0.1
elseif ex == 2
	thetannn=4.8
end
[conv_t, sprd_t] = SolveNEPv(A, B, D, D1, P, 2, thetannn, 0, mtol*100, 0, V00, 1000)


return
