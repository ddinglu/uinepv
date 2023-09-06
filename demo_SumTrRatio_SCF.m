% Test the first NEPv; 3-by-3 example, k = 1 or 2

close all; clear all;
set(0,'defaultTextInterpreter','latex');
rng(0); 

ex = 2; % ex = 1, single vector D; 2, two vectors D

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
	ssigma = -100;
	nntt = 21;
elseif ex == 2
	D =[ 
	  -1.430   2.768
	  -0.120  -0.630
	   1.098   2.229];
	n = 3; k = 2;
	ssigma = -50;
	nntt = 19;
end

corder = colororder;
maxit = 100000; mtol = 1.0E-13;
nsample = 200; 
alphan = linspace(0,1,nsample+1);


% ------------------------------------------------------------------------
% 2. Solve NEPv by SCF
% ------------------------------------------------------------------------
[V0,E0] = eig(A,B);
[~,idx1] = sort(real(diag(E0)), 'descend'); 
V00 = orth( V0(:,idx1(1:k)) ); % initial vectors

[conv_hist, sprd_hist] = SolveNEPv(A, B, D, D, eye(k), 1, alphan, ssigma, mtol, 0, V00);


% ------------------------------------------------------------------------
% 3. Estimate rates of convergence at alphas where convergence happened
% ------------------------------------------------------------------------
conv_idx = find( sprd_hist < 0.99 & sprd_hist > 0.06); % region of convergence
tt = linspace(0,1,nntt); % get evenly spaced subssamples
idxx = round((length(conv_idx) - 1)* tt) + 1;
idxx = unique(idxx); idxx = idxx(2:end); % get rid of the first entry;
alpha_idx = conv_idx(idxx); 
sample_alpha = alphan(alpha_idx);

[sample_conv, sample_sprd, res_hist] = SolveNEPv(A, B, D, D, eye(k), 1, sample_alpha, 0, mtol, 1, V00);


% ------------------------------------------------------------------------------- 
% 4. Plot convergence history and curve of spectral radius
% ------------------------------------------------------------------------------- 
% 4.1. Plot spectral radius
figure(2); 
h1 = plot(alphan, sprd_hist, '-k', 'LineWidth', 2); hold on;
h2 = plot(sample_alpha, sample_conv, 'or','MarkerSize',8,'LineWidth',2); hold on;
plot(alphan, ones(size(alphan)),':k')
xlabel('$\alpha$');
ylabel('convergence rate');
Lgd = legend([h1,h2], 'spectral radius', 'observed rate','Location','northwest');
Lgd.AutoUpdate = 'off'; % do not update legend 

% 4.2. Plot convergence history of SCF at manually selected alpha
if ex == 1
	% manually selected lines for k=1
	idxxx = [1,7,12,13,15]; 
	idxmk = [1,2,3,5,6]; % alpha index
	%txtx = [1,30,40,40,18]; for text 
	%txty = [0.5e-10, 0.5E-12,0.5E-5,0.5E-7,0.5E-12];
	txtx = [6,25,40,40,17];
	txty = [1.0e-11, 1.0E-11,0.5E-4,0.5E-7,1.0E-11];
	ysftmk = [-0.05,0,0,0,0.05]; % marker shift for sampled alpha

	alpha = 0.6; % diverged alpha
	[~, idx_nc] = min(abs(alphan-alpha));
	alpha = alphan(idx_nc);
	sprd = sprd_hist(idx_nc);
	[~, ~, Res_hist] = SolveNEPv(A, B, D, D, eye(k), 1, alpha, 0, mtol, 0, V00, 100);
	Res = Res_hist{1};

	figure(3)
	semilogy(Res, '-', 'linewidth', 2); hold on;
	%text(30,1.0E-2,['$\alpha_4=$',num2str(alpha,'%4.3f')]);
	text(40,0.5E-1,'$\alpha_4$','Fontsize',15);

	figure(2)
   	plot(alpha, sprd, 'xr','MarkerSize', 10, 'LineWidth',2);
	text(alpha, sprd - .05,'$\alpha_4$','Fontsize',15);

elseif ex == 2
	% manually selected lines for k=2
	idxxx = [1,7,9,10,12]; 
	idxmk = [1,2,3,5,6]; % alpha index
	%txtx = [3,30,43,38,20];
	%txty = [0.5e-10, 1.0E-8,0.1E-2,1.0E-4,1.0E-10];
	txtx = [13,45,45,45,23];
	txty = [1.0e-11, 1.0E-10, 3.0E-3,1.0E-4,1.0E-11];
	ysftmk = [-0.05,0,0,0,0.05]; % marker shift for sampled alpha

	alpha = 0.5; % diverged alpha
	[~, idx_nc] = min(abs(alphan-alpha));
	alpha = alphan(idx_nc);
	sprd = sprd_hist(idx_nc);
	[~, ~, Res_hist] = SolveNEPv(A, B, D, D, eye(k), 1, alpha, 0, mtol, 0, V00, 100);
	Res = Res_hist{1};

	figure(3)
	semilogy(Res, '-', 'linewidth', 2); hold on;
	%text(30,0.5E-1,['$\alpha_4=$',num2str(alpha,'%4.3f')]);
	text(45, 2.0E-1,'$\alpha_4$','Fontsize',15);

	figure(2)
   	plot(alpha, sprd, 'xr','MarkerSize', 10, 'LineWidth',2);
	text(alpha+.02, sprd,'$\alpha_4$','Fontsize',15);
end


for i = 1:length(idxxx)
   alpha = sample_alpha(idxxx(i));
   convsp = sample_conv(idxxx(i));
   Res = res_hist{idxxx(i)};
   figure(2) % mark solid circle in spectral radius plot
   plot(alpha, convsp, 'or','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','r');
   text(alpha+.02, convsp+ysftmk(i),['$\alpha_', num2str(idxmk(i)),'$'],'Fontsize',15);
   figure(3) % draw convergence history
   semilogy(Res, '-', 'linewidth', 2); hold on;
   %text(txtx(i),txty(i),['$\alpha_', num2str(idxmk(i)),'=$',num2str(alpha,'%4.3f')]);
   text(txtx(i),txty(i),['$\alpha_', num2str(idxmk(i)),'$'],'Fontsize',15);
end

xlabel('SCF iteration $i$');
ylabel('$\mbox{NRes}(X_i)$');

if ex==1
	xlim([0,50])
    ylim([1.0E-13,1.0E3]) % for k=1
else
	xlim([0,50])
    ylim([1.0E-13,1.0E2]) % for k=2
end


% ------------------------------------------------------------------------
% 5. MISC
% ------------------------------------------------------------------------
% 5.1 Show sampled convergence rates
disp('sampled alpha, sprd, conv_est')
if ex == 1
    alphannn = .46;
else
    alphannn = 0.305;
end
[conv_t, sprd_t] = SolveNEPv(A, B, D, D, eye(k), 1, alphannn, 0, mtol*10, 1, V00)

% 5.2 Show divergence interval
[~,sp_tmp] = find(sprd_hist>1);
if ~isempty(sp_tmp) 
	disp('divergence interval')
   	[alphan(sp_tmp(1)), alphan(sp_tmp(end))]
end

return;
% END