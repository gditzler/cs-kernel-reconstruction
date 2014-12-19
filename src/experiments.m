%%
%  AUTHORS
%    Belhassan Bayer, Nidhal Bouynaya, and Gregory Ditzler
%
%  MAINTAINER
%    Gregory Ditzler (gregory.ditzler@gmail.com)
%
%  LICENSE
%    MIT
%% 
clc;
clear;
close all;

n_avg = 20;
n = 20;
k = 7;
M = 16;

errs = zeros(3,M);
timez = zeros(3,M);

opts.printEvery = 10000000;
errFcn = [];

delete(gcp('nocreate'));
parpool(50);

for i = 1:n_avg
  disp(['Running trial ',num2str(i), ' of ', num2str(n_avg)]);
  
  for m = 1:M
    [A, x, y] = cs_model(m, n, k);
    
    % KR
    tic;
    errs(1, m) = errs(1, m) + l0_exact_reconstruction(A, x, y);
    timez(1, m) = timez(1, m) + toc;
    
    % CoSamp
    tic;
    x_hat = cosamp(A, y, k, errFcn, opts);
    timez(2, m) = timez(2, m) + toc;
    errs(2, m) = errs(2, m) + per_error(x_hat, x);
    
    % OMP
    tic;
    x_omp = omp(A, y, k, errFcn, opts);
    timez(3, m) = timez(3, m) + toc;
    errs(3, m) = errs(3, m) + per_error(x_omp, x);
  end
end
errs = errs/n_avg;
timez = timez/n_avg;

h = figure; 
hold on;
box on;
plot(errs(1,:), 'ro-', 'LineWidth', 2);
plot(errs(2,:), 'bs-', 'LineWidth', 2);
plot(errs(3,:), 'k^-', 'LineWidth', 2);
xlim([1,M])
legend('KR', 'CoSamp', 'OMP', 'Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('reconstruction error', 'FontSize', 20);
set(gca, 'fontsize', 20);

save('../mat/gaussian_reconstruction_n20k7.mat');
delete(gcp('nocreate'));
