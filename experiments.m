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

addpath('src/')

n_avg = 10;
n = 50;
k = 7;
M = n-1;
k_alg = 13;

errs = zeros(4, M);
timez = zeros(4, M);
sparsity = zeros(4, M);

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
    [~, ~, x_sparsest] = l0_exact_reconstruction(A, x, y);
    timez(1, m) = timez(1, m) + toc;
    errs(1, m) = errs(1, m) + per_error(x_sparsest/norm(x_sparsest), x/norm(x));
    sparsity(1, m) = sparsity(1, m) + sum(abs(x_sparsest) <= sqrt(eps))/numel(x);
    
    % CoSamp
    tic;
    x_hat = cosamp(A, y, k_alg, errFcn, opts);
    timez(2, m) = timez(2, m) + toc;
    errs(2, m) = errs(2, m) + per_error(x_hat/norm(x_hat), x/norm(x));
    sparsity(2, m) = sparsity(2, m) + sum(abs(x_hat) >= sqrt(eps))/numel(x);
    
    % OMP
    tic;
    x_omp = omp(A, y, k_alg, errFcn, opts);
    timez(3, m) = timez(3, m) + toc;
    errs(3, m) = errs(3, m) + per_error(x_omp/norm(x_omp), x/norm(x));
    sparsity(3, m) = sparsity(3, m) + sum(abs(x_omp) >= sqrt(eps))/numel(x);
    
    % L1-Approx of KR
    tic;
    x_l1kr = l1_approximate_reconstruction(A, y);
    timez(4, m) = timez(4, m) + toc;
    errs(4, m) = errs(4, m) + per_error(x_l1kr/norm(x_l1kr), x/norm(x));
    sparsity(4, m) = sparsity(4, m) + sum(abs(x_l1kr) >= sqrt(eps))/numel(x);
  end
end
errs = errs/n_avg;
timez = timez/n_avg;
sparsity = sparsity/n_avg;


n_avg = 10;
n = 50;
k = 7;
M = n-1;
k_alg = 13;

save(['mat/gaussian_reconstruction_n', num2str(n), 'k', num2str(k), ...
  'ka',num2str(k_alg),'.mat']);

h = figure; 
hold on;
box on;
plot(errs(1,:), 'ro-', 'LineWidth', 2);
plot(errs(2,:), 'bs-', 'LineWidth', 2);
plot(errs(3,:), 'k^-', 'LineWidth', 2);
plot(errs(4,:), 'cp-', 'LineWidth', 2);
xlim([1, M])
legend('KR', 'CoSamp', 'OMP', 'L1KR','Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('reconstruction error', 'FontSize', 20);
set(gca, 'fontsize', 20);

g = figure; 
hold on;
box on;
plot(timez(1,:), 'ro-', 'LineWidth', 2);
plot(timez(2,:), 'bs-', 'LineWidth', 2);
plot(timez(3,:), 'k^-', 'LineWidth', 2);
plot(timez(4,:), 'cp-', 'LineWidth', 2);
xlim([1, M])
legend('KR', 'CoSamp', 'OMP', 'L1KR','Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('time (s)', 'FontSize', 20);
set(gca, 'fontsize', 20);

q = figure; 
hold on;
box on;
plot(sparsity(1,:), 'ro-', 'LineWidth', 2);
plot(sparsity(2,:), 'bs-', 'LineWidth', 2);
plot(sparsity(3,:), 'k^-', 'LineWidth', 2);
plot(sparsity(4,:), 'cp-', 'LineWidth', 2);
xlim([1, M])
legend('KR', 'CoSamp', 'OMP', 'L1KR','Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('sparsity', 'FontSize', 20);
set(gca, 'fontsize', 20);

delete(gcp('nocreate'));
