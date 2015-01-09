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
n = 20;
k = 7;
M = 16;
k_alg = 7;

errs = zeros(6, M);
timez = zeros(6, M);
sparsity = zeros(6, M);

opts.printEvery = 10000000;
errFcn = [];

delete(gcp('nocreate'));
parpool(6);

for i = 1:n_avg
  disp(['Running trial ',num2str(i), ' of ', num2str(n_avg)]);
  
  for m = 1:M
    [A, x, y] = cs_model(m, n, k);
    q = 1;
    
    % KR
    tic;
    [~, x_kr_err, x_kr_spar] = l0_exact_reconstruction(A, x, y);
    timez(q, m) = timez(q, m) + toc;
    errs(q, m) = errs(q, m) + per_error(x_kr_err/norm(x_kr_err), x/norm(x));
    sparsity(q, m) = sparsity(q, m) + sum(abs(x_kr_err) <= sqrt(eps))/numel(x);
    q = q+1;
    
    errs(q, m) = errs(q, m) + per_error(x_kr_spar/norm(x_kr_spar), x/norm(x));
    sparsity(q, m) = sparsity(q, m) + sum(abs(x_kr_spar) <= sqrt(eps))/numel(x);
    q = q+1;
    
    % CoSamp
    tic;
    x_hat = cosamp(A, y, k_alg, errFcn, opts);
    timez(q, m) = timez(q, m) + toc;
    errs(q, m) = errs(q, m) + per_error(x_hat/norm(x_hat), x/norm(x));
    sparsity(q, m) = sparsity(q, m) + sum(abs(x_hat) >= sqrt(eps))/numel(x);
    q = q+1;
    
    % OMP
    tic;
    x_omp = omp(A, y, k_alg, errFcn, opts);
    timez(q, m) = timez(q, m) + toc;
    errs(q, m) = errs(q, m) + per_error(x_omp/norm(x_omp), x/norm(x));
    sparsity(q, m) = sparsity(q, m) + sum(abs(x_omp) >= sqrt(eps))/numel(x);
    q = q+1;
    
    % L1-Approx of KR
    tic;
    [x_l1kr, x_l1] = l1_approximate_reconstruction(A, y);
    timez(q, m) = timez(q, m) + toc;
    errs(q, m) = errs(q, m) + per_error(x_l1kr/norm(x_l1kr), x/norm(x));
    sparsity(q, m) = sparsity(q, m) + sum(abs(x_l1kr) >= sqrt(eps))/numel(x);
    q = q+1;
    
    errs(q, m) = errs(q, m) + per_error(x_l1/norm(x_l1), x/norm(x));
    sparsity(q, m) = sparsity(q, m) + sum(abs(x_l1) >= sqrt(eps))/numel(x);
  end
end
errs = errs/n_avg;
timez = timez/n_avg;
sparsity = sparsity/n_avg;


save(['mat/gaussian_reconstruction_n', num2str(n), 'k', num2str(k), ...
  'ka',num2str(k_alg),'.mat']);

h = figure; 
hold on;
box on;
plot(errs(1,:), 'ro-', 'LineWidth', 2);
plot(errs(3,:), 'bs-', 'LineWidth', 2);
plot(errs(4,:), 'k^-', 'LineWidth', 2);
plot(errs(5,:), 'cp-', 'LineWidth', 2);
xlim([1, M])
legend('KR', 'CoSamp', 'OMP', 'L1KR','Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('reconstruction error', 'FontSize', 20);
set(gca, 'fontsize', 20);

g = figure; 
hold on;
box on;
plot(timez(1,:), 'ro-', 'LineWidth', 2);
plot(timez(3,:), 'bs-', 'LineWidth', 2);
plot(timez(4,:), 'k^-', 'LineWidth', 2);
plot(timez(5,:), 'cp-', 'LineWidth', 2);
xlim([1, M])
legend('KR', 'CoSamp', 'OMP', 'L1KR','Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('time (s)', 'FontSize', 20);
set(gca, 'fontsize', 20);

q = figure; 
hold on;
box on;
plot(sparsity(1,:), 'ro-', 'LineWidth', 2);
plot(sparsity(3,:), 'bs-', 'LineWidth', 2);
plot(sparsity(4,:), 'k^-', 'LineWidth', 2);
plot(sparsity(5,:), 'cp-', 'LineWidth', 2);
xlim([1, M])
legend('KR', 'CoSamp', 'OMP', 'L1KR','Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('sparsity', 'FontSize', 20);
set(gca, 'fontsize', 20);

delete(gcp('nocreate'));
