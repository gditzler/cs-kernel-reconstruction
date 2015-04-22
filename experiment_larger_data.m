%%
%  AUTHORS
%    Belhassan Bayer, Nidhal Bouynaya, and Gregory Ditzler
%
%  MAINTAINER
%    Gregory Ditzler (gregory.ditzler@gmail.coj)
%
%  LICENSE
%    MIT
%% 
clc;
clear;
close all;

addpath('src/')

n_avg = 250;
n_set = 100:20:2000;
k_set = floor(.05*n_set);
M = 20;
m = 30;
k_alg_set = floor(.1*n_set);
types = {'Gaussian', 'Uni',   'Bernoulli'};
% delete(gcp('nocreate'));
% parpool(50);

for t = 1:length(types)
  errs = zeros(7, length(n_set));
  errs_no_norm = zeros(7, length(n_set));
  timez = zeros(7, length(n_set));
  sparsity = zeros(7, length(n_set));

  opts.printEvery = 10000000;
  errFcn = [];

  for i = 1:n_avg
    disp(['Running trial ',num2str(i), ' of ', num2str(n_avg)]);
    
    for j = 1:length(n_set)
      n = n_set(j);
      k = k_set(j);
      k_alg = k_alg_set(j);
      
      [A, x, y] = cs_model(m, n, k, types{t});
      q = 1;

      % CoSamp
      tic;
      x_hat = cosamp(A, y, k_alg, errFcn, opts);
      timez(q, j) = timez(q, j) + toc;
      errs(q, j) = errs(q, j) + per_error(x_hat/norm(x_hat), x/norm(x));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x_hat, x);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_hat) >= sqrt(eps))/numel(x);
      q = q+1;

      % OMP
      tic;
      x_omp = omp(A, y, k_alg, errFcn, opts);
      timez(q, j) = timez(q, j) + toc;
      errs(q, j) = errs(q, j) + per_error(x_omp/norm(x_omp), x/norm(x));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x_omp, x);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_omp) >= sqrt(eps))/numel(x);
      q = q+1;

      % L1-Approx of KR
      tic;
      [x_l1kr, x_l1] = l1kr(A, y);
      timez(q, j) = timez(q, j) + toc;
      errs(q, j) = errs(q, j) + per_error(x_l1kr/norm(x_l1kr), x/norm(x));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x_l1kr, x);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_l1kr) >= sqrt(eps))/numel(x);
      q = q+1;
      
      errs(q, j) = errs(q, j) + per_error(x_l1/norm(x_l1), x/norm(x));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x_l1, x);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_l1) >= sqrt(eps))/numel(x);
      
    end
    
  end
  
  save(['mat/large_',types{t}, 'k', num2str(k), '.mat']);
end

delete(gcp('nocreate'));
