%%
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

set(0,'RecursionLimit', 10000);

addpath('src/')

n_avg = 250;
n_set = [50:10:100];
k_set = floor(.05*n_set);
M = 20;
m = 30;
k_alg_set = floor(.1*n_set);
types = {'Gaussian', 'Uni',   'Bernoulli'};
opts.printEvery = 10000000;
errFcn = [];
epsilon = 0.05;

delete(gcp('nocreate'));
parpool(20);

for t = 1:length(types)
  errs = zeros(7, length(n_set));
  errs_no_norm = zeros(7, length(n_set));
  timez = zeros(7, length(n_set));
  sparsity = zeros(7, length(n_set));

  

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
      errs(q, j) = errs(q, j) + per_error(x/norm(x), x_hat/norm(x_hat));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x, x_hat);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_hat) >= sqrt(eps))/numel(x);
      q = q+1;

      % OMP
      tic;
      x_omp = omp(A, y, k_alg, errFcn, opts);
      timez(q, j) = timez(q, j) + toc;
      errs(q, j) = errs(q, j) + per_error(x/norm(x), x_omp/norm(x_omp));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x, x_omp);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_omp) >= sqrt(eps))/numel(x);
      q = q+1;

      % L1-Approx of KR
      tic;
      [x_l1kr, x_l1] = l1kr(A, y);
      timez(q, j) = timez(q, j) + toc;
      errs(q, j) = errs(q, j) + per_error(x/norm(x), x_l1kr/norm(x_l1kr));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x, x_l1kr);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_l1kr) >= sqrt(eps))/numel(x);
      q = q+1;
      
      errs(q, j) = errs(q, j) + per_error(x/norm(x), x_l1/norm(x_l1));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x, x_l1);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_l1) >= sqrt(eps))/numel(x);
      q = q+1;
      
      % L1Noise-Approx of KR
      tic;
      [x_l1kr, x_l1] = l1kr_noise(A, y, epsilon);
      timez(q, j) = timez(q, j) + toc;
      errs(q, j) = errs(q, j) + per_error(x/norm(x), x_l1kr/norm(x_l1kr));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x,x_l1kr);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_l1kr) >= sqrt(eps))/numel(x);
      q = q+1;
      
      errs(q, j) = errs(q, j) + per_error(x/norm(x), x_l1/norm(x_l1));
      errs_no_norm(q, j) = errs_no_norm(q, j) + per_error(x, x_l1);
      sparsity(q, j) = sparsity(q, j) + sum(abs(x_l1) >= sqrt(eps))/numel(x);
      
    end
    
  end
  
  errs = errs/n_avg;
  errs_no_norm = errs_no_norm/n_avg;
  sparsity = sparsity/n_avg;
  timez = timez/n_avg;
  
  save(['mat/large_',types{t}, 'k', num2str(k), '.mat']);
end

delete(gcp('nocreate'));
