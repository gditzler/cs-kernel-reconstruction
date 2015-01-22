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

n_avg = 30;
n = 20;
k = 7;
M = n-3;
k_alg = 15;
types = {'Gaussian', 'GaussianShift'};


for t = 1:length(types)
  errs = zeros(6, M);
  errs_no_norm = zeros(6, M);
  timez = zeros(6, M);
  sparsity = zeros(6, M);

  opts.printEvery = 10000000;
  errFcn = [];

  delete(gcp('nocreate'));
  parpool(60);

  for i = 1:n_avg
    disp(['Running trial ',num2str(i), ' of ', num2str(n_avg)]);

    for m = 1:M
      [A, x, y] = cs_model(m, n, k, types{t});
      q = 1;

      % KR
      tic;
      [~, x_kr_err, x_kr_spar] = l0_exact_reconstruction(A, x, y);
      timez(q, m) = timez(q, m) + toc;
      errs(q, m) = errs(q, m) + per_error(x_kr_err/norm(x_kr_err), x/norm(x));
      errs_no_norm(q, m) = errs_no_norm(q, m) + per_error(x_kr_err, x);
      sparsity(q, m) = sparsity(q, m) + sum(abs(x_kr_err) <= sqrt(eps))/numel(x);
      q = q+1;

      errs(q, m) = errs(q, m) + per_error(x_kr_spar/norm(x_kr_spar), x/norm(x));
      errs_no_norm(q, m) = errs_no_norm(q, m) + per_error(x_kr_spar, x);
      sparsity(q, m) = sparsity(q, m) + sum(abs(x_kr_spar) <= sqrt(eps))/numel(x);
      q = q+1;

      % CoSamp
      tic;
      x_hat = cosamp(A, y, k_alg, errFcn, opts);
      timez(q, m) = timez(q, m) + toc;
      errs(q, m) = errs(q, m) + per_error(x_hat/norm(x_hat), x/norm(x));
      errs_no_norm(q, m) = errs_no_norm(q, m) + per_error(x_hat, x);
      sparsity(q, m) = sparsity(q, m) + sum(abs(x_hat) >= sqrt(eps))/numel(x);
      q = q+1;

      % OMP
      tic;
      x_omp = omp(A, y, k_alg, errFcn, opts);
      timez(q, m) = timez(q, m) + toc;
      errs(q, m) = errs(q, m) + per_error(x_omp/norm(x_omp), x/norm(x));
      errs_no_norm(q, m) = errs_no_norm(q, m) + per_error(x_omp, x);
      sparsity(q, m) = sparsity(q, m) + sum(abs(x_omp) >= sqrt(eps))/numel(x);
      q = q+1;

      % L1-Approx of KR
      tic;
      x_l1kr = l1kr(A, y);
      timez(q, m) = timez(q, m) + toc;
      errs(q, m) = errs(q, m) + per_error(x_l1kr/norm(x_l1kr), x/norm(x));
      errs_no_norm(q, m) = errs_no_norm(q, m) + per_error(x_l1kr, x);
      sparsity(q, m) = sparsity(q, m) + sum(abs(x_l1kr) >= sqrt(eps))/numel(x);
    end
  end
  errs = errs/n_avg;
  timez = timez/n_avg;
  sparsity = sparsity/n_avg;


  save(['mat/',types{t}, '_n', num2str(n), 'k', num2str(k), ...
    'ka',num2str(k_alg),'.mat']);
end

delete(gcp('nocreate'));
