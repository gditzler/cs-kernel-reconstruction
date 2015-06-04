clc;
clear;
close all;

addpath('src/')

n_avg = 250;
n = 20;
k = 9;
M = n-3;
k_alg = 9;
epsilon = 0.05;
types = {'Gaussian', 'Uni', 'Bernoulli'};
delete(gcp('nocreate'));
parpool(50);


errs_clean = zeros(7, M);
errs_noise = zeros(7, M);
opts.printEvery = 10000000;
errFcn = [];

for i = 1:n_avg
  disp(['Running trial ',num2str(i), ' of ', num2str(n_avg)]);

  mm = 1;
  for m = 5:M


    [A, x, y] = cs_model(m, n, k, 'Gaussian');
    [~, ~, x_kr_spar] = l0_exact_reconstruction(A, x, y);
    x_cosamp = cosamp(A, y, k_alg, errFcn, opts);
    x_omp = omp(A, y, k_alg, errFcn, opts);
    [x_l1kr, x_l1] = l1kr(A, y);
    [x_l1krn, x_l1n] = l1kr_noise(A, y, epsilon);

    errs_clean(1, mm) = errs_clean(1, mm) + per_error(x, x_kr_spar);
    errs_clean(2, mm) = errs_clean(2, mm) + per_error(x, x_cosamp);
    errs_clean(3, mm) = errs_clean(3, mm) + per_error(x, x_omp);
    errs_clean(4, mm) = errs_clean(4, mm) + per_error(x, x_l1kr);
    errs_clean(5, mm) = errs_clean(5, mm) + per_error(x, x_l1krn);
    errs_clean(6, mm) = errs_clean(6, mm) + per_error(x, x_l1);
    errs_clean(7, mm) = errs_clean(7, mm) + per_error(x, x_l1n);

    [A, x, y] = cs_model(m, n, k, 'GaussianNoise');
    [~, x_kr_err, x_kr_spar] = l0_exact_reconstruction(A, x, y);
    x_hat = cosamp(A, y, k_alg, errFcn, opts);
    x_omp = omp(A, y, k_alg, errFcn, opts);
    [x_l1kr, x_l1] = l1kr(A, y);
    [x_l1krn, x_l1n] = l1kr_noise(A, y, epsilon);

    errs_noise(1, mm) = errs_noise(1, mm) + per_error(x, x_kr_spar);
    errs_noise(2, mm) = errs_noise(2, mm) + per_error(x, x_cosamp);
    errs_noise(3, mm) = errs_noise(3, mm) + per_error(x, x_omp);
    errs_noise(4, mm) = errs_noise(4, mm) + per_error(x, x_l1kr);
    errs_noise(5, mm) = errs_noise(5, mm) + per_error(x, x_l1krn);
    errs_noise(6, mm) = errs_noise(6, mm) + per_error(x, x_l1);
    errs_noise(7, mm) = errs_noise(7, mm) + per_error(x, x_l1n);

    mm = mm+1;
  end
end
errs_clean = errs_clean/n_avg;
errs_noise = errs_noise/n_avg;

save('mat/noise_experiments.mat');
