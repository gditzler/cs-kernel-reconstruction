clc;
clear;
close all;

addpath('src/')

n_avg = 250;
n = 50;
k = 10;
M = n-3;
shifts = [0,1,2,3,4];
types = {'Gaussian'};
delete(gcp('nocreate'));
parpool(50);

for t = 1:length(types)
  timez = zeros(length(shifts), M);
  errs = zeros(length(shifts), M);
  errs_no_norm = zeros(length(shifts), M);

  for i = 1:n_avg
    disp(['Running trial ',num2str(i), ' of ', num2str(n_avg)]);
    
    for q = 1:length(shifts)
      for m = 5:M
        [A, x, y] = cs_model(m, n, k, types{t});

        tic;
        [x_l1kr, x_l1] = l1kr(A, y, shifts(q));
        timez(q, m) = timez(q, m) + toc;
        errs(q, m) = errs(q, m) + per_error(x/norm(x), x_l1kr/norm(x_l1kr));
        errs_no_norm(q, m) = errs_no_norm(q, m) + per_error(x, x_l1kr);       
      end
    end
  end
  errs = errs/n_avg;
  errs_no_norm = errs_no_norm/n_avg;
  timez = timez/n_avg;

  save(['mat/',types{t}, '_shifts_expr.mat']);
end

delete(gcp('nocreate'));
