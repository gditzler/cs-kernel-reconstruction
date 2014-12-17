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
errs = zeros(2,16);

opts.printEvery = 10000000;
errFcn = [];

delete(gcp('nocreate'));
parpool(25);

for i = 1:n_avg
  disp(['Running trial ',num2str(i), ' of ', num2str(n_avg)]);
  
  for m = 1:16
    [A, x, y] = cs_model(m, n, k);
    
    errs(1, m) = errs(1, m) + l0_exact_reconstruction(A, x, y);
    
    x_hat = cosamp(A, y, k, errFcn, opts);
    errs(2, m) = errs(2, m) + per_error(x_hat, x);
  end
end
errs = errs/n_avg;
save('../mat/gaussian_reconstruction_n20k7.mat');
delete(gcp('nocreate'));
