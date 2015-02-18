%% plot gaussian_reconstruction_n20k7ka7.mat
clc;
clear;
close all;

load mat/Gaussian_n20k7ka9.mat
% load mat/GaussianShift_n20k7ka9.mat

I = 2:size(errs_no_norm,2);

% 1) kr (min error)
% 2) kr (min sparse)
% 3) cosamp
% 4) omp
% 5) l1kr
% 6) l1

h = figure; 
hold on;
box on;
plot(I, errs_no_norm(2, I), 'ro-', 'LineWidth', 2);
plot(I, errs_no_norm(3, I), 'bs-', 'LineWidth', 2);
plot(I, errs_no_norm(4, I), 'k^-', 'LineWidth', 2);
plot(I, errs_no_norm(5, I), 'cp-', 'LineWidth', 2);
plot(I, errs_no_norm(6, I), 'ms-', 'LineWidth', 2);
xlim([1, M])
ylim([0, 1.5*max(errs_no_norm(2,I))])
axis tight;
legend('KR', 'CoSamp', 'OMP', 'L1KR', 'L1', 'Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('reconstruction error', 'FontSize', 20);
set(gca, 'fontsize', 20);
saveas(h, ['eps/error_n',num2str(n),'k',num2str(k),'ka',num2str(k_alg),'.eps'], 'eps2c')

g = figure; 
hold on;
box on;
plot(I, timez(1, I), 'ro-', 'LineWidth', 2);
plot(I, timez(3, I), 'bs-', 'LineWidth', 2);
plot(I, timez(4, I), 'k^-', 'LineWidth', 2);
plot(I, timez(5, I), 'cp-', 'LineWidth', 2);
xlim([1, M])
legend('KR', 'CoSamp', 'OMP', 'L1KR/L1', 'Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('time (s)', 'FontSize', 20);
axis tight;
set(gca, 'fontsize', 20);
saveas(g, ['eps/times_n',num2str(n),'k',num2str(k),'ka',num2str(k_alg),'.eps'], 'eps2c')

q = figure; 
hold on;
box on;
plot(I, sparsity(2, I), 'ro-', 'LineWidth', 2);
plot(I, sparsity(3, I), 'bs-', 'LineWidth', 2);
plot(I, sparsity(4, I), 'k^-', 'LineWidth', 2);
plot(I, sparsity(5, I), 'cp-', 'LineWidth', 2);
plot(I, sparsity(6, I), 'ms-', 'LineWidth', 2);
xlim([1, M])
axis tight;
set(gca, 'fontsize', 20);
legend('KR', 'CoSamp', 'OMP', 'L1KR', 'L1', 'Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('sparsity', 'FontSize', 20);
saveas(q, ['eps/sparsity_n',num2str(n),'k',num2str(k),'ka',num2str(k_alg),'.eps'], 'eps2c')
