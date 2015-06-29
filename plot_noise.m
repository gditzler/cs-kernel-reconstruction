clc
clear 
close all

load mat/noise_experiments.mat
I = 5:M; %2:size(errs_no_norm,2);


% errs_noise(1, mm) = errs_noise(1, mm) + per_error(x, x_kr_spar);
% errs_noise(2, mm) = errs_noise(2, mm) + per_error(x, x_cosamp);
% errs_noise(3, mm) = errs_noise(3, mm) + per_error(x, x_omp);
% errs_noise(4, mm) = errs_noise(4, mm) + per_error(x, x_l1kr);
% errs_noise(5, mm) = errs_noise(5, mm) + per_error(x, x_l1krn);
% errs_noise(6, mm) = errs_noise(6, mm) + per_error(x, x_l1);
% errs_noise(7, mm) = errs_noise(7, mm) + per_error(x, x_l1n);

    
h = figure; 
hold on;
box on;
plot(I, errs_clean(1, 1:numel(I)), 'ro-', 'LineWidth', 2);
plot(I, errs_clean(2, 1:numel(I)), 'bs-', 'LineWidth', 2);
plot(I, errs_clean(3, 1:numel(I)), 'k^-', 'LineWidth', 2);
plot(I, errs_clean(4, 1:numel(I)), 'cp-', 'LineWidth', 2);
plot(I, errs_clean(5, 1:numel(I)), 'ms-', 'LineWidth', 2);
plot(I, errs_clean(6, 1:numel(I)), 'rs-', 'LineWidth', 2);
plot(I, errs_clean(7, 1:numel(I)), 'bo-', 'LineWidth', 2);
axis tight;
xlim([5, M])
ylim([0, 2])
legend('KR', 'CoSamp', 'OMP', 'AKRON', 'AKRONoi', 'L1', 'L1n', 'Location', 'best');
% legend('KR', 'CoSamp', 'OMP', 'AKRON', 'AKRONoi', 'Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('reconstruction error', 'FontSize', 20);
set(gca, 'fontsize', 20);
saveas(h, 'eps/noise_exp_clean.eps', 'eps2c')


g = figure; 
hold on;
box on;
plot(I, errs_noise(1, 1:numel(I)), 'ro-', 'LineWidth', 2);
plot(I, errs_noise(2, 1:numel(I)), 'bs-', 'LineWidth', 2);
plot(I, errs_noise(3, 1:numel(I)), 'k^-', 'LineWidth', 2);
plot(I, errs_noise(4, 1:numel(I)), 'cp-', 'LineWidth', 2);
plot(I, errs_noise(5, 1:numel(I)), 'ms-', 'LineWidth', 2);
plot(I, errs_noise(6, 1:numel(I)), 'rs-', 'LineWidth', 2);
plot(I, errs_noise(7, 1:numel(I)), 'bo-', 'LineWidth', 2);
xlim([5, M])
ylim([0, 2])
legend('KR', 'CoSamp', 'OMP', 'AKRON', 'AKRONoi', 'L1', 'L1n', 'Location', 'best');
% legend('KR', 'CoSamp', 'OMP', 'AKRON', 'AKRONoi', 'Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('reconstruction error', 'FontSize', 20);
set(gca, 'fontsize', 20);
saveas(g, 'eps/noise_exp_noise.eps', 'eps2c')
