%% plot gaussian_reconstruction_n20k7ka7.mat
clc;
clear;
close all;

load mat/gaussian_reconstruction_n20k7ka7.mat

h = figure; 
hold on;
box on;
plot(errs_no_norm(2,:), 'ro-', 'LineWidth', 2);
plot(errs_no_norm(3,:), 'bs-', 'LineWidth', 2);
plot(errs_no_norm(4,:), 'k^-', 'LineWidth', 2);
plot(errs_no_norm(5,:), 'cp-', 'LineWidth', 2);
xlim([1, M])
ylim([0, 1.5*max(errs_no_norm(2,:))])
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
plot(sparsity(2,:), 'ro-', 'LineWidth', 2);
plot(sparsity(3,:), 'bs-', 'LineWidth', 2);
plot(sparsity(4,:), 'k^-', 'LineWidth', 2);
plot(sparsity(5,:), 'cp-', 'LineWidth', 2);
xlim([1, M])
legend('KR', 'CoSamp', 'OMP', 'L1KR','Location', 'best');
xlabel('m', 'FontSize', 20);
ylabel('sparsity', 'FontSize', 20);
set(gca, 'fontsize', 20);
