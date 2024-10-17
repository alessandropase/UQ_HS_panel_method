% This script perform the sampling and then a MC simulation to evaluate the
% impact of the uncertainty on the output. Then the plots are generated.
% The final part is related to the computation of Sobol' indices with out
% the usage of external library

%% Intro
clc;
clear;
close all;

%%
addpath("Hess-Smith");
addpath("UQLab integration");

%% Montecarlo simulation
N_samples = 3000;
L_vect = zeros(1, N_samples);
mean_L = zeros(1, N_samples);
sigma_L = zeros(1, N_samples);
coeff_var = zeros(1, N_samples);
Samples = zeros(2, N_samples);
U_inf = 20;
alpha = 5*pi/180;
sigma_u = 0.01*U_inf;
sigma_a = 0.1*pi/180;

a = U_inf^2/sigma_u^2; 
b = sigma_u^2/U_inf;

L_nominal = model_fun(U_inf, rad2deg(alpha));
fprintf("\nNominal value = %.2f N/m \n", L_nominal);

%%
tic
for i = 1 : N_samples
    clc;
    perc = i/N_samples*100;
    fprintf('Percentage of the samples: %.2f%%\n', perc);
    U_inf_s = gamrnd(a,b);
    alpha_s = 180/pi*circ_vmrnd(alpha, 1/(sigma_a^2), 1);
    Samples(1, i) = U_inf_s;
    Samples(2, i) = alpha_s;
    L = model_fun(U_inf_s, alpha_s);
    L_vect(i) = L;
    mean_L(i) = mean(L_vect(1:i));
    sigma_L(i) = var(L_vect(1:i));
    coeff_var(i) = std(L_vect(1:i))/abs(mean_L(i));

    if i == 50
        [f_50,xi_50] = ksdensity(L_vect(1:i));
    elseif i == 300
        [f_200,xi_200] = ksdensity(L_vect(1:i)); 
    elseif i == 500
        [f_500,xi_500] = ksdensity(L_vect(1:i)); 
    elseif i == 1000
        [f_1000,xi_1000] = ksdensity(L_vect(1:i));
    elseif i == 1500
        [f_1500,xi_1500] = ksdensity(L_vect(1:i));
    elseif i == 2000
        [f_2000,xi_2000] = ksdensity(L_vect(1:i));
    elseif i == 2500
        [f_2500,xi_2500] = ksdensity(L_vect(1:i));
    end

end
toc;

%% Plot for the first 150 samples
U_nom = U_inf;
U_vect_plot = (linspace(18, 22, 1001));
pdf_U = gampdf(U_vect_plot, a, b); 

alpha_vect_plot = (linspace(3.5, 6.5, 1001));
pdf_alpha = normpdf(alpha_vect_plot, alpha*180/pi, sigma_a*180/pi); 
%pdf_alpha = circ_vmpdf(alpha_vect_plot./180.*pi, alpha, sigma_a); 
[xx,yy] = meshgrid(U_vect_plot,alpha_vect_plot);

[pdf_uu,pdf_aa] = meshgrid(pdf_U, pdf_alpha);

pdf_au = pdf_uu.*pdf_aa;

figure; 
ax = gca;
ax.FontSize = 16; 
contour(xx,yy,pdf_au);
hold on
for i = 1:100
    plot(Samples(1,i), Samples(2,i), 'bo', 'MarkerSize', 3, 'LineWidth', 1.3);
end

ylabel("$\alpha$ [deg]", 'Interpreter', "latex", 'FontSize', 14)
xlabel("$U_\infty$", 'Interpreter', 'latex', 'FontSize', 14)
xlim([19.25 20.75]);
ylim([4.5 5.5])
grid off;

%% Plot for MC simulations
close all;
figure;
plot(1 : N_samples, L_vect, 'LineWidth',1.3);
grid on;
box on;
xlabel('Iterations [-]', 'Interpreter', 'latex','fontsize', 14)
ylabel('Lift sample [N/m]','Interpreter', 'latex', 'fontsize', 14)

figure;
plot(1 : N_samples, mean_L, 'LineWidth',1.3);
xlabel('Iterations [-]', 'Interpreter', 'latex', 'fontsize', 14)
ylabel('Mean of the samples [N/m]','Interpreter', 'latex', 'fontsize', 14)
grid on;

figure;
plot(1 : N_samples, sigma_L, 'LineWidth',1.3);
xlabel('Iterations [-]', 'Interpreter', 'latex','fontsize', 14);
ylabel('Variance of the sample [N/m]','Interpreter', 'latex', 'fontsize', 14);
grid on;

figure;
plot(1 : N_samples, coeff_var, 'LineWidth',1.3);
xlabel('Iterations [-]', 'Interpreter', 'latex','fontsize', 14)
ylabel('Coefficient of variation','Interpreter', 'latex', 'fontsize', 14)
grid on;


figure;
histogram(L_vect, 50, 'Normalization','pdf', 'FaceAlpha',0.3);
xlabel('Lift [N/m]', 'Interpreter', 'latex','fontsize', 14)
ylabel('Kernel density estimation','Interpreter', 'latex', 'fontsize', 14)
[f,xi] = ksdensity(L_vect); 
hold on;
plot(xi,f, 'LineWidth', 2, 'Color', '#000080');

%% Kernel density estimation convergence
figure;
plot(50*f_50 + 0.5, xi_50, 'LineWidth', 1.3);
hold on;
plot(50*f_200 + 2, xi_200, 'LineWidth', 1.3);
plot(50*f_500 + 5, xi_500, 'LineWidth', 1.3);
plot(50*f_1000 + 10, xi_1000, 'LineWidth', 1.3)
plot(50*f_1500 + 15, xi_1500, 'LineWidth', 1.3)
plot(50*f_2000 + 20, xi_2000, 'LineWidth', 1.3)
plot(50*f_2500 + 25, xi_2500, 'LineWidth', 1.3)
plot(50*f + 30, xi, 'LineWidth', 1.3, 'Color','blu')
%plot([1 : N_samples]./100, mean_L, 'LineWidth',1.3);
grid on;
grid minor;
ylabel('Lift [N/m]', 'FontSize',14, 'Interpreter','latex')
xlabel('Hundreds of iterations [-]', 'Interpreter','latex', 'fontsize', 14)
%legend('PDF at 50 iterations', 'PDF at 200 iterations', 'PDF at 500 iterations', 'PDF at 1000 iterations', 'PDF at 1500 iterations', 'PDF at 2000 iterations')

%% Sobol' indices
% Montecarlo sampling
% Input sampling
N_samples = 3000;
U_vect = zeros(1, N_samples);
alpha_vec = zeros(1, N_samples);
U_tilde_vec = zeros(1, N_samples);
alpha_tilde_vec = zeros(1, N_samples);
L_vect = zeros(1, N_samples);
mean_L = zeros(1, N_samples);
sigma_L = zeros(1, N_samples);
coeff_var = zeros(1, N_samples);
U_inf = 20;
alpha = 5*pi/180;
sigma_u = 0.01*U_inf;
sigma_a = 0.1*pi/180;
a = U_inf^2/sigma_u^2; 
b = sigma_u^2/U_inf;
for i = 1 : N_samples
    clc;
    perc = i/N_samples*100;
    fprintf('Percentage of the samples: %.2f%%\n', perc);
    U_inf_s = gamrnd(a,b);
    U_vec(i) = U_inf_s;
    U_inf_s_bis = gamrnd(a,b);
    U_tilde_vec(i) = U_inf_s_bis;
    alpha_s = 180/pi*circ_vmrnd(alpha, 1/(sigma_a^2), 1);
    alpha_vec(i) = alpha_s;
    alpha_s_bis = 180/pi*circ_vmrnd(alpha, 1/(sigma_a^2), 1);
    alpha_tilde_vec(i) = alpha_s_bis;
 end
% Output sampling
G_U_alpha = zeros(1, N_samples);
G_Utilde_alpha = zeros(1, N_samples);
G_U_alphatilde = zeros(1, N_samples);
for i = 1:N_samples
    clc;
    perc = i/N_samples*100;
    fprintf('Percentage of the output sampling: %.2f%%\n', perc);
    G_U_alpha(i) = model_fun(U_vec(i), alpha_vec(i));
    G_Utilde_alpha(i) = model_fun(U_tilde_vec(i), alpha_vec(i));
    G_U_alphatilde(i) = model_fun(U_vec(i), alpha_tilde_vec(i));
end
% Sobolev indices integration
mean1 = mean(G_U_alpha); 
mean2 = mean(G_Utilde_alpha);
mean3 = mean(G_U_alphatilde);
s_U = 0;
s_alpha = 0;
nu = N_samples;
for l = 1:nu
    s_U = s_U + 1/nu*(G_U_alpha(l)-mean1)*(G_U_alphatilde(l)-mean2);
    s_alpha = s_alpha + 1/nu*(G_U_alpha(l)-mean1)*(G_Utilde_alpha(l)-mean3);
    s_U_vec(l) = s_U;
    s_alpha_vec(l) = s_alpha;
end

sigma_sqrd = var(G_U_alpha);
s_U_vec = s_U_vec./sigma_sqrd;
s_alpha_vec = s_alpha_vec./sigma_sqrd;
S_U = s_U_vec(end);
S_alpha = s_alpha_vec(end);

