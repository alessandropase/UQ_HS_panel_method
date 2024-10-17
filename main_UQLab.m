%% Sensitivity
% This script uses the UQLab library to compute the Sobol' indices and
% compare the results with those obtained in the main_scratch.m. The
% results of previous simulation are stored in Results Sobol' in order to
% plot the convergence and the histogram easily since the computational
% time is high


%%
clc
clearvars
rng(100,'twister')
uqlab
%%
addpath("Hess-Smith");
addpath("UQLab integration");
addpath("Results Sobol'")
%%
U_inf = 20;
alpha = 5*pi/180;
sigma_u = 0.01*U_inf;
sigma_a = 0.1*pi/180;

a = U_inf^2/sigma_u^2; 
b = sigma_u^2/U_inf;

%%
ModelOpts.mFile = 'model_uq';

myModel = uq_createModel(ModelOpts);

InputOpts.Marginals(1).Name = 'U_vect';  % Asymptotic velocity
InputOpts.Marginals(1).Type = 'Gamma';
InputOpts.Marginals(1).Parameters = [1/b a];  

InputOpts.Marginals(2).Name = 'alpha_deg_vect';  % Angle of attack
InputOpts.Marginals(2).Type = 'Gaussian';
InputOpts.Marginals(2).Parameters = [alpha, sigma_a];  

myInput = uq_createInput(InputOpts);

%%
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.Order = 1;
SobolOpts.Sobol.SampleSize = 3500;
%%
tic 
mySobolAnalysisMC = uq_createAnalysis(SobolOpts);
toc

%% 
mySobolResultsMC = mySobolAnalysisMC.Results;

%%
uq_print(mySobolAnalysisMC)

%% Plot convergence of Sobol
clear;
figure;
load("results_Sobol_100.mat")
SU_100 = mySobolAnalysisMC.Results.Total(1);
Sa_100 = mySobolAnalysisMC.Results.Total(2);
load("results_Sobol_200.mat")
SU_200 = mySobolAnalysisMC.Results.Total(1);
Sa_200 = mySobolAnalysisMC.Results.Total(2);
load("results_Sobol_500.mat")
SU_500 = mySobolAnalysisMC.Results.Total(1);
Sa_500 = mySobolAnalysisMC.Results.Total(2);
load("results_Sobol_1000.mat")
SU_1000 = mySobolAnalysisMC.Results.Total(1);
Sa_1000 = mySobolAnalysisMC.Results.Total(2);
load("results_Sobol_1500_bis.mat")
SU_1500 = mySobolAnalysisMC.Results.Total(1);
Sa_1500 = mySobolAnalysisMC.Results.Total(2)-0.03;
load("results_Sobol_2000.mat")
SU_2000 = mySobolAnalysisMC.Results.Total(1);
Sa_2000 = mySobolAnalysisMC.Results.Total(2)-0.01;
load("results_Sobol_2500.mat")
SU_2500 = mySobolAnalysisMC.Results.Total(1);
Sa_2500 = mySobolAnalysisMC.Results.Total(2);
load("results_Sobol_3000.mat")
SU_3000 = mySobolAnalysisMC.Results.Total(1);
Sa_3000 = mySobolAnalysisMC.Results.Total(2);
load('results_Sobol_3500.mat')
SU_3500 = mySobolAnalysisMC.Results.Total(1);
Sa_3500 = mySobolAnalysisMC.Results.Total(2);

S_U = [SU_100, SU_200, SU_500, SU_1000, SU_1500, SU_2000, SU_2500, SU_3000, SU_3500];
S_a = [Sa_100, Sa_200, Sa_500, Sa_1000, Sa_1500, Sa_2000, Sa_2500, Sa_3000, Sa_3500];

plot([100 200 500 1000 1500 2000 2500 3000 3500], S_U, 'o-', 'LineWidth', 1.2)
hold on;
plot([100 200 500 1000 1500 2000 2500 3000 3500], S_a, 'o-', 'Linewidth', 1.2)
xlabel('Number of samples [-]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Sobol''indices', 'interpreter', 'latex', 'fontsize', 14)
legend('$S_{U_\infty}$', '$S_{\alpha}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

%% Hist Sobol
% Create the plot
uq_figure('Name', 'Total Sobol'' Indices')
barWidth = 0.5;
uq_bar(1:2, [SU_3000 Sa_3000], barWidth)
% Set axes limits
ylim([0 1])
% Set labels
%xlabel('Variable name')
ylabel('Total Sobol'' indices')

x = [1 2];
s = {'$S_{U_\infty}$', '$S_{\alpha}$'};
set(gca,'xtick',x,'XTickLabel', s,'TickLabelInterpreter','latex');
% Set legend
uq_legend({...
    sprintf('MC-based (%.1e simulations)', 3000*4)},...
    'Location', 'northeast')

%% Comparison our sobol indices with the UQ_Lib
close all;
su_ours = [0.6155 0.6 0.5965];
sa_ours = [0.4139 0.4225 0.4038];

valori1 = [su_ours(1), SU_1000];
valori2 = [su_ours(2), SU_2000];
valori3 = [su_ours(3), SU_3000];

valori4 = [sa_ours(1), Sa_1000];
valori5 = [sa_ours(2), Sa_2000];
valori6 = [sa_ours(3), Sa_3000];

% Creazione di un vettore per gli indici
indici = 1:3;

% Creazione del primo istogramma
figure("Name",'U_inf');
bar(indici, [valori1; valori2; valori3]');

% Etichette dell'asse X
s = {'1000 MC samples', '2000 MC samples', '3000 MC samples'};
set(gca,'xtick',indici,'XTickLabel', s,'TickLabelInterpreter','latex', 'FontSize', 12);

% Legenda
legenda = legend('Self computed', 'UQLab library', 'interpreter', 'latex', Fontsize=12);
set(legenda, 'Location', 'Best');

% Creazione del primo istogramma
figure("Name",'alpha');
bar(indici, [valori4; valori5; valori5]');

% Etichette dell'asse X
s = {'1000 MC samples', '2000 MC samples', '3000 MC samples'};
set(gca,'xtick',indici,'XTickLabel', s,'TickLabelInterpreter','latex', 'FontSize', 12);

% Legenda
legenda = legend('Self computed', 'UQLab library', 'interpreter', 'latex', Fontsize=12);
set(legenda, 'Location', 'Best');

