%% Clear Workspace
close all
clear
clc

%% Range
x_max = 130;
x = linspace(0, x_max, 100)';

%% Ideal Heaviside Function
r = 91.7;
alpha_h = 3.333e-12;
y = alpha_h * (x - r) .* heaviside(x - r);
u_max = alpha_h * (x_max - r);

%% Hill Function
n = 4;
kappa = 300;
alpha = u_max*(1 + (kappa/x_max)^n);
% u = HillFunction(x, alpha, kappa, n);

%% Fit
ft = fittype('HillFunction(x, alpha, kappa, n)');
f = fit(x, y, ft, 'StartPoint', [alpha, kappa, n]);

%% Plot
figure();
plot(f, x, y); 
hold on

n = 4;
kappa = 500;
alpha = 8e-9;
u = HillFunction(x, alpha, kappa, n);

plot(x, u);