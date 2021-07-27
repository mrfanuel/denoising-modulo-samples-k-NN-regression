%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for the modulo denoising 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
addpath(genpath('lib'))
addpath(genpath('data'))

%%%%%%%%% Pre-processing %%%%%%%%%%%%%%%%%% 
load('Vesuvius_zoom.mat');

%% removing some part of the sea
x0 = x0(41:end);
y0 = y0(41:end);
alt0 = alt0(41:end,41:end);

alt0(alt0<0)=0; % removing negative entries

% change the scale
alt = double(alt0)/500.; 


%%%%%%%%% Defining the grid %%%%%%%%%%%%%%%%%% 

%%%% change latitude and longitude to [0,1] grid %%%%

m = length(x0);
a = 0;    b = 1; 
x0 = (a:((b-a)/(m-1)):b)'; % mx1 vector 

[X,Y] = meshgrid(x0,x0);
x = [X(:), Y(:)];

n = m^2; %% number of samples is m^2
sigma = 0.1; %% noise level 
k = 40;     %% number of neigbours

%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%% 

f_clean_grid = alt;
f_clean = reshape(f_clean_grid,[n,1]);

f_mod1_clean = mod(f_clean,1);
f_mod1_clean_grid = mod(f_clean_grid,1);
f_noise = f_clean + sigma * randn(n,1); 

y = mod(f_noise,1);  
y_grid = reshape(y,[m,m]);
y_clean = reshape(f_mod1_clean,[m,m]);

z = exp(1i*2*pi*y);  % nx1 vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In-sample knn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gest_kNN = kNN_denoise(z,x,k);
gest_kNN_proj = project_manifold(gest_kNN);
f_mod1_denoised = extract_modulo(gest_kNN_proj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unwrapping denoised
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_mod1_denoised_grid = reshape(f_mod1_denoised,[m,m]);

f_unwrapped_grid = unwrap_2D(f_mod1_denoised_grid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unwrapping noisy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_unwrapped_grid = unwrap_2D(y_grid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unwrapping clean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unwrapped_clean = unwrap_2D(y_clean);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
surf(X,Y,unwrapped_clean)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_unwrapped_clean.png')

figure;
surf(X,Y,f_clean_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_clean.png')

figure;
surf(X,Y,y_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_noisy_mod1.png')

figure;
surf(X,Y,f_unwrapped_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_denoised_unwrapped.png')

figure;
surf(X,Y,y_unwrapped_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_noisy_unwrapped.png')

figure;
surf(X,Y,y_clean)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_mod1clean_unwrapped.png')

%% Contour plots

figure;
contour(X,Y,f_clean_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_clean_contour.png')

figure;
contour(X,Y,f_unwrapped_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_denoised_unwrapped_contour.png')

figure;
contour(X,Y,y_unwrapped_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_noisy_unwrapped_contour.png')


figure;
contour(X,Y,f_mod1_denoised_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_denoised_mod1_contour.png')

figure;
contour(X,Y,y_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_noisy_mod1_contour.png')

figure;
contour(X,Y,f_mod1_clean_grid)
xlabel('','Interpreter','latex', 'FontSize', 25)
ylabel('','Interpreter','latex', 'FontSize', 25)
saveas(gcf,'figures/2D/Vesuvius_clean_mod1_contour.png')

