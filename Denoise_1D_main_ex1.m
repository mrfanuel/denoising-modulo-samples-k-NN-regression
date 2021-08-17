%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for the modulo denoising via kNN and unwrapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

addpath(genpath('lib'))  

% choose the method below
method = 'kNN';%'TRS';%'UCQP';
disp(method)
disp('Started')

%% Parameters
sigma = 0.12; % noise level 
n_MC = 50; % number of Monte-Carlo runs

range_n = 100:50:1000; % range of n values
c = 0.04;%0.07;
if strcmp(method,'kNN')
    C_kNN = 0.09; %0.07;%
    k = ceil(C_kNN*(range_n.^(2/3)).*(log(range_n).^(1/3))); %% number of neigbours
elseif strcmp(method,'UCQP')
    C_UCQP =  c;
    lambda = C_UCQP*(range_n.^(10/3)).^(1/4);
elseif strcmp(method,'TRS')
    C_TRS =  c;
    lambda = C_TRS*(range_n.^(10/3)).^(1/4);
end




%% Initialization mean errors and std's
err_wrap_around = zeros(size(range_n)); std_err_wrap_around = zeros(size(range_n));
err_wrap_around_noisy = zeros(size(range_n));std_err_wrap_around_noisy = zeros(size(range_n));
err_unwrapped = zeros(size(range_n)); std_err_unwrapped = zeros(size(range_n));
err_unwrapped_noisy = zeros(size(range_n)); std_err_unwrapped_noisy = zeros(size(range_n));


err_wrap_around_temp = ones(n_MC,1);
err_wrap_around_noisy_temp = ones(n_MC,1);
err_unwrapped_temp = ones(n_MC,1);
err_unwrapped_noisy_temp = ones(n_MC,1);


%% Ground truth
ff = @(x) sin(4*pi*x);

a = 0;b = 1; 

for index = 1:length(range_n)

    n = range_n(index);

    x = (a:((b-a)/(n-1)):b)'; % nx1 vector  
    f_clean = ff(x); 
    f_mod1_clean = mod(f_clean,1); 
    
        
    if strcmp(method,'kNN')
        nb_nb = k(index);
    else
        %% For graph-based methods
        r = 1;
        % Form graph G with "connectivity radius" = r and its Laplacian
        A = zeros(n);
        for i=1:n
            for j=1:n
                if i~=j & abs(i-j) <= r
                    A(i,j)=1;
                end
            end
        end
        d = A*ones(n,1);
        L = diag(d) - A;
    end

    for iter= 1:n_MC

        f_noise = f_clean + sigma * randn(n,1); 
        y = mod(f_noise,1);  
        z = exp(1i*2*pi*y);  % nx1 vector
        
        %% Denoising
        if strcmp(method,'kNN')
            gest_kNN = kNN_denoise(z,x,nb_nb);
            gest_kNN_proj = project_manifold(gest_kNN);
            f_mod1_denoised = extract_modulo(gest_kNN_proj);
            
        elseif strcmp(method,'UCQP')
            reg_param = lambda(index);
            gest_ucqp = UCQP_denoise(z,L,reg_param,n);
            gest_ucqp_proj = project_manifold(gest_ucqp);
            f_mod1_denoised = extract_modulo(gest_ucqp_proj);      
        elseif strcmp(method,'TRS')
            reg_param = lambda(index);
            gest_trs = TRS_denoise(z,L,reg_param,n);
            gest_trs_proj = project_manifold(gest_trs);
            f_mod1_denoised = extract_modulo(gest_trs_proj); 
        end
     
        err_wrap_around_temp(iter) = MS_wrap_around_error(f_mod1_denoised, f_mod1_clean);
        err_wrap_around_noisy_temp(iter) = MS_wrap_around_error(y, f_mod1_clean);

        %% Unwrapping
        f_unwrapped = unwrap_1D(f_mod1_denoised);

        err_unwrapped_temp(iter) = MS_error(f_unwrapped, f_clean);

        y_noisy_unwrapped = unwrap_1D(y);
        err_unwrapped_noisy_temp(iter) = MS_error(y_noisy_unwrapped, f_clean);

    end


    err_wrap_around(index) = mean(err_wrap_around_temp);
    std_err_wrap_around(index) = std(err_wrap_around_temp);

    err_wrap_around_noisy(index) = mean(err_wrap_around_noisy_temp);
    std_err_wrap_around_noisy(index) = std(err_wrap_around_noisy_temp);

    err_unwrapped(index) = mean(err_unwrapped_temp);
    std_err_unwrapped(index) = std(err_unwrapped_temp);

    err_unwrapped_noisy(index) = mean(err_unwrapped_noisy_temp);
    std_err_unwrapped_noisy(index) = std(err_unwrapped_noisy_temp);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot the results
subplot(6,1,1);
plot(x,f_clean); % Unwrapped Clean mod 1 samples
title('Clean ground truth','Interpreter','latex') 

subplot(6,1,2);
plot(x,f_mod1_clean); % Clean mod 1 samples
title('Clean $\bmod 1$','Interpreter','latex') 

subplot(6,1,3);
plot(x,y); % Noisy mod 1 samples
title('Noisy $\bmod 1$','Interpreter','latex') 

subplot(6,1,4);
plot(x,f_mod1_denoised); % denoised mod 1 samples
title('Denoised $\bmod 1$','Interpreter','latex') 

subplot(6,1,5);
plot(x,f_unwrapped); % Unwrapped Clean mod 1 samples
title('Unwrapped denoised','Interpreter','latex') 

subplot(6,1,6);
plot(x,y_noisy_unwrapped); % Unwrapped Clean mod 1 samples
title('Unwrapped noisy','Interpreter','latex') 

place = strcat(strcat('/ex1_Comparison_Denoising_UnWrapping_1D_',method),'.png');
folder = strcat('figures/',method);
place = strcat(folder,place);
saveas(gcf,place)



%% Plot the results for the paper

%%%%%%%%%%%%%%% MSE Wrap around error %%%%%%%%%%%%%%% 
close all;
figure;
plot(range_n, err_wrap_around, '-bo', 'MarkerSize', 4, 'markerfacecolor','b');hold on;
errorbar(range_n, err_wrap_around, std_err_wrap_around,'b' )
xlim([95,1005])

hold on;
plot(range_n, err_wrap_around_noisy, '-k<', 'MarkerSize', 4, 'markerfacecolor','k')
errorbar(range_n, err_wrap_around_noisy, std_err_wrap_around_noisy,'k' )

xlabel('$n$','Interpreter','latex', 'FontSize', 25)
ylabel('Wrap around RMSE','Interpreter','latex', 'FontSize', 25)

place = strcat(strcat('/ex1_paper_err_mod1_',method),'.png');
folder = strcat('figures/',method);
place = strcat(folder,place);
saveas(gcf,place)

%%%%%%%%%%%%%%% MSE %%%%%%%%%%%%%%% 

close all;
figure;
plot(range_n, err_unwrapped, '-bo', 'MarkerSize', 4, 'markerfacecolor','b');hold on;
errorbar(range_n, err_unwrapped, std_err_unwrapped ,'b')
xlim([95,1005])

hold on;
plot(range_n, err_unwrapped_noisy, '-r<', 'MarkerSize', 4, 'markerfacecolor','r');hold on;
errorbar(range_n, err_unwrapped_noisy, std_err_unwrapped_noisy,'r' )

xlabel('$n$','Interpreter','latex', 'FontSize', 25)
ylabel('RMSE','Interpreter','latex', 'FontSize', 25)
place = strcat(strcat('/ex1_paper_err_unwrapped_',method),'.png');
folder = strcat('figures/',method);
place = strcat(folder,place);
saveas(gcf,place)

disp('ended')

