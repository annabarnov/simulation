%% 
clear all; close all;
col={'b', 'g', 'r', 'k', 'm'};

%% compare theory with monte carlo
t60arr = [0 0.2 0.5 0.8]; 
theta_desired_arr = [0 30 60 90];
fig = 0;
for ttt = 1:length(t60arr)
    fig = fig + 1;
    f = figure(fig);
   
    for hhh = 1:length(theta_desired_arr)
        
        t60             = t60arr(ttt);
        theta_desired   = theta_desired_arr(hhh);
        M_arr           = 1:9;
        nstart          = 20;
        nfft            = 1024;
        neff = nfft/2 + 1;
        store_str = ['B_power_' num2str(t60*1000) '_' num2str(theta_desired) '_' num2str(M_arr(1)) '_' num2str(length(M_arr))];
        store_str_theory = ['B_power_theory_' num2str(t60*1000) '_' num2str(theta_desired) '_' num2str(M_arr(1)) '_' num2str(length(M_arr))];

        dtheta = 2;
        theta_arr           = (0:dtheta:360)*pi/180; theta_arr = theta_arr(1:end-1);
        % load theory
        load([store_str_theory '.mat'])
        % load monte carlo
        load([store_str '.mat'])

        figure(fig);
%         suptitle(['reverbaration = ' num2str(t60) 'sec'])
        subplot(2,2,hhh)
        plot(theta_arr*180/pi - 180, fftshift(squeeze(10*log10(real(mean(B_power_theory(nstart:neff,:),1))))))
        hold on;
        plot(theta_arr*180/pi - 180, fftshift(squeeze(10*log10(real(mean(B_power(nstart:neff,:),1))))),'r')
        title(['constrined \theta = ' num2str(theta_desired)])
        ylim([-15, 5])
        legend('Theoretical Model', 'Monte Carlo')
    end
    
end
%% compare reverbaration times

t60arr = [0 0.2 0.5 0.8]; 
theta_desired_arr = [0 30 60 90];
fig = 0;
for ttt = 1:length(theta_desired_arr)
    fig = fig + 1;
   
   
    for hhh = 1:length(t60arr)
        
        t60             = t60arr(hhh);
        theta_desired   = theta_desired_arr(ttt);
        M_arr           = 1:9;
        nstart          = 20;
        nfft            = 1024;
        neff = nfft/2 + 1;
        store_str = ['B_power_' num2str(t60*1000) '_' num2str(theta_desired) '_' num2str(M_arr(1)) '_' num2str(length(M_arr))];
        store_str_theory = ['B_power_theory_' num2str(t60*1000) '_' num2str(theta_desired) '_' num2str(M_arr(1)) '_' num2str(length(M_arr))];

        dtheta = 2;
        theta_arr           = (0:dtheta:360)*pi/180; theta_arr = theta_arr(1:end-1);
        % load theory
        load([store_str_theory '.mat'])
        % load monte carlo
        load([store_str '.mat'])

        figure(fig);
        hold on;
        plot(theta_arr*180/pi - 180, fftshift(squeeze(10*log10(real(mean(B_power(nstart:neff,:),1))))),col{hhh})
        
%         title(['T60 = ' num2str(t60) 'sec'])
     
    end
    figure(fig)
    str = strtrim(cellstr(num2str(t60arr.')));
    legend(str{:})
    title(['Constrained \theta =  ' num2str(theta_desired) 'deg'])
end

%% compare 3 mics aliasing

t60arr = [0 0.8];%[0 0.2 0.5 0.8]; 
theta_desired_arr = [0 30];%[0 30 60 90];
M_arr_cell = {[1:9], [1:2:9]};
[spcing1, spcing2] = M_arr_cell{:};
spacing(1) = -1*(spcing1(1) - spcing1(2))*0.5;
spacing(2) = -1*(spcing2(1) - spcing2(2))*0.5;
fig = 0;
for ttt = 1:length(t60arr)
    fig = fig + 1;
    f = figure(fig);
   
    for hhh = 1:length(theta_desired_arr)
        
        for mmm = 1:length(M_arr_cell)
        
        t60             = t60arr(ttt);
        theta_desired   = theta_desired_arr(hhh);
        M_arr           = M_arr_cell{mmm};
        nstart          = 20;
        nfft            = 1024;
        neff = nfft/2 + 1;
        store_str = ['B_power_' num2str(t60*1000) '_' num2str(theta_desired) '_' num2str(M_arr(1)) '_' num2str(length(M_arr))];
        store_str_theory = ['B_power_theory_' num2str(t60*1000) '_' num2str(theta_desired) '_' num2str(M_arr(1)) '_' num2str(length(M_arr))];

        dtheta = 2;
        theta_arr           = (0:dtheta:360)*pi/180; theta_arr = theta_arr(1:end-1);
        % load theory
        load([store_str_theory '.mat'])
        % load monte carlo
        load([store_str '.mat'])

        figure(fig);
       
        
        subplot(2,1,hhh)
        hold on
        fff = 128;
        plot(theta_arr*180/pi - 180, fftshift(squeeze(10*log10(real(mean(B_power(fff:fff,:),1))))), col{mmm})
        title(['Constrained \theta = ' num2str(theta_desired)])
        legend(['interspacing of ' abs(num2str(spacing(1))) ' cm'], ['interspacing of ' abs(num2str(spacing(2))) ' cm'])
        end
        
    end
    
end