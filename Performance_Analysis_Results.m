
%% default system parameters
user_count = 10;
antenna_count = 5;
ap_density = 40; %should be an integer
W_c = 20*10^6;  % Bandwidth
f_0 = 1.9*10^9; % Carrier Frequency
rho_tr = 0.1;   % Transmit Power per Pilot Symbol
rho_d = 0.2;    % Downlink Transmit Power
alpha = 3.5;    % Path Loss Exponent
B_c = 200*10^3; % Coherence Bandwidth
T_c = 10^(-3);  % Coherence Time
tau_c=B_c*T_c;
tau_tr = 10;    % Uplink Training duration (Samples)
k_b = 1.381*10^(-23);   % Boltzmann Constant
T_0 = 290;      % Noise Temperature
N_F = 10^(0.9); % Noise Figure
NoisePower= N_F*T_0*W_c*k_b;
ap_iterations = 1000;
SINR_iterations = 1000;


%%  Coverage Probability with various AP densities
ap_density_range=[40,60,100,120];
figure(); hold on;
title('Coverage Probability vs Target SINR with various AP densities')
for c1=1:length(ap_density_range)
    ap_density=ap_density_range(c1);
    [SINR_Deterministic,SINR_Simulation]=cell_free_with_PPP(user_count,antenna_count ,ap_density,rho_tr ,rho_d,alpha,tau_tr ,NoisePower,ap_iterations,SINR_iterations);
    
    P_c_Deterministic_vec=[];
    P_c_Simulation_vec=[];
    for T=10.^((0:1:15)/10)
        P_c_Deterministic_vec=[P_c_Deterministic_vec,sum((SINR_Deterministic >T))/length(SINR_Deterministic)];
        P_c_Simulation_vec=[P_c_Simulation_vec,sum((SINR_Simulation >T))/length(SINR_Simulation)];
    end
    plot(0:1:15,P_c_Deterministic_vec,'LineStyle','--','LineWidth',2)
    plot(0:1:15,P_c_Simulation_vec,'LineWidth',2)
end
legend('\lambda=40, deterministic','\lambda=40, simulation','\lambda=60, deterministic','\lambda=60, simulation','\lambda=100, deterministic','\lambda=100, simulation','\lambda=120, deterministic','\lambda=120, simulation');
xlabel('Target SINR(T) (in dB)')
ylabel('coverage probability')
hold off

ap_density=40;%default

%%  average downlink achievable SE vs number of users
user_count_range=10:30;
tau_tr_range=[10,20];
figure(); 

title('average downlink achievable SE vs number of users various \tau_{tr}')
ave_SE_simulation=zeros(length(user_count_range),length(tau_tr_range));
ave_SE_deterministic=zeros(length(user_count_range),length(tau_tr_range));
for c1=1:length(user_count_range)
    user_count=user_count_range(c1);
    for c2=1:length(tau_tr_range)
        tau_tr=tau_tr_range(c2);
        [SINR_Deterministic,SINR_Simulation]=cell_free_with_PPP(user_count,antenna_count ,ap_density,rho_tr ,rho_d,alpha,tau_tr ,NoisePower,ap_iterations,SINR_iterations);
        ave_SE_simulation(c1,c2)=mean(log2(1+SINR_Simulation));
        ave_SE_deterministic(c1,c2)=mean(log2(1+SINR_Deterministic));
    end
end
hold on;
plot(user_count_range,ave_SE_deterministic(:,1),'LineStyle','--','LineWidth',2)
plot(user_count_range,ave_SE_simulation(:,1),'LineWidth',2)
plot(user_count_range,ave_SE_deterministic(:,2),'LineStyle','--','LineWidth',2)
plot(user_count_range,ave_SE_simulation(:,2),'LineWidth',2)
legend('\tau_{tr}=10, deterministic','\tau_{tr}=10, simulation','\tau_{tr}=20, deterministic','\tau_{tr}=20, simulation');
xlabel('Number of users(K)')
ylabel('average dowmlink SE (bits/s/Hz)')
hold off

user_count=10;
tau_tr=10;

%%  average downlink achievable SE vs AP density
ap_density_range=40:10:110;
tau_tr_range=[10,20];
figure(); 
hold on;
title('average downlink achievable SE vs AP density various \tau_{tr}')
ave_SE_simulation=zeros(length(ap_density_range),length(tau_tr_range));
ave_SE_deterministic=zeros(length(ap_density_range),length(tau_tr_range));
for c1=1:length(ap_density_range)
    ap_density =ap_density_range(c1);
    for c2=1:length(tau_tr_range)
        tau_tr=tau_tr_range(c2);
        user_count=tau_tr;%assuming orthogonal pilots
        [SINR_Deterministic,SINR_Simulation]=cell_free_with_PPP(user_count,antenna_count ,ap_density,rho_tr ,rho_d,alpha,tau_tr ,NoisePower,ap_iterations,SINR_iterations);
        ave_SE_simulation(c1,c2)=mean(log2(1+SINR_Simulation));
        ave_SE_deterministic(c1,c2)=mean(log2(1+SINR_Deterministic));
    end
end
plot(ap_density_range,ave_SE_deterministic(:,1),'LineStyle','--','LineWidth',2)
plot(ap_density_range,ave_SE_simulation(:,1),'LineWidth',2)
plot(ap_density_range,ave_SE_deterministic(:,2),'LineStyle','--','LineWidth',2)
plot(ap_density_range,ave_SE_simulation(:,2),'LineWidth',2)
legend('\tau_{tr}=10, deterministic','\tau_{tr}=10, simulation','\tau_{tr}=20, deterministic','\tau_{tr}=20, simulation');
xlabel('AP density (\lambda)')
ylabel('average dowmlink SE (bits/s/Hz)')
hold off

ap_density=40;
tau_tr=10;
user_count=10;

%%  average downlink achievable SE vs path loss exponent (alpha)
alpha_range=2:0.25:5;
ap_density_range=[60,120];
figure(); 
hold on;
title('average downlink achievable SE vs path loss exponent various AP density')
ave_SE_simulation=zeros(length(alpha_range),length(ap_density_range));
ave_SE_deterministic=zeros(length(alpha_range),length(ap_density_range));
for c1=1:length(alpha_range)
    alpha =alpha_range(c1);
    for c2=1:length(ap_density_range)
        ap_density=ap_density_range(c2);
        [SINR_Deterministic,SINR_Simulation]=cell_free_with_PPP(user_count,antenna_count ,ap_density,rho_tr ,rho_d,alpha,tau_tr ,NoisePower,ap_iterations,SINR_iterations);
        ave_SE_simulation(c1,c2)=mean(log2(1+SINR_Simulation));
        ave_SE_deterministic(c1,c2)=mean(log2(1+SINR_Deterministic));
    end
end
plot(alpha_range,ave_SE_deterministic(:,1),'LineStyle','--','LineWidth',2)
plot(alpha_range,ave_SE_simulation(:,1),'LineWidth',2)
plot(alpha_range,ave_SE_deterministic(:,2),'LineStyle','--','LineWidth',2)
plot(alpha_range,ave_SE_simulation(:,2),'LineWidth',2)
legend('\lambda=60, deterministic','\lambda=60, simulation','\lambda=120, deterministic','\lambda=120, simulation');
xlabel('path loss exponent (\alpha)')
ylabel('average dowmlink SE (bits/s/Hz)')
hold off

ap_density=40;
alpha=3.5;


