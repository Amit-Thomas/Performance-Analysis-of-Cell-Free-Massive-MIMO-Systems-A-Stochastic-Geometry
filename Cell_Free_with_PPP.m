function [SINR_Deterministic, SINR_Simulation]=cell_free_with_PPP(user_count,antenna_count ,ap_density,rho_tr ,rho_d,alpha,tau_tr ,NoisePower,ap_iterations,SINR_iterations)

    % Creating Pilot Sequences and Distributing Amongst Users
    available_pilots=orth(rand(tau_tr, tau_tr));
    pilots = zeros(tau_tr,user_count);
    for k=1:user_count
        pilot_number=mod(k-1,tau_tr)+1;
        pilots(:,k)=available_pilots(:,pilot_number);
    end
    SINR_Simulation = zeros(ap_iterations*user_count,1);
    SINR_Deterministic = zeros(ap_iterations*user_count,1);
    
    for i = 1:ap_iterations

        % Generating Different AP and User Positions
        ap_count = poissrnd(ap_density);%We assume the area is a 1kmx1km square
        ap_pos = rand(ap_count, 2)*1000;
        user_pos = rand(user_count, 2)*1000;
        lmk = zeros(ap_count, user_count);

        % Calculating Wrapped Minimum Distance between APs and Users
        for m = 1:ap_count
           for k = 1:user_count
               min_distance=min(vecnorm(abs(ap_pos(m, :) - user_pos(k, :))-[0,0;1000,0;0,1000;1000,1000],2,2));%APs are wrapped around due to finite area
               lmk(m, k) = min((min_distance)^(-alpha), 1);
           end
        end
        
        % Defining covariance matrices of channel estimate vector and
        % corresponding inverse
        dmk =lmk*abs(pilots'*pilots).^2+1/(tau_tr*rho_tr/NoisePower);
    
        phi_mk = lmk.^2./dmk;
        c_mk = phi_mk.^(-1);
        SINR_Num = zeros(user_count, 1);
        SINR_Den = zeros(user_count, 1);
        
        for j = 1:SINR_iterations

            %% Uplink Training
            %Channel Generation
            gmk = normrnd(0, sqrt(0.5), ap_count*antenna_count, user_count) + 1j*normrnd(0, sqrt(0.5), ap_count*antenna_count, user_count);
            hmk=kron(sqrt(lmk), ones(antenna_count, 1)).*gmk;
            y_tr_m=zeros(ap_count*antenna_count,tau_tr); %received sigmal at the APs
    
            for k =1:user_count
               y_tr_m=y_tr_m+sqrt(tau_tr*rho_tr/NoisePower)*hmk(:,k)*pilots(:,k)'; % Adding Pilots
            end         
               
            % Uplink Training Receive Signal
            y_tr_m=y_tr_m+sqrt(0.5)*(randn(ap_count*antenna_count,tau_tr)+1j*randn(ap_count*antenna_count,tau_tr)); % Adding Noise
            
            % Post Matched Filter Receive Signal
            y_mk = (1/sqrt(tau_tr*rho_tr/NoisePower))*y_tr_m*pilots; % Projection onto Pilot
        
            %LMMSE estimation of channel
            h_hat_mk =kron(lmk./dmk ,ones(antenna_count,1)).*y_mk ;
             
            mu = ap_count/(sum(c_mk, 'all')); %normalising factor         
                 
        
            %% SINR Numerator and Denominator Averages Calculation
            SINR_Num = SINR_Num + sum(conj(hmk).*kron(c_mk,ones(antenna_count,1)).*h_hat_mk,1).';
            SINR_Den = SINR_Den + sum(abs(hmk'*(kron(c_mk,ones(antenna_count,1)).*h_hat_mk)).^2,2);
        
        end
        
        % SINR calculation using simulation
        SINR_Simulation((i-1)*user_count+1:i*user_count) =  abs(SINR_Num).^2./(SINR_iterations*SINR_Den - abs(SINR_Num).^2 + SINR_iterations^2*NoisePower/(mu*rho_d));
        
        % SINR calculation using Deterministic Equivalent technique
        SINR_Deterministic((i-1)*user_count+1:i*user_count)=ap_count*antenna_count*ones(user_count,1)./(1/ap_count*sum((lmk+ap_count*antenna_count*NoisePower/rho_d)'*(dmk.*lmk.^(-2)),2)-1);
    end
end