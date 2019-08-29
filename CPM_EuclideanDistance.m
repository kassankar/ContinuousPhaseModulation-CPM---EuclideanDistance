clc;
clear;
close all;

%% Pulse shape & Variable ini
MAIN          = CPM_Main_Functions_EuclideanDistance;
pulse         = 1;   % 1 -> lorentzian pulse
                     % 2 -> GMSK pulse BT = 0.3
                     % 3 -> LRC pulse
                     % 4 -> LREC pulse                        
L             = 4;  % Pulse length            
                     % 1  -> Full response
                     % >1 -> Partial response                   
Fs            = 128; % Sampling frequency 
Ts            = 1/Fs;% Sampling Time
M             = 4;   % M_array symbols used
                     % 2 -> Binary
h_min  = 0.02;       % hmin should be taked higher than 0 (it can be qual to 0 for the Upper Bound), 
                     % so the calculation of dmin don't take all combination for h =0; (for h=0 the simulation will take all the ram)
h_max  = 2.5;
deltah = 0.01;
%% Main code

width = 0.37;        % This variable is used for Lorentzian Pulse only. (Not be used for pulse # 1)
[g_t,q_t] = MAIN.CREATECPMPULSE(pulse,L,width,Fs); % Function return the CPM pulse and phase.
                                                   % g_t = g(t) is the CPM pulse shape.
                                                   % q_t -> is the phase, integral of g_t.
                                                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euclidean Distance Upper Bound -> (Book: Digital Phase Modulation page: 71 - Paper: CPM--Part II: Partial Response Page: 212-213)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H      = h_min:deltah:h_max;
gamma0 = 2*(M-1);


for i = 1:M
    gamma(i) = 2*(M-i);                    % Create gamma_i (see page 212 from Paper)
end

%%%%%%%%%%%%%%%%%%%
% For L > 1
%%%%%%%%%%%%%%%%%%%
if (L ~= 1)
    gamma_pn  = [gamma -gamma(1:end-1)];   % take all positive and negative values from gamma 
    v         = sort(gamma_pn);            
    I         = L;                         % I is used to increase the number of observation symbol (To tighten the Upper Bound)
    if (I>5)
    I=5;                                   % I=5 is good enough.
    end
    N         = L+I;                       % Number of observation symbols.
    [comb,~]  = MAIN.PERMUREPET(v,I);      % PERMUREPET is a function used from the Main_Function file to create all possible combination of gamma_i numbers.
    gamma     = nonzeros(gamma);           % Remove zeros from gamma_i, cz it'll be used in the next step to create all combinations that start with different gamma_0.
    for i = 1:length(gamma)
    gamm0_cb          = gamma(i).*ones(length(comb),1);
    comb_cell(1,i)    = {[gamm0_cb comb]};    % always start with gamma0, with M>2 gamma0 take more than one value
    end
comb          = cat(1,cell2mat(comb_cell.')); % concatenate all the values in the comb_cell in one matrix comb to start with all the combinations of gamma0
j = 0;
comb_s0 = [];

for i = 1:length(comb)
    if (sum(comb(i,1:end))==0)
        j=j+1;
        comb_s0(j,:) = comb(i,1:end);         % Remove all combination with sum #0 (Book page 74)
    end
end

add_zero = zeros(length(comb_s0(1:end,1)),N-(I+1)); % add zeros at the end for t>L+1 we have zeros
comb_s0  = [comb_s0 add_zero];
for j= 1:length(comb_s0(1:end,1))
    for i= 1:length(H)
        
    gamma_s       = upsample (comb_s0(j,1:end),Fs);
    t_seq         = 0:Ts:length(comb_s0(j,1:end));           
    S_N           = conv(gamma_s,g_t);
    S_N           = S_N(1:length(t_seq));
    Phi_N         = cumsum(S_N)*Ts;
    Phasee        = 2*pi*H(i)*Phi_N(1:end);
    dB(j,i)       = log2(M)*(N-Ts*sum(cos(Phasee))); 
    end
end                                               % Calculate the Upper Bound distance for all combinations created, 
                                                  % using the formula in the Book page 90 (same formula will be used after to calculate the minimum distance)

%%%%%%%%%%%%%%%%%%%
% For L = 1
%%%%%%%%%%%%%%%%%%%
else
    for i = 1:M
    gamma(i) = 2*(M-i);
    end 
gamma = nonzeros(gamma);
for j = 1:length(gamma)
for i = 1:length(H)
    gamma_s       = upsample ([gamma(j) -gamma(j)],Fs);
    t_seq         = 0:Ts:2;           
    S_N           = conv(gamma_s,g_t);
    S_N           = S_N(1:length(t_seq));
    Phi_N         = cumsum(S_N)*Ts;
    Phasee        = 2*pi*H(i)*Phi_N(1:end);
    dB(j,i)       = log2(M)*(2-Ts*sum(cos(Phasee)));
end
end
end
dB_min = min(dB,[],1);                         % Take the minimum dB from all combination created, see Figure 3.11 Book page 75.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum Euclidien Distnace -> (Book: Digital Phase Modulation page: 463-464 - Paper: CPM--Part II: Partial Response Page: 215)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmax                 = 10;                       % Maximum number of observation symbols
dmin                 = 10^5*ones(1,length(H));  
gamma_0              = 2;
gamma_1              = 0;
gamma_N_1            = 0;
d1                   = zeros(1,length(H));
dN                   = 0;
i = 1;
j = 1;
N                    = 1;
scratch_pad_old      = zeros(1,N);
scratch_pad_New      = zeros(1,N);
while(h_min<h_max  ) %% add deltah for h go to 1.5
    %%%%%%%%%%%%%%%%%%%Betta%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma_s       = upsample (gamma_0,Fs);
    t_seq         = 0:Ts:length(gamma_0)+Ts;           
    S_N           = conv(gamma_s,g_t);
    S_N           = S_N(1:length(t_seq));
    Phi_N         = cumsum(S_N)*Ts;
    Phasee        = 2*pi*h_min*Phi_N(1:end);
    d1(i)         = log2(M)*(1-Ts*sum(cos(Phasee))); 
    dmin(i)       = min(dmin(i),d1(i));
    
    if(d1(i)>dB_min(i))
    else
        scratch_pad_old(j,1:end) = [gamma_0];
        j                        = j+1;
    end
    
    if(gamma_0==2*(M-1))
            sprintf('N is %d, h is %f, dmin is %f \n',N,h_min,dmin(i))
            if(Nmax==1)
                h_min                = h_min +deltah;
                i                    = i+1;
                scratch_pad_old      = zeros(1,N);
                j                    = 1;
                gamma_0              = 2; 
            else
    N = N +1;
    dmin(i) = 10^5;
    scratch_pad_New = scratch_pad_old;
    scratch_pad_old = zeros(1,N);
    gamma_N_1 = -2*(M-1);
    j               = 1;
    spdn_y    = 1;
    kk        = 0;
    N_loop    = 1;
    %%%%%%%%%%%%%%%%%%%alpha%%%%%%%%%%%%%%%%%%%%%%%%% 
    while N_loop > 0    
    gamma_1       = scratch_pad_New(spdn_y,1:end);
    gamma_s       = upsample ([gamma_1 gamma_N_1],Fs);
    t_seq         = 0:Ts:(length(gamma_1)+1)+Ts;       % + Ts optimization for the sum or trapz (so we don't get dmin = 0)     
    S_N           = conv(gamma_s,g_t);
    S_N           = S_N(1:length(t_seq));
    Phi_N         = cumsum(S_N)*Ts;
    Phasee        = 2*pi*h_min*Phi_N(1:end);
    dN            = log2(M)*((N)-Ts*sum(cos(Phasee))); 
    
    
    if(dN>dB_min(i))
    else
        scratch_pad_old(j,1:end) = [gamma_1 gamma_N_1];
        j                        = j+1;
    end
    
    dmin(i) = min(dmin(i),dN);
    
    if(gamma_N_1==2*(M-1))
        if(spdn_y==size(scratch_pad_New,1))
            sprintf('N is %d, h is %f, dmin is %f \n',N,h_min,dmin(i))
            if(N<Nmax)
                if(sum(scratch_pad_old==0))
%                     gamma_N_1 = -2*(M-1);
%                     break;
                end
            dmin(i) = 10^5;
            N = N + 1;
            end
            kk     = kk+1;
            j      = 1;
            spdn_y = 1;
            gamma_N_1 = -2*(M-1);
            scratch_pad_New = scratch_pad_old;
            scratch_pad_old = zeros(1,N);
            
        else
            spdn_y = spdn_y+1;
            gamma_N_1 = -2*(M-1);
        end
    else
        gamma_N_1 = gamma_N_1+2;
    end
    
    if(kk==(Nmax-1))
        N_loop = 0;
    else
        
    end
    
    end
    h_min                = h_min +deltah;
    i                    = i+1;
    N                    = 1;
    scratch_pad_old      = zeros(1,N);
    scratch_pad_New      = zeros(1,N);
    j                    = 1;
    gamma_0              = 2;            
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            end
    else
    gamma_0 = gamma_0 +2;
    end
end


dmin    = [0 dmin];                             % Correction to start the plot from zero.
dB_min  = [0 dB_min];
H       = [0 H];

%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%
figure(1)

hold on;
plot(H,dB_min,'k--','LineWidth',1);
plot(H,dmin,'k','LineWidth',2);
xlabel('modulation index (h)','FontName','Arial','FontSize',12)
ylabel('d^2(h)','FontName','Arial','FontSize',12)
legend('d_{B}^2(h)','dmin^2(h)','FontName','Arial','FontSize',12)
M_string = num2str(M);
title(strcat('M-array =',M_string));
box on;
set(gca,'GridAlpha',0.35);
