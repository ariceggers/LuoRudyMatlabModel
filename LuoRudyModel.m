%single cell simulation
clear;

figure(1); clf

for s = 1:7
    
del_t=0.04; %intergration time step

%Initial value of membrane potential 
v=-84.0; %mV

%MODEL PARAMETERS 
g_Na=23;    %Max Na conductance (units: mO^-1.cm^-2) 
g_Ca=.3.*.09;    %Max Ca conductance 
g_K=0.282;    %Max K (time dep.) conductance 
g_K1=0.6047;  %Max K (time indep.) conductance 
g_to=.5.*(.1652);
g_Kur=.5.*(0.005+((0.05)./(1+exp(-1.*((v-15)./(13))))));

Ko = 5.4;%mM
Ki = 145;%mM
      
Nao = 140;
Nai = 18;

PRNaK = 0.01833;  % Na - K exchanger rate

KQ10=3;

R = 8.315; %(J/(mol*K)
T = 310;   %(K)  
F = 96.49; %(kC/mol where kC = J/mV);


E_Na=54.4;    %Reversal potential of Na (units: mV) 
E_K=R*T/F*log((Ko + PRNaK*Nao)/(Ki + PRNaK*Nai));      %Reversal potential of K (time dep.) 
E_K1=R*T/F*log(Ko/Ki);  %Reversal potential of K (time-indep.) 
E_Kp=E_K1;  %Reversal potential of K (plateau) 
 



%voltage gating parameters 
alpha_m=0.32*(v+47.13)./(1-exp(-0.1*(v+47.13))); 
beta_m=0.08*exp(-v/11.0); 
alpha_h=0.5*(1-sign(v+40+eps))*0.135.*exp(-1*(80+v)/6.8); 
beta_h=(3.56*exp(0.079*v)+3.1*1e5*exp(0.35*v))*0.5.*(1-sign(v+40+eps))+(1./(0.13*(1+exp(-1*(v+10.66)/11.1))))*0.5.*(1+sign(v+40+eps)); 
alpha_j=0.5*(1-sign(v+40+eps)).*(-1.2714*1e5*exp(0.2444*v)-3.474*1e-5*exp(-0.04391*v)).*(v+37.78)./(1+exp(0.311*(v+79.23))); 
beta_j=0.5*(1-sign(v+40+eps)).*0.1212.*exp(-0.01052*v)./(1+exp(-0.1378*(v+40.14)))+0.5*(1+sign(v+40+eps))*0.3.*exp(-2.535*1e-7*v)./(1+exp(-0.1*(v+32))); 
alpha_d=0.095*exp(-0.01*(v-5))./(1+exp(-0.072*(v-5))); 
beta_d=0.07*exp(-0.017*(v+44))./(1+exp(0.05*(v+44))); 
alpha_f=0.012*exp(-0.008*(v+28))./(1+exp(0.15*(v+28))); 
beta_f=0.0065*exp(-0.02*(v+30))./(1+exp(-0.2*(v+30))); 
alpha_x=0.0005*exp(0.083*(v+50))./(1+exp(0.057*(v+50))); 
beta_x=0.0013*exp(-0.06*(v+20))./(1+exp(-0.04*(v+20))); 
alpha_k1=1.02./(1+exp(0.2385*(v-E_K1-59.2915))); 
beta_k1=(0.49124*exp(0.08032*(v-E_K1+5.476))+exp(0.06175*(v-E_K1-594.31)))./(1+exp(-0.5143*(v-E_K1+4.753)));

alpha_oa=0.65.*(exp(-1.*((v+10)/(8.5)))+exp(-1.*((v-30)/(59)))).^(-1);
beta_oa=0.65.*(2.5+exp((v+82)/(17))).^(-1);

alpha_oi=(18.53+exp(((v+113.7)/(10.95)))).^(-1);
beta_oi=(35.56+exp(-1.*((v+1.26)/(7.44)))).^(-1);

alpha_ua=0.65.*(exp(-1.*((v+10)./(8.5)))+exp(-1.*((v-30)/(59)))).^(-1);
beta_ua=0.65.*(2.5+exp((v+82)/17)).^(-1);

%Initial conditions: gating variables & Ca concentration 
m=alpha_m./(alpha_m+beta_m);      %Steady state values 
h=alpha_h./(alpha_h+beta_h); 
j=alpha_j./(alpha_j+beta_j); 
d=alpha_d./(alpha_d+beta_d); 
f=alpha_f./(alpha_f+beta_f); 
x=alpha_x./(alpha_x+beta_x); 
oi=(1+exp((v+43.1)/(5.3))).^(-1);
oa=(1+exp(-1.*((v+20.47)/(17.54)))).^(-1);
ua= (1+exp((v-99.45)/27.48)).^(-1);


Ca =2e-4; %mmol/L - initial (resting) intracellular Ca concentration 
 
t=0;
for i=1:15000   %50000 iterations = 500 msec

 %Stimulation current 
 I_stim =0; 
% current stimulation of 1 msec applied at 20 msec
if i>500 & i<(525)
      I_stim=-40.0;
    else
      I_stim=0;
end

  
 E_Ca = 7.7 - 13.0287*log(Ca);   % Calcium reversal potential 
 x_i=0.5*(1+sign(v+100-eps))*2.837.*(exp(0.04*(v+77))-1)./((v+77).*exp(0.04*(v+35)))+0.5*(1-sign(v+100-eps)); 
 


dX = zeros(2,1);
if s == 7 
    D = (2.55.*10.^(-5)).*(s)./1.4;
else
    D = (2.55.*10.^(-5)).*(s-1)./5;
end
kr = 0;
ka = 26000;
ki = 0;
lr = .05;
la = 2.7;
li = 6.*10.^(-5);

alpha_hd=0.5*(1-sign((v+40)+40+eps))*0.135.*exp(-1*(80+(v+40))/6.8); 
beta_hd=(3.56*exp(0.079*(v+40))+3.1*1e5*exp(0.35*(v+40)))*0.5.*(1-sign((v+40)+40+eps))+(1./(0.13*(1+exp(-1*((v+40)+10.66)/11.1))))*0.5.*(1+sign((v+40)+40+eps));

if i == 1
    X = [0 0];
else
    [T,X] = ode15s(@(t,X) [alpha_hd.*(X(2) - X(1)) + (kr.*(h-((m.^3)*h)) + ka.*((m.^3)*h)).*D - beta_hd.*((X(1) - ((m.^3).*X(1))) + ((m.^3).*X(1))) - (lr.*(X(1) - ((m.^3).*X(1))) + la.*((m.^3).*X(1))); (kr.*(h-((m.^3)*h)) + ka.*((m.^3)*h) + ki.*(1 - X(2) - h)).*D - (lr.*(X(1) - ((m.^3).*X(1))) + la.*((m.^3).*X(1)) + li.*(X(2) - X(1)))],[0 t],[0 0]);
end

[Max slsls] = max(X(:,2));

 %Ionic currents 
 I_Na = (g_Na*(m.^3)*h*j*(v-E_Na)).*(1-Max);;     %Fast sodium current (inward) 
 I_Ca = g_Ca*d*f*(v-E_Ca); %Slow calcium current (inward)
 I_K = g_K*x*x_i*(v-E_K);    %Time-dependent potassium current (outward) 
 I_K1 = g_K1*(alpha_k1/(alpha_k1+beta_k1))*(v-E_K1);    %Time-independent potassium current (outward) 
 I_Kp = 0.0183*(v-E_Kp)/(1+exp((7.488-v)/5.98));    %Plateau potassium current (outward) 
 I_b = 0.03921*(v+59.87);    %Background current (outward)
 I_to=g_to.*(oa.^3).*oi.*(v-E_K);
 I_Kur=g_Kur.*(ua.^3).*ua.*(v-E_K);
 I_ion = I_Na + I_Ca + I_K + I_K1 + I_Kp + I_b + I_to + I_Kur;   %Total ionic current

  
    %TAU_X AND X_INFINITY (variables for the gate ODEs) 
    tau_m=1./(alpha_m+beta_m);m_inf=alpha_m./(alpha_m+beta_m); 
    tau_h=1./(alpha_h+beta_h);h_inf=alpha_h./(alpha_h+beta_h); 
    tau_j=1./(alpha_j+beta_j);j_inf=alpha_j./(alpha_j+beta_j); 
    tau_d=1./(alpha_d+beta_d);d_inf=alpha_d./(alpha_d+beta_d); 
    tau_f=1./(alpha_f+beta_f);f_inf=alpha_f./(alpha_f+beta_f); 
    tau_x=1./(alpha_x+beta_x);x_inf=alpha_x./(alpha_x+beta_x);
    tau_oa=((alpha_oa+beta_oa).^(-1))./KQ10; oa_inf=(1+exp(-1.*((v+20.47)/(17.54)))).^(-1);
    tau_oi=((alpha_oi+beta_oi).^(-1))./KQ10; oi_inf=(1+exp((v+43.1)/(5.3))).^(-1);
    tau_ua=((alpha_ua+beta_ua).^(-1))./KQ10;   ua_inf=(1+exp((v-99.45)/27.48)).^(-1);
 del_v = -I_ion-I_stim; 
 
 %del_m = (m_inf-m)./tau_m;   %alpha_m*(1-m)-beta_m*m; 
 %del_h = (h_inf-h)./tau_h;   %alpha_h*(1-h)-beta_h*h; 
 %del_j = (j_inf-j)./tau_j;   %alpha_j*(1-j)-beta_j*j; 
 
 m = m_inf - (m_inf-m)*exp(-del_t/tau_m);     % the hybrid method
 h = h_inf - (h_inf-h)*exp(-del_t/tau_h);
 j = j_inf - (j_inf-j)*exp(-del_t/tau_j); 
 
 oi= oi_inf - (oi_inf-oi)*exp(-del_t/tau_oi); 
 oa= oa_inf - (oa_inf-oa)*exp(-del_t/tau_oa); 
 
 ua= ua_inf - (ua_inf-ua)*exp(-del_t/tau_ua);
 
 
 del_d = (d_inf-d)./tau_d;   %alpha_d*(1-d)-beta_d*d; 
 del_f = (f_inf-f)./tau_f;   %alpha_f*(1-f)-beta_f*f; 
 del_x = (x_inf-x)./tau_x;   %alpha_x*(1-x)-beta_x*x; 
 del_Ca = -1e-4*I_Ca+0.07*(1e-4-Ca); 
  
   %UPDATE (forward euler) 
    %v=v+del_t*del_v;m=m+del_t*del_m;h=h+del_t*del_h;j=j+del_t*del_j; 
    v=v+del_t*del_v;
    d=d+del_t*del_d;
    f=f+del_t*del_f;
    x=x+del_t*del_x;
    Ca=Ca+del_t*del_Ca;

  
    %updating voltage gating variables 
    alpha_m = 0.32*(v+47.13)./(1-exp(-0.1*(v+47.13))); 
    beta_m = 0.08*exp(-v/11.0); 
    alpha_h = 0.5*(1-sign(v+40+eps))*0.135.*exp(-1.0*(80+v)/6.8); 
    beta_h=(3.56*exp(0.079*v)+3.1*1e5*exp(0.35*v))*0.5.*(1-sign(v+40+eps))+ (1./(0.13*(1+exp(-1*(v+10.66)/11.1))))*0.5.*(1+sign(v+40+eps)); 
    alpha_j = 0.5*(1-sign(v+40+eps)).*(-1.2714*1e5*exp(0.2444*v) - 3.474*1e-5*exp(-0.04391*v)).*(v+37.78)./(1+exp(0.311*(v+79.23))); 
    beta_j = 0.5*(1-sign(v+40+eps)).*0.1212.*exp(-0.01052*v)./(1+exp(-0.1378*(v+40.14))) + 0.5*(1+sign(v+40+eps))*0.3.*exp(-2.535*1e-7*v)./(1+exp(-0.1*(v+32))); 
    alpha_d = 0.095*exp(-0.01*(v-5))./(1+exp(-0.072*(v-5))); 
    beta_d = 0.07*exp(-0.017*(v+44))./(1+exp(0.05*(v+44))); 
    alpha_f = 0.012*exp(-0.008*(v+28))./(1+exp(0.15*(v+28))); 
    beta_f = 0.0065*exp(-0.02*(v+30))./(1+exp(-0.2*(v+30))); 
    alpha_x = 0.0005*exp(0.083*(v+50))./(1+exp(0.057*(v+50))); 
    beta_x = 0.0013*exp(-0.06*(v+20))./(1+exp(-0.04*(v+20))); 
    alpha_k1 = 1.02./(1+exp(0.2385*(v-E_K1-59.2915))); 
    beta_k1 = (0.49124*exp(0.08032*(v-E_K1+5.476)) + exp(0.06175*(v-E_K1-594.31)))./(1+exp(-0.5143*(v-E_K1+4.753)));
    alpha_oa=0.65.*(exp(-1.*((v+10)/(8.5)))+exp(-1.*((v-30)/(59)))).^(-1);
    beta_oa=0.65.*(2.5+exp((v+82)/(17))).^(-1);
    alpha_oi=(18.53+exp(((v+113.7)/(10.95)))).^(-1);
    beta_oi=(35.56+exp(-1.*((v+1.26)/(7.44)))).^(-1);    
    alpha_ua=0.65.*(exp(-1.*((v+10)./(8.5)))+exp(-1.*((v-30)/(59)))).^(-1);
    beta_ua=0.65.*(2.5+exp((v+82)/17)).^(-1);

    
 %Recording the time series of the voltage and currents 
    % To plot membrane potential, "plot(rec_v)" 
    t=t+del_t;
     rec_I_Na(i)=I_Na; 
     rec_I_Ca(i)=I_Ca; 
     rec_I_K(i)=I_K; 
     rec_I_K1(i)=I_K1; 
     rec_I_Kp(i)=I_Kp; 
     rec_I_b(i)=I_b; 
     rec_I_ion(i)=I_ion; 
     rec_I_to(i)=I_to;
     rec_I_Kur(i)=I_Kur;
     rec_Ca(i)=Ca; 
     rec_v(i)=v; 
     rec_t(i)=t;

end

figure(1);
hold on;
t=linspace(0,150,15000);
if s == 1
    plot(rec_v,'k-')
elseif s == 2
    plot(rec_v,'b--')
elseif s == 3
    plot(rec_v,'g--')
elseif s == 4
    plot(rec_v,'r--')
elseif s == 5
    plot(rec_v,'m--')
elseif s == 6
    plot(rec_v,'k--')
elseif s == 7
    plot(rec_v,'c--')
end
title('Membrane Potential (mV) vs. Time (ms) For Varying Quinidine Concentrations')
xlabel('Time (ms)')
ylabel('Membrane Voltage (mV)')
legend('0x','.2x','.4x','.6x','.8x','1x','5x')

end

hold off
grid on

