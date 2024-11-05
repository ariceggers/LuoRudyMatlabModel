Modeling the Effect of Antiarrhythmic Drugs on Atrial Fibrillation
Alex Wilson, Josh Engel, Aric Eggers, Uriel Mandl
BME260
Abstract: Atrial fibrillation (AF) is the most common sustained arrhythmia in the United States, with an estimated 2.7-6.1 million people in the US currently living with the disease (Center for Disease Control). AF is characterized by a decrease in strength and duration of atrial action potential. With antiarrhythmic drugs being the most popular treatment option, modeling the AF action potential and how it is affected by antiarrhythmic drugs is especially clinically relevant. In this report, we use the Hodgkin-Huxley model of membrane capacitance to develop a model that reflects the pathophysiology of the AF action potential. We then incorporate a model for ion channel binding of the Class 1A drug quinidine to examine the effect of the antiarrhythmic drug on ion currents and the overall action potential. Comparison of the results of our model and published experimental data shows that the model successfully incorporates the effect of quinidine on the AF action potential. The model is a viable method for developing new antiarrhythmic drugs and provides insights for new therapeutic options.

Introduction

Atrial fibrillation (AF) is a common heart arrhythmia characterized by quivers that occur due to abnormal electrical activity within the upper chambers of the heart.  Exact causation varies from patient to patient, but most cases occur due to a gradual enlargement of the left atrium. As a result, it is more difficult for the heart to maintain the proper sinus rhythm, leading to disorganized electrical activity within the atria.  This leads to atria quivers, as well as the possibility of blood pooling and clotting.  
In order to combat this chaotic electrical activity, there are two commonly utilized methods of treatment: antiarrhythmic drugs and catheter ablation.  While catheter ablation has proven to be a highly successful and safe option to eliminate some forms of AF, its invasive nature prevents it from serving as the primary treatment option for many AF patients. As the noninvasive treatment option, antiarrhythmic drugs play a large role in regulating basal heart activity in a vast majority of patients, most often by altering the chemical kinetics of different ionic channels in cardiac tissue (Kowey, 1998). These drugs are split up into three classes, Class I (Sodium Channel Blockers), Class II (beta-blockers), Class III (Potassium Channel Blockers), and Class IV (Calcium Channel Blockers). Here, we analyze the action of a Class IA drug, quinidine, and its effects on the electrical potential of the heart.  It is important to note that quinidine is often prescribed in conjunction with certain shock therapies of the heart, thereby reducing the possibility for erratic electrical impulses and allowing for a safe elongation of an AF patients action potentials. Quinidine has demonstrated the ability to lengthen AP duration, thereby allowing the heart to return to normal rhythm. Past efforts have been made to accurately model the electrical potential of the heart of those affected with AFib, and those currently taking antiarrhythmic drugs. Courtemanche et al. in 1998 sought to further understand the ionic mechanisms that govern the human atrial action potential by modeling the basic underlying factors behind the various determinants of AFib through experimentally determined kinetic coefficients and differential equations (Courtemanche et al., 1998). This model was then extended to produce a precise model of an atrial fibrillation action potential after additional data was gathered (Courtemanche et al., 1999). The precision of this atrial fibrillation action potential model opens an opportunity to analyze how chemical alterations of ionic currents affect the atrial fibrillation action potential. We use this mathematical model to answer the question of how ion channel targeting antiarrhythmic drugs affect the atrial fibrillation action potential over time. Results from this model will give insight into the therapeutic mechanisms of antiarrhythmic drugs. In doing so, this model can help make more accurate predictions on the effects of future antiarrhythmic drugs and better treatments  for AF.


Model Development

Action Potential Model:
Using the basis of the circuit model created by Hodgkin and Huxley, we model the cell membrane as a capacitor connected in parallel with variable resistors and current source elements that represent ionic channels and electrical stimuli in cardiac muscle tissue. We model membrane electrical potential according to the first order differential equation:
dVdt=-(Iion+Istim)Cm
Where Istim represents the initiating stimulus current and Iion represents the total ionic current. Cm is the total capacitance of the cellular membrane. Iion is comprised of the sum:
Iion=INa+IK+Ito+IKur+IK1+IKp+ICa,L+Ip,Ca+INaK+INaCa+Ib,Na+Ib,Ca                  (2)
To simplify the execution of our model, we excluded a portion of these individual ionic currents and simplified some of the ionic current equations. Ib,Na and Ib,Ca represent background ion leakage currents and have maximal conductance values three to four orders of magnitude smaller than the other ionic currents; we combined the two currents into one background current to simplify calculations (Courtemanche et al., 1998). We also alter the expressions for the ionic currents from the Na+/K+ pump (INaK), the Na+/Ca2+ exchange (INaCa) and the sarcolemmal Ca2+ pump (Ip,Ca). Because we are examining single action potentials, we assume changes in ion concentrations to be negligible and these restorative pumps to be active at only a low level. We evaluate these currents using a constant ion exchange rate accordingly. Descriptions of each membrane current are listed below and definitions of all variables and constants can be found in the appendix.
Membrane Currents
	Fast Na+ current: The model implements ionic current according the expressions proposed by Luo and Rudy (1994) and Courtemanche and his colleagues (1998). 
INa=gNam3hjV-ENa
GNa represents the maximum sodium conductance, m, h, and j are the gating variables for sodium channels, and ENa is the membrane potential of sodium.
	Time-dependent K+ current: The following three potassium currents represent various phases of potassium channels allowing ion flow to repolarize the cell membrane.
IK=gKxxiV-EK
GK represents the maximum potassium conductance, the x terms are gating variables and EK is the membrane potential of potassium.
	Time-independent K+ current: 
IK1=gK1K1K1+K1(V-EK1)
GK1 represents the maximum potassium conductance, the alpha and beta terms are voltage gating parameters and EK1 is the membrane potential of potassium.
	Plateau K+ current: 
IKp=.0183(V-EKp)1+exp⁡(5.987.488-V)
EKp represents the plateau potassium membrane potential.

	Ultrarapid rectifier K+ current:
IKur=gKurua3uiV-EK
gKur=.005+.051+exp⁡(-13V-15)
Ua is the activation gating parameter, ui is the inactivation gating parameter, and EK is the time dependent potassium membrane potential. The maximal conductance for the ultra-rapid rectifier K+ current varies with voltage.
	L-Type Ca2+ channel current: 
ICa,L=gCadfV-ECa
GCa represents the maximum calcium conductance, d and f are the gating parameters, and ECa is the membrane potential of calcium.
	Transient outward current: Transient outward current is taken from and developed in Courtemanche’s 1998 publication modeling the atrial action potential.
Ito=gtooa3oiV-EK
Gto is the maximum transient outward conductance, oa is the activation gating parameter, oi is the inactivation gating parameter, and EK is the time-dependent potassium membrane potential.
	Background current: 
Ib=.03921V+59.87
All other equations governing the calculation of gating parameters can be found in the appendix.
Modeling the Effect of Drug Binding: We incorporate the model of ion channel gate binding formulated by Hondeghem and Katzung (1977). This model assumes that drugs can bind to channels in either a resting (R), activated (A), or inactivated (I) state. Conservation of mass is assumed for all mass balances. Mass balance analysis on drug binding to ion channels yields:
R+A+I+R'+A'+I'=1
R, A, and I represent the fraction of unbound channels, R’, A’, and I’ represent the fraction of bound channels, and the total amount of ion channels is set to 1. Furthermore, a mass balance on bound channels yields:
B=R'+A'+I'
where B is the fraction of bound ion channels. In this experiment, we study the effect of a sodium channel inhibitor; the R, A and I states, bound and unbound, relate to the m, h and j gating parameters for sodium channels mentioned above. The following differential equations are used as a basis for modeling drug binding:
dBdt=kRR+kAA+kII-lRR'+lAA'+lII
This represents the chemical kinetics equation for the rate of drug binding, with that rate constants for association being notated with k and the rate constants for dissociation being notated with l. Drug binding also affects action at the h gate in sodium channels according to:
dh'dt=h'I'+kRR+kAAD-h'R'+A'-lRR'+lAA'
This represents the rate of binding to sodium h gates, with the alpha and beta terms being voltage-dependent functions governing gate states, and D representing the concentration of the drug.
Numerical Methods: Calculations were performed using a combination of forward Euler method numerical integration with a time step fixed at .04 ms and the MATLAB ode15s differential equation solver, which uses the Runge-Kutta numerical integration method. All cell conditions and parameters are included in the Appendix, along with MATLAB code showing the numerical methods of the model along with more underlying equations.
Results
In Fig. 1, a normal atrial action potential is compared to an atrial action potential during fibrillation. In order to produce an action potential during fibrillation, the conductances for the Calcium(L) Channels were decreased by 70%, and the conductances for the Transient Outward Current and the Ultra-Rapid Potassium current were reduced by 50%. This is representative of the pathology of atrial fibrillation as it is propagated by non-standard conductances, and has been linked to decreases in these specific conductance values (Courtemanche et. al). 
Figure 1	
	
In Fig. 2, the effect of varying Quinidine concentrations on atrial action potentials are demonstrated. Increasing concentrations produce action potentials with longer effective refractory periods and slower phase 0 depolarization.
Figure 2
	Fig. 3 represents the effects of Quinidine on an atrial action potential during fibrillation. This plot shows similar results to Fig. 2, where the length of the action potential increases and the peak of the action potential decreases with increasing concentrations of Quinidine. Only concentrations up to .4x were plotted, for concentrations greater than this amount completely nullified the action potential.

Figure 3
Fig. 4 shows the effects of varying CaL Channel conductance. As the conductance increases, the effective refractory period increases. An array of CaL conductances was used because there was insufficient data on the drug binding rate coefficients for Class IV antiarrhythmic drugs.  

Figure 4

Discussion

The developed model successfully answered the question of how quinidine binding affects atrial fibrillation action potential over time. The model was able to incorporate drug concentration and binding coefficients in an ordinary differential equation to model the effect of this drug on the heart. Importantly, it showed the clinically relevant effect of dampening and raising the duration of the action potential of an heart undergoing AF. 
The model demonstrated numerous strengths as a method of modeling drug binding and AF action potentials. The first of these is the ability to accurately reproduce an AF action potential. Alterations in transient outward current, inward rectifier current, and L-type calcium current successfully mimicked the effect of decreasing the action potential in AFIB. This also follows physiological sense, as chaotic signals reduce the overall conductance of these ion channels. The second strength of the model lies in the validation of the quinidine-binding action potential with literature examples (Figures 5-6).

                
Figures 5-6. Comparison of literature-defined effect of class IA drugs and results of mathematical model6
	As seen in the comparison of these graphs, the mathematical model effectively mimicked the increase in action potential duration and slight dampening effect of quinidine on a normal AP. This supports the viability of the model.
	Despite the literature-verified strengths, there were some inherent weaknesses in the model. The first of these is the fact that the model only incorporates binding and inactivation of sodium channels. Quinidine has exhibited effects on potassium and calcium channel gates as well, which the model fails to account for. If the binding coefficients and model equations for potassium, and calcium channels were experimentally determined, the model would describe total effect of drug binding significantly more accurately. The second weakness of the model is that it assumes total absorption of the dose of quinidine by the body. In realistic ingestion of the drug, the entire dosage would not be absorbed. This leads to overestimation in the effect of quinidine on the AF action potential and excessive dampening of the action potential. 
	This model matches well with similarly developed models in the literature. Shown below is a different model for AF action potential after treatment with quinidine. This model undergoes similar increase in AP duration and dampening of AP that our model successfully mimics (Figures 7-8).

	
Figures 7-8: Literature AFIB action potential after quinidine treatment5
A possible extension of this model exists in the form of incorporating other ionic currents into the model of drug binding. As mentioned earlier, quinidine has demonstrated the ability to block L-type calcium and potassium channels in addition to sodium channels. In order to create this model extension, binding coefficients and gate equations for these channels would be required. As an example, to simulate the effect of calcium channel binding, we varied the conductance of calcium channels and measured the effect on action potential(Figure in Results). This matched the literature-established plot for changing calcium channel conductance, validating the potential extension and opening the window for this model to be substantially improved. 

Figure 9. Effect of L-type calcium channel conductance on AP

	Using the results of these experiments and the incorporation of new ion channels, this model can serve as a method of testing the effect of current antiarrythmic drugs over time and developing new drugs to combat the effects of AFIB.

References

Courtemanche, M., Ramirez, R.J., & Nattel, S. (1998). Ionic Mechanisms Underlying Human Atrial Action Potential Properties: Insights from a Mathematical Model. American Journal of Physiology, 275(4), 301-321.

Courtemanche, M., Ramirez, R.J., & Nattel, S. (1999). Ionic Targets for Drug Therapy and Atrial-Fibrillation Induced Electrical Remodeling: Insights from a Mathematical Model. Cardiovascular Research, 42, 477-489. 

Healio. Atrial Fibrillation ECG Review. Retrieved from: https://www.healio.com/cardiology/learn-the-heart/ecg-review/ecg-topic-reviews-and-criteria/atrial-fibrillation-review.Medscape

Hondeghem, L. M., & Katzung, B. G. (1977). Time- and voltage-dependent interactions of antiarrhythmic drugs with cardiac sodium channels. Biochimica et Biophysica Acta (BBA), 72(3-4), 373-398. https://doi.org/10.1016/0304-4157(77)90003-X

Jain, D. K., Arya, R. K., & Jain, A. K. (2010). A mathematical model to understand the mechanisms of action of class 1 antiarrhymthmic drugs. Indian Journal of Pharmacology, 42(3), 195-196.

Kowey, P.R. (1998). Pharmacological Effects of Antiarrhythmic Drugs: Review and Update. JAMA Internal Medicine, 158(4), 325-332.

Luo, C. H. & Rudy, Y. (1994). A Dynamic Model of the Cardiac Ventricular Action Potential. I. Simulations of Ionic Currents and Concentration Changes. Circulation Research, 74, 1071-1096.

Medscape. Atrial Fibrillation Medication. Retrieved from: https://emedicine.medscape.com/article/151066-medication.

Odutayo, Ayodele; Wong, Christopher X; Hsiao, Allan J; Hopewell, Sally; Altman, Douglas G; Emdin, Connor A (6 September 2016). "Atrial fibrillation and risks of cardiovascular disease, renal disease, and death: systematic review and meta-analysis". BMJ: i4482. doi:10.1136/bmj.i4482.

Quinidine: MedlinePlus Drug Information. (2018, October 22). Retrieved from https://medlineplus.gov/druginfo/meds/a682396.html

Rivard, L; Khairy, P (December 2017). "Mechanisms, Clinical Significance, and Prevention of Cognitive Impairment in Patients With Atrial Fibrillation". Canadian Journal of Cardiology (Review). 33 (12): 1556–64. doi:10.1016/j.cjca.2017.09.024.PMID 29173598.

Xu, Lizhi et al. 3D multifunctional integumentary membranes for spatiotemporal cardiac measurements and stimulation across the entire epicardium. Nature Communications volume 5, Article number: 3329 (2014) doi:10.1038/ncomms4329




Appendix A: Glossary
Physiological Parameters: 
Cell modeled as a cylinder with L=100μm, r=11 μm
Cell Volume, V=38e-6 μL
[K+]o=5.4 mM, [K+]i=145 mM, 
[Na+]o=140 mM, [Na+]o=18 mM, 
Na-K exchange rate=.01833
[Ca2+]resting=2e-4 mM
ENa=54.4
Vinitial=-84 mV
Rate Constants:
kR=0
kA=26000
kI=0
lR=.05
lA=2.7
lI=6e-5



Appendix B: Matlab Code
%bme260model.m
%constants
del_t=0.04; 
v=-84.0;
 
g_Na=23;    
g_Ca=.3.*.09;    
g_K=0.282;   
g_K1=0.6047;  
g_to=.5*(.1652);
g_Kur=.5*(0.005+((0.05)./(1+exp(-1.*((v-15)./(13))))));
 
Ko = 5.4;
Ki = 145;     
Nao = 140;
Nai = 18;
Ca =2e-4;  
PRNaK = 0.01833;  
 
KQ10=3;
R = 8.315; 
T = 310;     
F = 96.49; 
 
E_Na=54.4;     
E_K=R*T/F*log((Ko + PRNaK*Nao)/(Ki + PRNaK*Nai));      
E_K1=R*T/F*log(Ko/Ki);  
E_Kp=E_K1;   
 
%gate equations 
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
 
%steady state gate equations
m=alpha_m./(alpha_m+beta_m);       
h=alpha_h./(alpha_h+beta_h); 
j=alpha_j./(alpha_j+beta_j); 
d=alpha_d./(alpha_d+beta_d); 
f=alpha_f./(alpha_f+beta_f); 
x=alpha_x./(alpha_x+beta_x); 
oi=(1+exp((v+43.1)/(5.3))).^(-1);
oa=(1+exp(-1.*((v+20.47)/(17.54)))).^(-1);
ua= (1+exp((v-99.45)/27.48)).^(-1);
 
t=0;
for i=1:15000   %15000 iterations = 150 ms
 
 I_stim =0; 
 
if i>500 & i<(525)
      I_stim=-40.0;
    else
      I_stim=0;
end
 
  
 E_Ca = 7.7 - 13.0287*log(Ca);    
 x_i=0.5*(1+sign(v+100-eps))*2.837.*(exp(0.04*(v+77))-1)./((v+77).*exp(0.04*(v+35)))+0.5*(1-sign(v+100-eps)); 

%Solve for drug binding effect over time
dX = zeros(2,1);
D = (2.55.*10.^(-5)).*(s-1)./5;
kr = 0;
ka = 26000;
ki = 0;
lr = 0.05;
la = 2.7;
li = 6.*10.^(-5);

alpha_hd=0.5*(1-sign((v+40)+40+eps))*0.135.*exp(-1*(80+(v+40))/6.8); 
beta_hd=(3.56*exp(0.079*(v+40))+3.1*1e5*exp(0.35*(v+40)))*0.5.*(1-sign((v+40)+40+eps))+(1./(0.13*(1+exp(-1*((v+40)+10.66)/11.1))))*0.5.*(1+sign((v+40)+40+eps));

if i == 1
    X = [0 0];
else
    [T,X] = ode15s(@(t,X) [alpha_hd.*(X(2) - X(1)) + (kr.*(h-((m.^3)*h)) + ka.*((m.^3)*h)).*D - beta_hd.*((X(1) - ((m.^3).*X(1))) + ((m.^3).*X(1))) - (lr.*(X(1) - ((m.^3).*X(1))) + la.*((m.^3).*X(1))); (kr.*(h-((m.^3)*h)) + ka.*((m.^3)*h) + ki.*(1 - X(2) - h)).*D - (lr.*(X(1) - ((m.^3).*X(1))) + la.*((m.^3).*X(1)) + li.*(X(2) - X(1)))],[0 t],[0 0]);
end

[MaxNa slsls] = max(X(:,2));

 %ionic currents 
 I_Na = g_Na*(m.^3)*h*j*(v-E_Na).*(1-MaxNa);     
 I_Ca = g_Ca*d*f*(v-E_Ca); 
 I_K = g_K*x*x_i*(v-E_K);    
 I_K1 = g_K1*(alpha_k1/(alpha_k1+beta_k1))*(v-E_K1);     
 I_Kp = 0.0183*(v-E_Kp)/(1+exp((7.488-v)/5.98));     
 I_b = 0.03921*(v+59.87);    
 I_to=g_to.*(oa.^3).*oi.*(v-E_K);
 I_Kur=g_Kur.*(ua.^3).*ua.*(v-E_K);
 
 %ion summation
 I_ion = I_Na + I_Ca + I_K + I_K1 + I_Kp + I_b + I_to + I_Kur;   
 
  
    %time constants 
    tau_m=1./(alpha_m+beta_m);
    m_inf=alpha_m./(alpha_m+beta_m); 
    tau_h=1./(alpha_h+beta_h);
    h_inf=alpha_h./(alpha_h+beta_h); 
    tau_j=1./(alpha_j+beta_j);
    j_inf=alpha_j./(alpha_j+beta_j); 
    tau_d=1./(alpha_d+beta_d);
    d_inf=alpha_d./(alpha_d+beta_d); 
    tau_f=1./(alpha_f+beta_f);
    f_inf=alpha_f./(alpha_f+beta_f); 
    tau_x=1./(alpha_x+beta_x);
    x_inf=alpha_x./(alpha_x+beta_x);
    tau_oa=((alpha_oa+beta_oa).^(-1))./KQ10; 
    oa_inf=(1+exp(-1.*((v+20.47)/(17.54)))).^(-1);
    tau_oi=((alpha_oi+beta_oi).^(-1))./KQ10;
    oi_inf=(1+exp((v+43.1)/(5.3))).^(-1);
    tau_ua=((alpha_ua+beta_ua).^(-1))./KQ10;
    ua_inf=(1+exp((v-99.45)/27.48)).^(-1);
    
    %diff eq. relating ion flow to membrane potential
   
    del_v = -I_ion-I_stim; 
 
 
 %gate fraction calculation
 m = m_inf - (m_inf-m)*exp(-del_t/tau_m);     
 h = h_inf - (h_inf-h)*exp(-del_t/tau_h);
 j = j_inf - (j_inf-j)*exp(-del_t/tau_j); 
 oi= oi_inf - (oi_inf-oi)*exp(-del_t/tau_oi); 
 oa= oa_inf - (oa_inf-oa)*exp(-del_t/tau_oa); 
 ua= ua_inf - (ua_inf-ua)*exp(-del_t/tau_ua);
 del_d = (d_inf-d)./tau_d;    
 del_f = (f_inf-f)./tau_f;    
 del_x = (x_inf-x)./tau_x;    
 del_Ca = -1e-4*I_Ca+0.07*(1e-4-Ca); 
  
  %forward euler
  v=v+del_t*del_v;
  d=d+del_t*del_d;
  f=f+del_t*del_f;
  x=x+del_t*del_x;
  Ca=Ca+del_t*del_Ca;
 
  
  %new voltage gating variable values
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
 
    
 %records values
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
clf;
hold on;
t=linspace(0,150,15000);
 
%plot(t,rec_v1,'r')
%plot(t,rec_v2,'r')
plot(t,rec_v,'k-')
%plot(t,rec_v3,'c')
%plot(t,rec_v4,'b')
%title('Membrane Potential (mV) vs. Time (ms) for Varying Calcium Channel Conductance (gCa)')
title('Membrane Potential (mV) vs. Time (ms)')
xlabel('Time (ms)')
ylabel('Membrane Voltage (mV)')
%legend('Atrial Action Potential During Fibrillation','Normal Atrial Action Potential')
%legend('50% Reduction in gTO','Normal Atrial Action Potential')
%legend('gCa=.03','gCa=.06','gCa=.09 (Normal)','gCa=.1.2','gCa=.1.5')



