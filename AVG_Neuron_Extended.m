%% Averaged-Neuron Extended Model %%
%
% Script for simulating membrane potential dynamics with changes in ion
% concentrations (K+, Ca2+ and Mg2+).
%
% Written by Mathias Heltberg, University of Copenhagen: mathias.heltberg@nbi.ku.dk. 
%
% Used the function 'goodplot.m' for generating the final plot.

clear all; close all; clc;

%% Perform Vm simulation %%

Le = 1000000; dt = 0.02;
time = linspace(dt,Le*dt,Le);

%%%% Currents
I = zeros(9,1);
Isav = zeros(9,Le);
Iextsav = zeros(3,Le);
Iext = zeros(3,1);

%%%% Variables
m2 = 0.01; h2 = 0.045; n3 = 0.54; m4 = 0.1; h4 = 0.045; 
m5 = 0.34; m6 = 0.01; m7 = 0.01; m8 = 0.01; h9 = 0.01;
s1 = 0.01; x2 = 0.01; s2 = 0.01; s3 = 0.01; Ca = 1.0; Kpl = 1.0;

%%%% Conductances of ion channels
g = zeros(9,1);
g1 = 0.03573;
g(1) = 0.03573; g(2) = 2.46; g(3) = 2.61868; g(4) = 1.79259; g(5) = 0.0350135; 
g(6) = 0.0256867; g(7) = 2.3; g(8) = 0.0717984; g(9) = 0.0166454;

gext = zeros(3,1);
gext(1) = 0.313425; gext(2) = 0.00434132; gext(3) = 0.00252916;

%%%% Potentials (in mV)
V = zeros(9,1); Vr = -55.0;
Vsav = zeros(Le,1); Casav = zeros(Le,1); Kplsav = zeros(Le,1);
V(1) = -60.95; V(2) = 55.0; V(3) = -100.0; V(6) = 120.0;
V(4) = V(3); V(5) = V(3); V(7) = V(3); V(8) = V(2); V(9) = V(3);
Vext = zeros(3,1); Vext(1) = 0; Vext(2) = 0; Vext(3) = -70.0;

%%%% Other Constants Needed
C = 1.0;
A = 0.02;
tau4 = 15.0; K = 30.0;
Ts1 = 2.0; Tx2 = 2.0; Ts2 = 100.0; Ts3 = 10.0; TCa = 121.4;
aCa = 0.5;

%%%% New Parameters
aKpl = 14.5;
TKpl = 1710.4;
Kons = 26.7137;

iterations = 1; % Number of iterations to perform

for runs = 1:iterations;
for i = 1:Le;    
    UpK = 1.0; UpNa = 0.08; UpCl = 0.1; 
    g(7) = 2.3*1.0; %Needed to put this here to get the right results when iterating multiple times

    %%%%% Here we change the ion concentrations
    
    if (i < Le/3) % Sleep state
        ccKo = 3.9*1000; ccKi = 140*1000; %%% Potassium
        ccNao = 140*1000; ccNai = 7.0*1000; %%% Sodium
        ccClo = 140*1000; ccCli = 7.0*1000; %%% Chloride
        ccCao = 1.35*1000; %%% Calcium
        ccMgo = 0.8; %%% Magnesium
        g(7) = 2.3*1.0;%%% Calcium-activated potassium channel

    elseif (i < 2*Le/3) % Awake state
        ccKo = 4.4*1000; ccKi = 140*1000; %%% Potassium
        ccNao = 140*1000; ccNai = 7.0*1000; %%% Sodium
        ccClo = 140*1000; ccCli = 7.0*1000; %%% Chloride
        ccCao = 1.2*1000; %%% Calcium
        ccMgo = 0.7; %%% Magnesium
        g(7) = 2.3*0.75; %%% Calcium-activated potassium channel

    else % Hyper awake state
        ccKo = 4.9*1000; ccKi = 140*1000; %%% Poassium
        ccNao = 140*1000; ccNai = 7.0*1000; %%% Sodium
        ccClo = 140*1000; ccCli = 7.0*1000; %%% Chloride
        ccCao = 1.05*1000; %%% Calcium 
        ccMgo = 0.7; %%% Magnesium
        g(7) = 2.3*0.75; %%% Calcium-activated potassium channel
    end
    
    ccCai = Ca;
    
    %%%%%%% Current 1
    V(1) = Kons*log( (UpK*ccKo + UpNa*ccNao + UpCl*ccCli) / (UpK*ccKi + UpNa*ccNai + UpCl*ccClo));
    I(1) = g(1)*(Vr - V(1));
    
    a = 0.1*(Vr + 33.)/(1 - exp(-(Vr + 33.)/10.0));
    b = 4.0*exp(-(Vr + 53.7)/12.0);
    m2 = a/(a + b);    
    a = 0.07*exp(-(Vr + 50.)/10.0);
    b = 1.0/(1.0 + exp(-(Vr + 20.0)/10.0));
    h2 = h2 + dt*( 4.*a*(1 - h2) - b*h2 ); 
    
    %%%%%%% Current 2
    V(2) = Kons*log( (ccNao) / (ccNai)); %%%% 
    I(2) = g(2)*(m2^3)*h2*(Vr-V(2));

    %%%%%%% Current 3
    a = 0.01 * (Vr + 34.0)/(1 - exp(-(Vr + 34.0)/10.0));
    b = 0.125*exp(-(Vr + 44.0)/25.0);
    n3 = n3 + dt*( 4.0*(a*(1.0 - n3) - b*n3));
    V(3) = Kons*log((ccKo) / (ccKi));     
    I(3) = g(3).*(n3.^4).*(Vr-V(3));
    
    %%%%%%% Current 4   
    m4 = 1.0/(1.0 + exp(-(Vr + 50.0 )/20.0 ));
    a = 1.0/(1.0 + exp( (Vr + 80.0)/ 6.0 ));
    h4 = h4 + dt*((a - h4)/tau4 );
    V(4) = Kons*log( (ccKo) / (ccKi));
    I(4) = g(4)*(m4^3)*h4*(Vr - V(4));
    
    %%%%%%% Current 5   
    a = 1.0/(1.0 + exp( -( Vr + 34.0 )/6.5));
    b = 8.0/(exp(-(Vr + 55.0)/30.0 ) + exp(( Vr + 55.0)/30.0));
    m5 = m5 + dt*( (a - m5)/b);
    V(5) = Kons*log((ccKo) / (ccKi));
    I(5) = g(5)*m5*(Vr - V(5));
 
    %%%%%%% Current 6     
    m6 = 1.0/(1.0 + exp(-( Vr + 20.0)/9.0));
    V(6) = Kons*log((ccCao) / (ccCai));
    I(6) = g(6)*(m6^2)*(Vr - V(6));

    %%%%%%% Current 7
    Kpl = Kpl + dt*(-aKpl*(10.0*A*(-I(3) - I(4) - I(5) - I(7) + I(9))) - Kpl/TKpl);
    Ca = Ca + dt*(-aCa*(10.0*A*I(6) + Iext(2)) - Ca/TCa);
    m7 = 1.0/( 1.0 + (K/Ca)^(3.5));
    V(7) = Kons*log((ccKo) / (ccKi));
    I(7) = g(7)*m7*(Vr - V(7));
 
    %%%%%%% Current 8     
    m8 = 1.0/(1.0 + exp(-(Vr + 55.7 )/7.7 ));
    V(8) = Kons*log((ccNao) / (ccNai));
    I(8) = g(8)*m8^3*(Vr-V(8));
     
    %%%%%%% Current 9    
    h9 = 1.0/( 1.0 + exp((Vr + 75.0)/4.0) );
    V(9) = Kons*log((ccKo) / (ccKi));
    I(9) = g(9)*h9*(Vr - V(9));
         
    %%%%%%%%%%%%%%  Extrinsic currents  
   
    fV = 1.0/(1.0 + exp(-(Vr - 20.0)/2.0));
    
    s1 = s1 + dt*(3.48*fV - s1/Ts1);
    UpK = 1.0; UpNa = 1.0; 
    Vext(1) = Kons*log((UpK*ccKo + UpNa*ccNao) / (UpK*ccKi + UpNa*ccNai));
    Iext(1) = gext(1)*s1*(Vr - Vext(1));
            
    x2 = x2 + dt*( 3.48*fV - x2/Tx2);
    s2 = s2 + dt*( 0.5*x2*(1.0 - s2) - s2/Ts2);
    UpK = 1.0; UpNa = 1.0; UpCa = 1.0;
    Vext(2) = Kons*log( (UpK*ccKo + UpNa*ccNao + UpCa*ccCao) / (UpK*ccKi + UpNa*ccNai + UpCa*ccCai));
    Iext(2) = 1.1/(1.0 + ccMgo/(0.8*1000))*gext(2)*s2*(Vr - Vext(2));        
    
    s3 = s3 + dt*( fV - s3/Ts3);
    Iext(3) = gext(3)*s3*(Vr - Vext(3));

    Vk1 = dt*(-10.0*A*sum(I)/(10.0*A*C) - sum(Iext)/(10.0*A*C));
    Vr = Vr + dt*(-10.0*A*sum(I)/(10.0*A*C) - sum(Iext)/(10.0*A*C) + normrnd(0.0,3.0));

    Vsav(i,runs) = Vr;    
    Casav(i,runs) = Ca;
end
end

%% Plotting simulated Vm dynamics

V = Vsav; % Raw membrane potential 
V_filtered = movmedian(V, 4000); % Spikes removed
Ca_i = Casav; % Intracellular Ca2+ 

figure(1);
plot(time/1000.,V, 'k', 'LineWidth', 0.3); hold on;
plot(time/1000.,V_filtered, 'r', 'LineWidth',2); hold on;
plot(time/1000.,Ca_i, 'b', 'LineWidth',2); hold on;
xlabel('Time [s]'); ylabel('Potential [mV]');
legend('Raw Vm', 'Filtered Vm', '[Ca2+]i','Orientation','horizontal');
goodplot;


