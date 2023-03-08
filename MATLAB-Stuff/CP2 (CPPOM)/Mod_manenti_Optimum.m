% This matlab's main purpose is to calculate the best operating conditions
% taking into account both opex and capex, it's basically a non linear
% optimisation, you got the output to be minimised at line 91, the weights
% are not "realistic", you need to modify them.
% please be aware that it'll take a while for ga to properly converge, and
% even then you'll not be sure if your results is truly the global minimum.
% nonlincon and fmincon work much faster than ga, but they reach usually an
% optimum configuration whose goal function is usually higher than ga.

% this model does not include the drop2particle hypothesis or anything else
% that you may consider as "advanced", it's the SuPER's model with slight
% modification to ensure that non physical behaviour do not occur (such as
% water exchange between the particle and air even when W is lower than 0)



clear
close all
Qmilk = 1750 / 3600;
Dsd = 5.5;
Gdry = 72000 / 3600;
P = 101325 * 1.3;
Tg0 = 378;
tend0 = 7;

x0 = [Qmilk,Dsd,Gdry,P,Tg0,tend0];


opts = optimoptions('ga','PlotFcn',@gaplotbestf);
lb = [1200/3600,0.5,15,50e+03,366.15,3];
ub = [3600/ 3600,6.5,30,500e+03,588.15,15];
[opt_values_ga,fvalga] = ga(@(Qmilk_Dsd_Gdry_P_Tg0_tend)fun2min(Qmilk_Dsd_Gdry_P_Tg0_tend),6,[],[],[],[],lb,ub,[],opts)

[opt_values_fmincon,fval] = fmincon(@(Qmilk_Dsd_Gdry_P_Tg0_tend)fun2min(Qmilk_Dsd_Gdry_P_Tg0_tend),opt_values_ga,[],[],[],[],lb,ub)

profits_ga = fun_eval_profits(opt_values_ga)

profits_fmincon = fun_eval_profits(opt_values_fmincon)

function out = fun2min(Qmilk_Dsd_Gdry_P_Tg0_tend)

Qmilk = Qmilk_Dsd_Gdry_P_Tg0_tend(1);
Dsd = Qmilk_Dsd_Gdry_P_Tg0_tend(2);
Gdry = Qmilk_Dsd_Gdry_P_Tg0_tend(3);
P = Qmilk_Dsd_Gdry_P_Tg0_tend(4);
Tg0 = Qmilk_Dsd_Gdry_P_Tg0_tend(5);
tend = Qmilk_Dsd_Gdry_P_Tg0_tend(6);
Qfat = Qmilk * 4.76 / 100;
vp0 = 0.3;
rhop = 1040;
Dp0 = 0.2 * 10 ^ (-3);

Tp0 = 303;

A = pi * Dsd ^ 2 / 4;
rhog0 = 0.867;
vg0 = Gdry / rhog0 / A;
vs0 = vp0 - vg0;
Hup0 = (Qmilk - Qfat) / Qfat;
mp0 = rhop * pi * Dp0 ^ 3 / 6;
mdry = mp0 / (1 + Hup0);
eta_A_vp = Qmilk / mp0;
eta_A_vp = floor(eta_A_vp);
y0 = [mp0,0,Tp0,Tg0,vs0,0];
[t,Y] = ode113(@(t,Y)SPRAYDRYER(t,Y,P,Gdry,A,eta_A_vp,mdry),[0,tend],y0);
L = Y(end,6);
Wend = Y(end,1) / mdry - 1;
Tmaxp = max(Y(:,3));
D_out = Dsd+0.05;
mass_SD = pi*(D_out^2-Dsd^2)/4*L*8000;
delp = P-101325;
w_vac = 0.0044*4;
w_hp = 0.0044;
if delp<0
    w_P = w_vac;
else
    w_P = w_hp;
end
Tbase_aft_compression = (513-303)*(delp)/(5e+05);
out = -mass_SD*10/(20*365*24*3600) - (Tg0 - Tbase_aft_compression)*Gdry *1.33e-07*1  + Qmilk* 4.76 / 100*4.5 - Qmilk*0.15 - abs(delp)/1e+05 * w_P;

out = -out+3;
if Wend>0.0005
    out = out+25;
elseif Tmaxp>345
    out = out+25;
elseif L>19
    out = out+25;
else
    out = out;
end

end


function out = fun_eval_profits(Qmilk_Dsd_Gdry_P_Tg0_tend)

Qmilk = Qmilk_Dsd_Gdry_P_Tg0_tend(1);
Dsd = Qmilk_Dsd_Gdry_P_Tg0_tend(2);
Gdry = Qmilk_Dsd_Gdry_P_Tg0_tend(3);
P = Qmilk_Dsd_Gdry_P_Tg0_tend(4);
Tg0 = Qmilk_Dsd_Gdry_P_Tg0_tend(5);
tend = Qmilk_Dsd_Gdry_P_Tg0_tend(6);
Qfat = Qmilk * 4.76 / 100;
vp0 = 0.3;
rhop = 1040;
Dp0 = 0.2 * 10 ^ (-3);

Tp0 = 303;

A = pi * Dsd ^ 2 / 4;
rhog0 = 0.867;
vg0 = Gdry / rhog0 / A;
vs0 = vp0 - vg0;
Hup0 = (Qmilk - Qfat) / Qfat;
mp0 = rhop * pi * Dp0 ^ 3 / 6;
mdry = mp0 / (1 + Hup0);
eta_A_vp = Qmilk / mp0;
eta_A_vp = floor(eta_A_vp);
y0 = [mp0,0,Tp0,Tg0,vs0,0];
[t,Y] = ode113(@(t,Y)SPRAYDRYER(t,Y,P,Gdry,A,eta_A_vp,mdry),[0,tend],y0);
L = Y(end,6);
Wend = Y(end,1) / mdry - 1;
Tmaxp = max(Y(:,3));
D_out = Dsd+0.05;
mass_SD = pi*(D_out^2-Dsd^2)/4*L*8000;
delp = P-101325;
w_vac = 0.0044*4;
w_hp = 0.0044;
if delp<0
    w_P = w_vac;
else
    w_P = w_hp;
end
Tbase_aft_compression = (513-303)*(delp)/(5e+05);
out = -mass_SD*10/(20*365*24*3600) - (Tg0 - Tbase_aft_compression)*Gdry *1.33e-07*1  +  Qmilk* 4.76 / 100*4 - Qmilk*0.15 - abs(delp)/1e+05 * w_P;

out = out*3600*24*365;
end


function dY = SPRAYDRYER(t,Y,P,Gdry,A,eta_A_vp,mdry)
global Vpstore
mp = Y(1);  Gi = Y(2); Tp = Y(3); Tg = Y(4); vs = Y(5); Z = Y(6);


MwH2O = 0.018;
g = 9.81;
rhoL = 1040;
R = 8.314;

W = mp/mdry - 1;

if W>=0.714
    Vp = mp/rhoL;
    Vpstore = Vp;
    rhop = - 8.2961733900239 * W + 1205;
elseif W<0.714
    Vp = Vpstore;
    rhop = mp/Vp;
end
if rhop<715
   rhop = 715;
else
end

% quantities changing at each iteration
Diff = 10^-4 * (0.242 + (0.399-0.242)*(Tg-293.15)/80);
cpg = 1006 + (1018-1006)*(Tg - 303.15)/(423.15 - 303.15);
mug = 10^-6 * ( 17.15 + (23.8-17.15)*(Tg - 273.15)/150);
DHev = 10^3*(2429.8 + (2113.7-2429.8)*(Tp - 303.15)/120);
cpmilk = (0.3828 * mdry/mp +  (mp-mdry)/mp) *4180;
kg = 10^-3 * (24.36 + (35 - 24.36)*(Tg - 273.15)/150);


Dp = (6*Vp/pi)^(1/3);
Sp = pi*Dp^2;
G = Gdry+Gi;
MMair = (Gdry/(G*0.02897) + Gi/G/0.018)^-1;
rhog = P * MMair/(R*Tg);%[kg/m3]
vg = G./rhog./A;%[m/s]
Re = rhog.*abs(vs).*Dp/mug;
f=24./Re.*(1+0.14*Re.^0.7);
Sc = mug/Diff/rhog;
Pr = mug*cpg/kg;
Sh = 2+0.4*Re.^0.5*Sc^(1/3);
Nu = 2+0.4*Re.^0.5*Pr^(1/3);

Kc = Sh.*Diff/Dp;%[m/s]
h = Nu.*kg/Dp;%[J/m2/s/K]
Kp = Kc.*MwH2O/R/((Tg+Tp)/2);%[kg/s/Pa/m2]


Pw = MMair/MwH2O*Gi*P/Gdry;%[Pa]
P0 = PoH20_T(Tp); %[Pa]


%ODE SYSTEM

dmp = heaviside(W-0.0005)*Kp*Sp*(Pw-P0);
dGi = -dmp*eta_A_vp;
dTp = 1/(mp*cpmilk)*(h*Sp*(Tg-Tp)+dmp*DHev);
dTg = 1/(G*cpg)*h*Sp*(Tp-Tg)*eta_A_vp;
dvs = 1/mp*((rhop-rhog)*Vp*g-1/2 * f * rhog * vs * abs(vs) * Dp^2 /4 *pi-vs*dmp);     
dZ = vs+vg;
dY = [dmp dGi dTp dTg dvs dZ]';

end

function out = PoH20_T(T)
% T in Kelvins
if T <=303
    out = 10.^(5.40221 - 1838.675./(T-31.737));
elseif T>303 && T<=333
    out = 10^(5.20389 - 1733.926/(T-39.485));
elseif T>333 && T<=363
    out = 10^(5.0768 - 1659.793/(T-45.854));
elseif T>364 && T<=373
    out = 10^(5.08354 - 1663.125/(T-45.622));
else 
    out = 10^(3.55959 - 643.748/(T-198.043));
end
out = out *10^5;
end