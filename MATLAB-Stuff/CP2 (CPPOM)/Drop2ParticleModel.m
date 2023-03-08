clear
clc
close all
tic
% best opex-capex from ga at [6.94376475e-01 6.42647653e+00 1.78114900e+01 5.78253102e+04
% 3.92502122e+02 1.49862419e+01] at x = [Qmilk,Dsd,Gdry,P,Tg0,tend0]
% best opex-capex from minimize at x_opt_minimise = [ 1.492e+01,  1.391e+01,  1.725e+02,  5.783e+04,  5.109e+02,
%            1.091e+01]
%Qmilk = 1750/3600;
Qmilk = 1.492e+01;
%Qmilk = 5.725e-01;

Qfat = Qmilk*4.76/100;

%Dsd = 5.5;
Dsd = 1.391e+01;
%Dsd = 6.378e+00;

%Gdry = 72000/3600*1;
Gdry = 1.725e+02;
%Gdry = 1.784e+01;

vp0 = 0.3;

rhop = 1040;

%P = 101325*1;
P =  5.783e+04;
%P = 5.783e+04;

Dp0 = 0.2*10^-3;

HupF = 0.005;

Tp0 = 303;

%Tg0 = 403;
Tg0 = 5.109e+02;
%Tg0 = 3.924e+02;

A = pi*Dsd.^2./4;
rhog0  =  0.867;%[kg/m3]
vg0  =  Gdry/rhog0./A;%[m/s]
vs0 = vp0-vg0;
Hup0 = (Qmilk-Qfat)/Qfat;%[kg/kg]
mp0 = rhop*pi*Dp0^3/6;%[kg]
mdry = mp0/(1+Hup0);%[kg]
mpF = mdry*(1+HupF);%[kg]
eta_A_vp = Qmilk/mp0;%[drops/s];
eta_A_vp = floor(eta_A_vp);

y0 = [mp0,0,Tp0,Tg0,vs0,0];
%tend = 6;
tend = 1.091e+01;
%tend = 1.499e+01;

[t,ys] = ode89(@(t,Y)SPRAYDRYER(t,Y,P,Gdry,A,eta_A_vp,mdry),[0,tend],y0);

mp = ys(:,1);
Gi = ys(:,2);
Tp = ys(:,3);
Tg = ys(:,4);
vs = ys(:,5);
Z = ys(:,6);

W = mp/mdry - 1;
toc
figure(1)
subplot(2,2,1)
plot(Z,W)
title('Humidity over spray dryer length')
xlabel('length [m]')
ylabel('humidity [kgh20/kgdry]')
ylim([-2,22])

subplot(2,2,2)
plot(Z,Tg)
title('Gas temperature over spray dryer length')
xlabel('length [m]')
ylabel('gas temperature [K]')
ylim([340,420])

subplot(2,2,3)
plot(Z,Tp)
title('Particle temperature over spray dryer length')
xlabel('length [m]')
ylabel('particle temperature [K]')
ylim([280,360])

subplot(2,2,4)
plot(Z,vs)
title('particle speed over spray dryer length')
xlabel('length [m]')
ylabel('entrainment speed [m/s]')
ylim([-1,1])




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