clear
clc
close all
tic 
Dp0 = 0.000200000000000000;
Qmilk = 1750/3600;
Gdry = 72000/3600;
Tg0 = 403;
tic
Dp0vct = linspace(log10(0.08*10^-3),log10(1*10^-3),40);
Dp0vct = 10.^Dp0vct;
t2zerovct = zeros(length(Dp0vct),1);
ZWendvct = [t2zerovct,t2zerovct,t2zerovct];

t2zerovct2 = zeros(length(Dp0vct),1);
Tg0vct = linspace(360,420,length(Dp0vct));
ZWendvct2 = [t2zerovct,t2zerovct,t2zerovct];
Tpendarray = zeros(length(Dp0vct),length(Dp0vct));
Wendarray = Tpendarray;




for i=1:length(Dp0vct)
    t2zerovct(i) = fzero(@(tend)f2zero(Dp0vct(i),tend,Qmilk,Gdry,Tg0),[0.0001,120]);
    ZWendvct(i,:) = Z_W_Tpend_finder(Dp0vct(i),t2zerovct(i),Qmilk,Gdry,Tg0);
end

Zend = 30;

for i=1:length(Dp0vct)
    t2zerovct2(i) = fzero(@(tend)f2zero2(Dp0vct(i),tend,Qmilk,Gdry,Tg0,Zend),[0.001,120]);
    ZWendvct2(i,:) = Z_W_Tpend_finder(Dp0vct(i),t2zerovct2(i),Qmilk,Gdry,Tg0);
end
for j=1:length(Tg0vct)
    for i=1:length(Dp0vct)
        t2zerovct3 = fzero(@(tend)f2zero2(Dp0vct(i),tend,Qmilk,Gdry,Tg0vct(j),Zend),[0.0001,120]);
        Z_W_Tpend = Z_W_Tpend_finder(Dp0vct(i),t2zerovct3,Qmilk,Gdry,Tg0vct(j));
        Wendarray(i,j) = Z_W_Tpend(2);
        Tpendarray(i,j) = Z_W_Tpend(3);
    end
end


figure(1)
plot(Dp0vct(:),ZWendvct(:,1))
title('Relationship between starting diameter and spray dryer length')
xlabel('Droplet intial diameter [m]')
ylabel('required length for the spray dryer [m]')

figure(2)
plot(Dp0vct(:),ZWendvct2(:,2))
title('Final humidity of the particle with a 30 m spray dryer')
xlabel('Droplet intial diameter [m]')
ylabel('particle humidity [kgH2O/kgparticle]')


figure(3)
plot(Dp0vct(:),ZWendvct2(:,3))
title('Final temperature of the particle with a 30 m spray dryer')
xlabel('Droplet intial diameter [m]')
ylabel('Particle final temperature [K]')


figure(4)
surf(Tg0vct,Dp0vct,Wendarray)
title('Final humidity of the particle with a 30 m spray dryer')
ylabel('Droplet intial diameter [m]')
xlabel('Air inlet temperature [K]')
zlabel('particle humidity [kgH2O/kgparticle]')

figure(5)
surf(Tg0vct,Dp0vct,Tpendarray)
title('Final temperature of the particle with a 30 m spray dryer')
ylabel('Droplet intial diameter [m]')
zlabel('particle temperature [K]')
xlabel('Air inlet temperature [K]')
toc

function out = f2zero(Dp0,tend,Qmilk,Gdry,Tg0)
Qfat = Qmilk*4.76/100;

Dsd = 5.5;

vp0 = 0.3;

rhop = 1040;

P = 101325;

Tp0 = 303;

A = pi*Dsd.^2./4;
rhog0  =  0.867;%[kg/m3]
vg0  =  Gdry/rhog0./A;%[m/s]
vs0 = vp0-vg0;
Hup0 = (Qmilk-Qfat)/Qfat;%[kg/kg]
mp0 = rhop*pi*Dp0^3/6;%[kg]
mdry = mp0/(1+Hup0);%[kg]
eta_A_vp = Qmilk/mp0;%[drops/s];
eta_A_vp = floor(eta_A_vp);

y0 = [mp0,0,Tp0,Tg0,vs0,0];

[t,ys] = ode15s(@(t,Y)SPRAYDRYER(t,Y,P,Gdry,A,eta_A_vp,mdry),[0,tend],y0);

mp = ys(:,1);

W = mp/mdry - 1;

out = W(end) - 0.005;

end

function out = f2zero2(Dp0,tend,Qmilk,Gdry,Tg0,Zend)
Qfat = Qmilk*4.76/100;

Dsd = 5.5;

vp0 = 0.3;

rhop = 1040;

P = 101325;

Tp0 = 303;

A = pi*Dsd.^2./4;
rhog0  =  0.867;%[kg/m3]
vg0  =  Gdry/rhog0./A;%[m/s]
vs0 = vp0-vg0;
Hup0 = (Qmilk-Qfat)/Qfat;%[kg/kg]
mp0 = rhop*pi*Dp0^3/6;%[kg]
mdry = mp0/(1+Hup0);%[kg]
eta_A_vp = Qmilk/mp0;%[drops/s];
eta_A_vp = floor(eta_A_vp);

y0 = [mp0,0,Tp0,Tg0,vs0,0];

[t,ys] = ode15s(@(t,Y)SPRAYDRYER(t,Y,P,Gdry,A,eta_A_vp,mdry),[0,tend],y0);

Z = ys(:,6);

out = Z(end) - Zend;

end


function out = Z_W_Tpend_finder(Dp0,tend,Qmilk,Gdry,Tg0)

Qfat = Qmilk*4.76/100;

Dsd = 5.5;

vp0 = 0.3;

rhop = 1040;

P = 101325;

HupF = 0.005;

Tp0 = 303;

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

[t,ys] = ode15s(@(t,Y)SPRAYDRYER(t,Y,P,Gdry,A,eta_A_vp,mdry),[0,tend],y0);

mp = ys(:,1);
Z = ys(:,6);
Tp = ys(:,3);

W = mp/mdry - 1;

out(1) = Z(end);
out(2) = W(end);
out(3) = Tp(end);

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
    rhop = - 8.2961733900239 * W + 1205;
    Vp = mp/rhoL;
    Vpstore = Vp;
else
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

dmp = Kp*Sp*(Pw-P0)*heaviside(W-0.0005);
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