%[CH4 H20 CO H2 CO2 N2]
clear
clc
close all
%%DATI CONDUCIBILITA'
metano_cond=[8.3983e-6 1.4268 -49.654 0];
H2O_cond=[6.2041e-6 1.3973 0 0];
CO_cond=[0.00059882 0.6863 57.13 501.92];
H2_cond=[0.002653 0.7452 12 0];
CO2_cond=[3.69 -0.3838 964 1860000];
N2_cond=[0.00033143 0.7722 16.323 373.72];
specie_cond=[ metano_cond; H2O_cond; CO_cond; H2_cond; CO2_cond; N2_cond];


%DATI VISCOSITA'
metano_viscosity=[5.2546e-7 0.59006 105.67 0];
H2O_viscosity=[1.7096e-8 1.1146 0 0];
CO_viscosity=[1.1127e-6 0.5338 94.7 0];
H2_viscosity=[1.797e-7 0.685 -0.59 140];
CO2_viscosity=[2.148e-6 0.46 290 0];
N2_viscosity=[6.5592e-7 0.6081 54.714 0];
specie_viscosity=[ metano_viscosity; H2O_viscosity; CO_viscosity; H2_viscosity; CO2_viscosity; N2_viscosity];


%%DATI QUALSIASI
N_tubi=80;
Q=18e6;%W
H_furnace=10;%m
n_in=1050;%kmol/h
x_in=[0.222 0.77 0 0 0.005 0.003];
PM=[16.043 18 28.01 2 44.01 14.0067];%kg/kmol
PM_mix=sum(PM.*x_in);
w_in=(x_in.*PM)/PM_mix;
nu=[-1 -1 1 3 0 0;0 -1 -1 1 1 0;-1 -2 0 4 1 0];
m_in=(n_in*PM_mix)/3600;

%%DATI CP
methane = [0.3330e5 0.7993e5 2.0869e3 0.4160e5 991.96];
H2 = [0.2896e5 0.0939e5 3.0120e3 0.0758e5 1484 ];
water=[0.33363e5 0.2679e5 2.6105e3 0.08896e5 1169];
CO=[0.29108e5 0.08773e5 3.0851e3 0.084553e5 1538.2];
CO2=[0.2937e5 0.3454e5 1.428e3 0.264e5 588];
N2=[0.29105e5 0.086149e5 1.7016e3 0.0010347e5 909.79];
matrix_cp=[methane;water;CO;H2;CO2;N2];
P_in = 21.2; % bar
type_tubo = 1;
type_cat =2;
type_furnace=1;
t = 3; %anni
y0=[w_in, 800,P_in];
choices=[type_tubo,type_cat,type_furnace];

[z,ys] = ode113(@(z,y)PFR_rea(z,y,type_tubo,type_cat,type_furnace,m_in,t),[0,10],y0);

mus = ys(:,1:6);
PM_mix=1./sum(mus./PM,2);
xs=(mus.*PM_mix)./PM;
x_out = xs(end,:);
n_out = m_in/PM_mix(end)*3600; %kmol/h
Ts = ys(:,7);
Ps = ys(:,8);
figure(1)
for i=1:6
    plot(z,xs(:,i));hold on
end
hold off
figure(2)
plot(z,Ts)
figure(3)
plot(z,Ps)
conv_metano = (n_in*x_in(1)-n_out*x_out(1))/(n_in*x_in(1))
prod_H2 = n_out*x_out(4)
disp('[CH4   H20   CO   H2   CO2   N2]  ')
disp(xs(end,:))

if type_tubo==1
    d_tubi_in=3;
elseif type_tubo==2
    d_tubi_in=4;
elseif type_tubo==3
    d_tubi_in=5;
else
end
%%CATALIZZATORE
if type_cat==1
    rho_bed=1058;%kg/m^3
    cost_cat=50;%$/kg
elseif type_cat==2
    rho_bed=556.3;%kg/m^3
    cost_cat=70;%$/kg
else
end
cost_ducts = (15000 + 1000*(d_tubi_in - 3))*N_tubi;
cost_cat = (pi*(d_tubi_in*0.0254)^2/4*rho_bed*10)*cost_cat*N_tubi;
m_cat = (pi*(d_tubi_in*0.0254)^2/4*rho_bed*10)*N_tubi
cost_tot = cost_cat+cost_ducts %$



function dydt=PFR_rea(z,Y,type_tubo,type_cat,type_furnace,m_in,t)
% input
PM=[16.043 18 28.01 2 44.01 14.0067];
wi(1:6)=Y(1:6);
T=Y(7);
P=Y(8);
PM_mix=1./sum(wi./PM);
x=(wi.*PM_mix)./PM;
% dati
% termodinamica
R = 8.314;
delHo = 1e+07.*[-7.5420, -24.1814, -11.0530, 0, -39.3510,  0];
delG0 = 1e+07*[-5.0496,-22.8590, -13.7150,0,-39.4370,0];
nu=[-1 -1 1 3 0 0;0 -1 -1 1 1 0;-1 -2 0 4 1 0];
% conducibilità
metano_cond=[8.3983e-6 1.4268 -49.654 0];
H2O_cond=[6.2041e-6 1.3973 0 0];
CO_cond=[0.00059882 0.6863 57.13 501.92];
H2_cond=[0.002653 0.7452 12 0];
CO2_cond=[3.69 -0.3838 964 1860000];
N2_cond=[0.00033143 0.7722 16.323 373.72];
specie_cond=[ metano_cond; H2O_cond; CO_cond; H2_cond; CO2_cond; N2_cond];
kt_tubo = 29.8; %W/(m*k)
% viscosità
metano_viscosity=[5.2546e-7 0.59006 105.67 0];
H2O_viscosity=[1.7096e-8 1.1146 0 0];
CO_viscosity=[1.1127e-6 0.5338 94.7 0];
H2_viscosity=[1.797e-7 0.685 -0.59 140];
CO2_viscosity=[2.148e-6 0.46 290 0];
N2_viscosity=[6.5592e-7 0.6081 54.714 0];
specie_viscosity=[ metano_viscosity; H2O_viscosity; CO_viscosity; H2_viscosity; CO2_viscosity; N2_viscosity];
% calori specifici
methane = [0.3330e5 0.7993e5 2.0869e3 0.4160e5 991.96];
H2 = [0.2896e5 0.0939e5 3.0120e3 0.0758e5 1484 ];
water=[0.33363e5 0.2679e5 2.6105e3 0.08896e5 1169];
CO=[0.29108e5 0.08773e5 3.0851e3 0.084553e5 1538.2];
CO2=[0.2937e5 0.3454e5 1.428e3 0.264e5 588];
N2=[0.29105e5 0.086149e5 1.7016e3 0.0010347e5 909.79];
matrix_cp=[methane;water;CO;H2;CO2;N2];
% calcolo proprietà 
%[CH4 1  H20 2 CO 3 H2 4 CO2 5 N2 6 ]
cp_mix=cpmix(matrix_cp,x,T)/PM_mix;
kt_miscela_1=cond_mix(specie_cond,T,x,PM);%W/mK
Xnew = [x(6), 0, x(5), x(2), 0, x(3), x(4), x(1), 0, 0, 0, 0, 0];
kt_miscela_2 = lambda(Xnew, T);
kt_miscela = (kt_miscela_1*0.85+kt_miscela_2*0.15);
visc_mix_2 = mu(Xnew,T);
visc_mix_1=viscosity_mix(specie_viscosity,T,x,PM);
visc_mix = (visc_mix_1+visc_mix_2)/2;
%%TUBI
if type_tubo==1
    d_tubi_in=0.0762;%m
    s=0.012;
    d_tubi_ext=d_tubi_in+2*s;
elseif type_tubo==2
    d_tubi_in=0.1016;%m
    s=0.015;
    d_tubi_ext=d_tubi_in+2*s;
elseif type_tubo==3
    d_tubi_in=0.127;%m
    s=0.018;
    d_tubi_ext=d_tubi_in+2*s;
else
end
Sv=4/d_tubi_in;
G=m_in/(pi*d_tubi_in^2/4*80);
rhomix = rhomixgp(PM,x,P*10^5,T);
v = G/rhomix;
%%CATALIZZATORE
if type_cat==1
    a0=70;%W/m2K
    dp=0.59e-2;%m
    eps=0.55;
    rho_bed=1058;%kg/m^3
    Re=G*dp/visc_mix;
    Pr=visc_mix*cp_mix/kt_miscela;
    f=(10.5*(1-eps)^1.2/eps^3)*Re^-0.3;
    Nu=(a0*dp/kt_miscela)+0.25*(Re^0.72)*(Pr^0.33);
    hcoeff = Nu*kt_miscela/dp;
    eff=0.01;
    cost_cat=50;%$/kg
elseif type_cat==2
    a0=70;%W/m2K
    dp=0.86e-2;%m
    eps=0.62;
    rho_bed=556.3;%kg/m^3
    Re=G*dp/visc_mix;
    Pr=visc_mix*cp_mix/kt_miscela;
    f=(4.63*(1-eps)^1.2/eps^3)*(Re^-0.16);
    Nu=(a0*dp/kt_miscela)+0.15*(Re^0.76)*(Pr^0.33);
    hcoeff = Nu*kt_miscela/dp;
    eff=0.02;
    cost_cat=70;%$/kg
else
end
%%FORNACE
if type_furnace==1
    f1=17.22;
    f2=0.0578;
    f3=3.930;
    f4=0.816;
elseif type_furnace == 2
    f1=5.819;
    f2=0.170;
    f3=2.791;
    f4=1.138;
else
end
L=10;%m
Qtot = 18*10^6; % W
N = 80;
Q_tubo=Qtot/(N*L*d_tubi_in*pi); %W/m^2
Q=Q_tubo*(f1*(z/L+f2)*exp(-f3*(z/L)^f4));%W/m^2

%%CINETICA
%[CO H2 CH4 H20]
Kj_0=[8.23e-5 6.12e-9 6.65e-4 1.77e5];
dHj=[-70.65 -82.9 -38.28 88.68]*1e+03;%J/mol
Kj=Kj_0.*exp(-dHj./(R*T));
K_CO=Kj(1);
K_H2=Kj(2);
K_Met=Kj(3);
K_H20=Kj(4);
Ea=[265 88 275]*1e+03;%J/mol
kj0=[13.2 136 1.68];
kj=kj0.*exp(-Ea./R*(1/T-1/873));
th = thermo_funs;
nrea = size(nu,1);
Keqs = zeros(1,nrea);
nspecies = size(nu,2);
g = zeros(1,nspecies);
for i=1:nspecies
    cps = [matrix_cp(i,:) 0 0 0 0 0 0 0 0 0 0 6000 6200];
    g(i) = th.delg(cps,delG0(i),delHo(i),T,'PERRY');
end
delGr = zeros(1,nrea);

for i=1:nrea
    delGr(i) = sum(nu(i,:).*g);
    Keqs(i) = exp(-delGr(i)./(R*1000*T));
end
Keq_1 = Keqs(1);
Keq_2 = Keqs(2);
Keq_3 = Keqs(3);
%[CH4 1  H20 2 CO 3 H2 4 CO2 5 N2 6 ]
if x(4)>1e-06
    rden = (1+K_CO*P*x(3)+K_H2*P*x(4)+K_Met*x(1)*P+K_H20*P*x(2)/(P*x(4)))^2;
elseif x(4)<1e-06
    rden = (1+K_CO*P*x(3)+K_H2*P*x(4)+K_Met*x(1)*P)^2;
end
r1num = kj(1)/((P*x(4)+1e-06)^2.5)*(P^2*x(1)*x(2)-(P*x(4))^3*P*x(3)/Keq_1);
%rden = (1+K_CO*P*x(3)+K_H2*P*x(4)+K_Met*x(1)*P+K_H20*P*x(2)/(P*x(4)+1e-06))^2;
r1 = r1num/rden;%kmol/(kgcat*h)
r2num = kj(2)/(P*x(4)+1e-06)*(P^2*x(3)*x(2)-P^2*x(4)*x(5)/Keq_2);
r2 = r2num/rden;%kmol/(kgcat*h)
r3num = kj(3)/(P*x(4)+1e-06)^3.5*(P^3*x(1)*x(2)^2-P^5*x(4)^4*x(5)/Keq_3);
r3 = r3num/rden;%kmol/(kgcat*h)
rj=[r1 r2 r3]/3600; %kmol/(kgcat*s)
%%EQUAZIONI
dH0_r=delhr(nu,matrix_cp,delHo,T);
if t==3
    eff = eff*0.6; % efficienza dopo 3 anni
else
end
dwidz=(rho_bed*sum(nu.*rj'.*eff.*PM))/G;
dTdz=(-rho_bed*sum(dH0_r'.*rj*eff)+(Q*Sv))/(G*cp_mix);
dPdz = (-f*v^2*rhomix/(dp))/10^5;
T_pelle = T+Q*(d_tubi_in*log(d_tubi_ext/d_tubi_in)/(2*kt_tubo)+1/hcoeff);
figure(5)
title('Temperatura di pelle')
plot(z,T_pelle,'r*');hold on;
figure(4)
title('Conducibilità gas W/(m*K)')
plot(z,kt_miscela_2,'g*'); hold on;
plot(z,kt_miscela_1,'r*');hold on;
plot(z,kt_miscela,'k*');hold on;
dydt=[dwidz dTdz dPdz]';
end

function output=cond_mix(specie_cond,T,x,PM)
cond = zeros(1,length(x));
for i=1:length(specie_cond)
    C1=specie_cond(:,1);
    C2=specie_cond(:,2);
    C3=specie_cond(:,3);
    C4=specie_cond(:,4);

    cond(i)=(C1(i).*T.^C2(i))./(1+(C3(i)./T)+(C4(i)./(T^2)));%W/mK
end
output=(sum(x.*cond(i).*PM.^(1/3)))/(sum(x.*PM.^(1/3)));
end

function output=viscosity_mix(specie_viscosity,T,x,PM)
visc = zeros(1,length(x));
for i=1:length(specie_viscosity)
C1=specie_viscosity(:,1);
C2=specie_viscosity(:,2);
C3=specie_viscosity(:,3);
C4=specie_viscosity(:,4);

visc(i)=(C1(i).*T.^C2(i))./(1+(C3(i)./T)+(C4(i)./(T^2)));%W/mK
end
output=(sum(x.*visc(i).*PM.^(1/3)))/(sum(x.*PM.^(1/3)));
end

function out=delhr(nu, cpmatrix, delHo,T)
ns = length(delHo);
ent = delHo;
for i=1:ns
    cps = [cpmatrix(i,:) 0 0 0 0 0 0 0 0 0 0 6000 6200];
    ent(i) = delh(cps,delHo(i),T,'PERRY');
end
nreas = size(nu,1);
delHr = zeros(nreas,1);
for i=1:nreas
    delHr(i) = sum(nu(i,:).*ent);
end
out = delHr;
end
function output=delh(species,delHof,Tend,fcorr)
% per dubbi/suggerimenti/bug mandatemi una mail a fcapelizza@gmail.com
% fcorr è opzionale come input
% delh: entalpia di GAS IDEALE in j/kmol
% species: vettore con le costanti per le correlazioni dei cp, sono
% supportati tre range di lavoro diversi per ciascuna correlazione, la
% prima correlazione può essere o come quella del perry (esponenziale) o come quella del
% NIST (polinomiale), le costanti da C1 a C5 sono quelle relative alla
% prima correlazione, i dati da inserire sono quelli dati dal testo o dal
% NIST, non vanno manipolate, la funzione fa già tutto.
% La seconda e la terza correlazione sono sempre del tipo NIST, sono
% richiesti i coefficienti dalla A alla E.
% gli ultimi due elementi del vettore species (da 17 elementi) sono le
% temperature massime di lavoro della prima e della seconda correlazione,
% se non ci sono 3 correlazioni sul NIST (mi sembra ci siano solo per i gas
% biatomici) mettete al posto dei 5 coefficienti della terza correlazione
% degli 0, in realtà potete mettere quello che volete, non cambia nulla.
% delHof: l'entalpia della specie a 298.15 K in j/kmol
% Tend: temperatura finale di integrazione in K
% fcorr: STRINGA che può essere o 'PERRY' o 'NIST' che specifica il tipo di correlazione
% usata nella prima integrazione, se non mettete niente usa quella del
% perry
if nargin<4
    fcorr = 'PERRY';
else
    % do nothing
end
C1 = species(1);
C2 = species(2);
C3 = species(3);
C4 = species(4);
C5 = species(5);
A =species(6);
B =species(7);
C =species(8);
D = species(9);
E = species(10);
Tmax1 = species(16);
Tmax2 = species(17);
A2 =species(11);
B2 =species(12);
C_2 =species(13);
D2 = species(14);
E2 = species(15);
Tref = 298.15;
sp = strcmp('PERRY',fcorr);
sn = strcmp('NIST',fcorr);
if sn==0 && sp==0
    error('problema in correlazione, riprova con "NIST" o "PERRY".')
else
    %do nothing
end
if sp==1
    if Tend > Tmax1 && Tend<Tmax2
        intcp1 = C1*Tmax1 - C1*Tref - (2*C2*C3)/(exp(-(2*C3)/Tmax1) - 1) - (2*C4*C5)/(exp(-(2*C5)/Tmax1) + 1) + (2*C2*C3)/(exp(-(2*C3)/Tref) - 1) + (2*C4*C5)/(exp(-(2*C5)/Tref) + 1);
        intcp2 = ((B*Tend.^2)/2000 - (B*Tmax1.^2)/2000+ (C*Tend.^3)/3000000 - (C*Tmax1.^3)/3000000 + (D*Tend.^4)/4000000000 - (D*Tmax1.^4)/4000000000 - (1000000*E)./Tend + (1000000*E)./Tmax1 + A*Tend - A*Tmax1)*10^3;
        intcp = intcp1+intcp2;
    elseif Tend > Tmax1 && Tend>Tmax2
        intcp1 = C1*Tmax1 - C1*Tref - (2*C2*C3)/(exp(-(2*C3)/Tmax1) - 1) - (2*C4*C5)/(exp(-(2*C5)/Tmax1) + 1) + (2*C2*C3)/(exp(-(2*C3)/Tref) - 1) + (2*C4*C5)/(exp(-(2*C5)/Tref) + 1);
        intcp2 = ((B*Tmax2.^2)/2000 - (B*Tmax1.^2)/2000+ (C*Tmax2.^3)/3000000 - (C*Tmax1.^3)/3000000 + (D*Tmax2.^4)/4000000000 - (D*Tmax1.^4)/4000000000 - (1000000*E)./Tmax2 + (1000000*E)./Tmax1 + A*Tmax2 - A*Tmax1)*10^3;
        intcp3 =((B2*Tend.^2)/2000 - (B2*Tmax2.^2)/2000+ (C_2*Tend.^3)/3000000 - (C_2*Tmax2.^3)/3000000 + (D2*Tend.^4)/4000000000 - (D2*Tmax2.^4)/4000000000 - (1000000*E2)./Tend + (1000000*E2)/Tmax2 + A2*Tend - A2*Tmax2)*10^3;
        intcp = intcp1+intcp2+intcp3;
    else
        intcp1 = C1*Tend - C1*Tref - (2*C2*C3)/(exp(-(2*C3)/Tend) - 1) - (2*C4*C5)/(exp(-(2*C5)/Tend) + 1) + (2*C2*C3)/(exp(-(2*C3)/Tref) - 1) + (2*C4*C5)/(exp(-(2*C5)/Tref) + 1);
        intcp = intcp1;
    end
elseif sn==1
    if Tend > Tmax1 && Tend<Tmax2
        intcp1 = ((C2*Tmax1.^2)/2000 + (C3*Tmax1.^3)/3000000 - (1000000*C5)/Tmax1 + (C4*Tmax1.^4)/4000000000 - (C2*Tref.^2)/2000 - (C3*Tref.^3)/3000000 + (1000000*C5)/Tref - (C4*Tref.^4)/4000000000 + C1*Tmax1 - C1*Tref)*10^3;
        intcp2 = ((B*Tend.^2)/2000 - (B*Tmax1.^2)/2000+ (C*Tend.^3)/3000000 - (C*Tmax1.^3)/3000000 + (D*Tend.^4)/4000000000 - (D*Tmax1.^4)/4000000000 - (1000000*E)/Tend + (1000000*E)/Tmax1 + A*Tend - A*Tmax1)*10^3;
        intcp = intcp1+intcp2;
    elseif Tend > Tmax1 && Tend>Tmax2
        intcp1 = ((C2*Tmax1.^2)/2000 + (C3*Tmax1.^3)/3000000 - (1000000*C5)/Tmax1 + (C4*Tmax1.^4)/4000000000 - (C2*Tref.^2)/2000 - (C3*Tref.^3)/3000000 + (1000000*C5)/Tref - (C4*Tref^4)/4000000000 + C1*Tmax1 - C1*Tref)*10^3;
        intcp2 = ((B*Tmax2.^2)/2000 - (B*Tmax1.^2)/2000+ (C*Tmax2.^3)/3000000 - (C*Tmax1.^3)/3000000 + (D*Tmax2.^4)/4000000000 - (D*Tmax1.^4)/4000000000 - (1000000*E)/Tmax2 + (1000000*E)/Tmax1 + A*Tmax2 - A*Tmax1)*10^3;
        intcp3 =((B2*Tend.^2)/2000 - (B2*Tmax2.^2)/2000+ (C_2*Tend.^3)/3000000 - (C_2*Tmax2.^3)/3000000 + (D2*Tend.^4)/4000000000 - (D2*Tmax2.^4)/4000000000 - (1000000*E2)/Tend + (1000000*E2)/Tmax2 + A2*Tend - A2*Tmax2)*10^3;
        intcp = intcp1+intcp2+intcp3;
    else
        intcp1 = ((C2*Tend.^2)/2000 + (C3*Tend.^3)/3000000 - (1000000*C5)/Tend + (C4*Tend.^4)/4000000000 - (C2*Tref.^2)/2000 - (C3*Tref.^3)/3000000 + (1000000*C5)/Tref - (C4*Tref.^4)/4000000000 + C1*Tend - C1*Tref)*10^3;
        intcp = intcp1;
    end
end
output = delHof + intcp;
end
function out =rhomixgp(PM,xs,P,T)
% rhomass in kg/m^3
% PMs in kg/kmol
% P in Pa
% T in K
R = 8314.5;
rhom = P/(R*T);
rhoms = ones(1,length(PM))*rhom;
rhomass = rhoms.*xs.*PM;
out = sum(rhomass);
end

function out=cpmix(cp_matrix,xs,T)
sz = size(cp_matrix);
cps = zeros(sz(2),1);
for i=1:sz(1)
    cps(i) = cp_matrix(i,1)+cp_matrix(i,2)*(cp_matrix(i,3)/T/sinh(cp_matrix(i,3)/T))^2+cp_matrix(i,4)*(cp_matrix(i,5)/T/cosh(cp_matrix(i,5)/T))^2;
end
out = sum(cps'.*xs);%J/kmolK
end

function out=delg(species,delGo,delHof,T,fcorr,npoints)
            % per dubbi/suggerimenti/bug mandatemi una mail a fcapelizza@gmail.com
            % fcorr e npoints sono argomenti opzionali
            % delg: energia libera di gibbs di una specie in fase GAS IDEALE a T in j/kmol
            % species: vettore con le costanti per le correlazioni dei cp, sono
            % supportati tre range di lavoro diversi per ciascuna correlazione, la
            % prima correlazione può essere o come quella del perry (esponenziale) o come quella del
            % NIST (polinomiale), le costanti da C1 a C5 sono quelle relative alla
            % prima correlazione, i dati da inserire sono quelli dati dal testo o dal
            % NIST, non vanno manipolate, la funzione fa già tutto.
            % La seconda e la terza correlazione sono sempre del tipo NIST, sono
            % richiesti i coefficienti dalla A alla E.
            % gli ultimi due elementi del vettore species (da 17 elementi) sono le
            % temperature massime di lavoro della prima e della seconda correlazione,
            % se non ci sono 3 correlazioni sul NIST (mi sembra ci siano solo per i gas
            % biatomici) mettete al posto dei 5 coefficienti della terza correlazione
            % degli 0, in realtà potete mettere quello che volete, non cambia nulla.
            % delGo: energia libera di gibbs della specie a 298.15 K
            % delHof: l'entalpia della specie a 298.15 K in j/kmol
            % T: temperatura finale di integrazione in K
            % fcorr: STRINGA che può essere o 'PERRY' o 'NIST' che specifica il tipo di correlazione
            % usata nella prima integrazione, se non mettete niente usa quella del
            % perry
            % npoints: numero di elementi usati per l'integrazione numerica della
            % funzione delh, se non impostato usa 200 punti di integrazione

            % BASI TEORICHE: integro l'equazione di Van't Hoff su un range di T facendo
            % variare l'intervallo per l'integrazione fatta da delh.
            
            %VALIDAZIONE: tramite diagrammi di francis per metano, etano, propano,
            %etilene e acetilene.

            %ATTENZIONE: se provate a vedere come varia delg in funzione della
            %temperatura vedrete che per i composti con delHof<0 sale e poi riscende,
            %con un andamento parabolico, mentre per composti elementari (delHof=0)
            %questo scende dopo i 298.15 K, è normale, è dato dal fatto che
            %d(delg/(R*T))/dT =-delH/(R*T^2) , se delH(T) è positivo delg/(R*T)
            %diminuisce col variare della temperatura.

            if nargin<6
                npoints = 200;
            elseif nargin<5
                npoints = 200;
                fcorr = 'PERRY';
            elseif nargin<4
                error('numero di input non sufficienti')
            else
                %do nothing
            end
            Tref = 298.15;
            R = 8314.5;
            Xs = linspace(Tref,T,npoints);
            Ys = zeros(1,npoints);
            th = thermo_funs;
            for i=1:length(Ys)
                Ys(i) = th.delh(species,delHof,Xs(i),fcorr)./(R*Xs(i).^2);
            end
            val_int = trapz(Xs,Ys);
            out = (delGo/(R*Tref) - val_int)*R*T;
        end
