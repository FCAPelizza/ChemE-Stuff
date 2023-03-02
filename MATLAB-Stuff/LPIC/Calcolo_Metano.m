clear
clc
%[CH4 CO2 O2 H20 N2]
global delGo delHo properties n0_CH4 T_fiamma
methane = [-0.703029 108.4773 -42.52157 5.862788  0.678565 85.81217 11.26467 -2.114146 0.138190 -26.42221 0 0 0 0 0 1300 6000];
CO2 = [24.99735 55.18696 -33.69137 7.948387 -0.136638 58.16639 2.720074 -0.492289 0.038844 -6.447293 0 0 0 0 0 1200 6000];
Water = [30.09200 6.832514 6.793435 -2.534480	 0.082139 41.96426 8.622053 -1.499780 0.098119 -11.15764 0 0 0 0 0 1700 6000];
Oxygen=[31.32234 -20.23531 57.86644 -36.50624 -0.007374 30.03235 8.772972 -3.988133 0.788313 -0.741599 23.245 10.72071 -2.020498 0.146449 9.245722 700 2000];
Nitrogen = [28.98641 1.853978  -9.647459 16.63537 0.000117 19.50583 19.88705 -8.598535 1.369784 0.527601 35.51872 1.128728 -0.196103 0.014662 -4.553760 500 2000];
delGo=1e+07.*[-5.0496, -39.4370,0, -22.8590,0];
delHo=1e+07.*[-7.4520, -39.3514,0,-24.1814,0];
properties = [methane;CO2;Oxygen;Water;Nitrogen];
beta=linspace(0,1,50);
alfa=linspace(0.75,2,50);
T_fiamma=zeros(50,50);
del = zeros(50,50);
n0_CH4 = 1;
for i=1:length(alfa)
    for j=1:length(beta)
        if j==1
            options = optimset('TolX',10^-10);
            T_fiamma(i,j) = fzero(@(T)fun_to_delTfiamma(T,alfa(i),beta(j)),2250,options);
            del(i,j) = fun_to_delTfiamma(T_fiamma(i,j),alfa(i),beta(j));
        else
            options = optimset('TolX',10^-10);
            T_fiamma(i,j) = fzero(@(T)fun_to_delTfiamma(T,alfa(i),beta(j)),T_fiamma(i,j-1),options);
            del(i,j) = fun_to_delTfiamma(T_fiamma(i,j),alfa(i),beta(j));
        end
    end
end
T_fiamma_max = fzero(@(T)fun_to_delTfiamma(T,1,1),3000)
delta_T_adiabatico=T_fiamma;
for i=1:50
    if alfa(i)<=1
        figure(1)
    plot(beta,delta_T_adiabatico(i,:),'k');hold on
    elseif alfa(i)>1
        plot(beta,delta_T_adiabatico(i,:),'r');hold on
    end
end
hold off

figure(2)
surf(beta,alfa,delta_T_adiabatico)
ylabel('alfa');
xlabel('beta');
zlabel('T adiabatica di fiamma')

n_fuel = zeros(50,50);
delta= zeros(50,50);
for i=1:length(alfa)
    for j=1:length(beta)
        n_fuel(i,j) = fzero(@(k)fun_for_fuel(alfa(i),beta(j),T_fiamma(i,j),k),2);
        delta(i,j) = fun_for_fuel(alfa(i),beta(j),T_fiamma(i,j),n_fuel(i,j));
    end
end
m_fuel = n_fuel*16.04*3600;
P_fuel = m_fuel*1.400; %price per hour
n_oxy = (beta.*alfa').*n_fuel*2*3600;
P_O2 = n_oxy.*32*0.50;
Ptot = P_fuel+P_O2;
price_min = min(min(Ptot))
[ialfa, ibeta] = find(Ptot==price_min);
cond_opt = [alfa(ialfa), beta(ibeta)]
figure(6)
surf(beta,alfa,Ptot)

figure(3)
surf(beta,alfa,m_fuel)
ylabel('alfa');
xlabel('beta');
zlabel('massa metano da immettere in kg/h')



function out=fun_to_delTfiamma(T_fiamma,alfa,beta)
global delHo properties n0_CH4

nu=[-1 1 -2 2 0]; %[CH4 CO2 O2 H20 N2]
%T_in_fuel=672.039;%K
T_in_oxydizer = 298.15; % K
T_in_fuel = 298.15;
nO2 = alfa*2*n0_CH4; % calcolo quantità moli ossigeno
n_oxy = beta*nO2; % moli ossigeno ceh arrivano dall'ossigeno puro
nO2_air = nO2-n_oxy; % moli ossigeno che arrivano dall'aria
n0_N2 = nO2_air/0.21*0.79; % moli azoto
n0=[n0_CH4 0 nO2 0 n0_N2]; % vettore moli iniziali
lambda=min([n0_CH4,nO2/2]); % determinazione grado di avanzamento 
n_out=n0+nu*lambda;
hi_in = zeros(1,5);
hi_out = zeros(1,5);
for i=1:5
    if i==1
        hi_in(i) = delH(properties(i,:),delHo(i),T_in_fuel,'NIST'); %entalpia in entrata per componente
    else
        hi_in(i) = delH(properties(i,:),delHo(i),T_in_oxydizer,'NIST');
    end
end
Hin = sum(hi_in.*n0); %entalpia totale in entrata
for i=1:5
    hi_out(i) = delH(properties(i,:),delHo(i),T_fiamma,'NIST');% entalpia in entrata per componente
end
Hout = sum(hi_out.*n_out); %entalpia totale in uscita
BE = Hin - Hout; % bilancio energetico
out = BE;
end


function out = fun_for_fuel(alfa,beta,T_flame,k)
global delHo properties n0_CH4
nu=[-1 1 -2 2 0]; %[CH4 CO2 O2 H20 N2]
T_out=1.1613e+03; %K
n_fuel=n0_CH4*k;
nO2 = alfa*2*n_fuel;
n_oxy = beta*nO2;
nO2_air = (1-beta)*nO2;
n0_N2 = nO2_air/0.21*0.79;
n0=[n_fuel 0 nO2 0 n0_N2];
lambda=min([n_fuel,nO2/2]);
n_out=n0+nu*lambda;
hi_in =zeros(1,5);
hi_out = zeros(1,5);
for i =1:5
    hi_in(i) = delH(properties(i,:),delHo(i),T_flame,'NIST');
    hi_out(i) = delH(properties(i,:),delHo(i),T_out,'NIST');
end
Hin = sum(hi_in.*n_out);
Hout = sum(hi_out.*n_out);
Qcalc = Hin-Hout;
Qwanted = 18.8438 * 10^6;
out = Qwanted-Qcalc;
end

function output=delH(species,delHof,Tend,fcorr)
% da C1 a C5 sono i coefficienti per la prima correlazione, che può essere
% o perry (con cosh e sinh) o nist (polinomiale), se volete la perry
% mettete come fcorr la stringa 'PERRY', tutto in joule/(kmol*k), le nist sono
% in j/mol, sono da mettere come presenti sul sito, non adattatele,
% l'output della funzione poi viene in j/kmol, poi da A a E sono i
% coefficienti della seconda correlazione sul nist, mentre da A2 a E2 sono
% i coefficienti della terza correlazione, utile solo se si ha gas
% biatomici che hanno 3 correlazioni diverse nel nist, Tmax1 è la
% temperatura massima del range di lavoro della prima correlazione, la
% Tmax2 è la temperatura massima del range di lavoro della seconda
% correlazione, se non c'è una terza correlazione mettete come Tmax2 6000K
% e mettete tutti i coefficienti della terza correlazione come 0, non li
% userà comunque. Il vettore species è [C1 C2 C3 C4 C5 A B C D E A2 B2 C_2
% D2 E2. La T è sempre in K. delHof è l'entalpia di formazione a 298.15K in
% joule/kmol.
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