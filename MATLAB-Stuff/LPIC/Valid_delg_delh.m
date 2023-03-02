clear all
close all
clc
% vettore come [met et prop metoh CO2 CO H2 H2O O2 N2 C]
% reazioni [SR metano; SR et; SR prop; WGS; form metoh]
%methane = [0.3330e+05 0.7993e+05 2.0869e+03 0.4160e+05 991.96 85.81217 11.26467 -2.114146 0.138190 -26.42221 0 0 0 0 0 1500 6000];
ethane = [0.4033e+05 1.3422e+05 1.6555e+03 0.7322e+05 752.87 0 0 0 0 0 0 0 0 0 0 6000 0];
propane = [0.5192e+05 1.9245e+05 1.6265e+03 1.1680e+05 723.6 0 0 0 0 0 0 0 0 0 0 6000 0];
methanol = [0.3925e+05 0.8790e+05 1.9165e+03 0.5365e+05 896.7 0 0 0 0 0 0 0 0 0 0 6000 0];
hydrogen = [0.2762e+05 0.0956e+05 2.4660e+03 0.0376e+05 567.6 18.563083 12.257357 -2.859786 0.268238 1.977990 43.413560 -4.293079 1.272428 -0.096876 -20.533862 1500 2500];
CO = [0.2911e+05 0.0877e+05 3.0851e+03 0.0846e+05 1538.2 35.15070 1.300095 -0.205921 0.013550 -3.282780 0 0 0 0 0 1500 6000];
CO2 = [0.2937e+05 0.3454e+05 1.4280e+03 0.2640e+05 588 58.16639 2.720074 -0.492289 0.038844 -6.447293 0 0 0 0 0 1500 6000];
Water = [0.3336e+05 0.2679e+05 2.6105e+03 0.0890e+05 1169 41.96426 8.622053 -1.499780 0.098119 -11.15764 0 0 0 0 0 1500 6000];
%Oxygen=[0.29103e+05 0.1004e+05 2.5265e+03 0.09356e+05 1153.8 30.03235 8.772972 -3.988133 0.788313 -0.741599 23.245 10.72071 -2.020498 0.146449 9.245722 1500 2000];
%Nitrogen = [0.2911e+05 0.0861e+05 1.7016e+03 0.0010e+05 909.79 19.50583 19.88705 -8.598535 1.369784 0.527601 35.51872 1.128728 -0.196103 0.014662 -4.553760 1500 2000];
methane = [-0.703029 108.4773 -42.52157 5.862788  0.678565 85.81217 11.26467 -2.114146 0.138190 -26.42221 0 0 0 0 0 1300 6000];
%CO2 = [24.99735 55.18696 -33.69137 7.948387 -0.136638 58.16639 2.720074 -0.492289 0.038844 -6.447293 0 0 0 0 0 1200 6000];
%Water = [30.09200 6.832514 6.793435 -2.534480	 0.082139 41.96426 8.622053 -1.499780 0.098119 -11.15764 0 0 0 0 0 1700 6000];
Oxygen=[31.32234 -20.23531 57.86644 -36.50624 -0.007374 30.03235 8.772972 -3.988133 0.788313 -0.741599 23.245 10.72071 -2.020498 0.146449 9.245722 700 2000];
Nitrogen = [28.98641 1.853978  -9.647459 16.63537 0.000117 19.50583 19.88705 -8.598535 1.369784 0.527601 35.51872 1.128728 -0.196103 0.014662 -4.553760 500 2000];
Carbon = [21.17510, -0.812428, 0.448537, -0.043256, -0.013103, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,6000, 10000];
Acetylene = [40.68697, 40.73279 -16.17840 3.669741 -0.658411 67.47244 11.75110 -2.021470 0.136195 -9.806418 0 0 0 0 0 1100 6000];
ethylene = [-6.387880 184.4019 -112.9718 28.49593 0.315540 106.5104 13.73260 -2.628481 0.174595 -26.14469 0 0 0 0 0 1200 6000];
delGo=1e+07*[-5.0496, -3.1926, -2.4390, -16.2320, -39.4370, -13.7150, 0, -22.8590,0,0,0,21.0680,6.8440];
delHo=1e+07.*[-7.5420, -8.3820, -10.4680, -20.0940, -39.3510, -11.0530, 0, -24.1814,0,0,0,22.8200,5.2510];
Ts = linspace(298.15,1200,100);

g=zeros(1,100);
gCH4 = g; gH2=g; gfCH4 = g; gC=g; geth = g; gfeth = g; gact = g; gfact = g;
gethy = g; gfethy = g; gprop = g; gfprop = g; gmeoh = g; gfmeoh = g; gCO = g;
gfCO = g; gCO2 = g; gfCO2 = g; gw = g; gfw = g; gO2 = g;
th = thermo_funs;
gs = delg(methane,delGo(1),delHo(1),298.15,'NIST');
for i=1:length(Ts)
    gCH4(i) = delg(methane,delGo(1),delHo(1),Ts(i),'NIST');
    gact(i) = delg(Acetylene,delGo(12),delHo(12),Ts(i),'NIST');
    geth(i) = delg(ethane,delGo(2),delHo(2),Ts(i),'PERRY');
    gprop(i) = delg(propane,delGo(3),delHo(3),Ts(i),'PERRY');
    gC(i) = delg(Carbon,delGo(11),delHo(11),Ts(i),'NIST');
    gH2(i) = delg(hydrogen,delGo(7),delHo(7),Ts(i),'PERRY');
    gethy(i) = delg(ethylene,delGo(13),delHo(13),Ts(i),'NIST');
    gmeoh(i) = delg(methanol,delGo(4),delHo(4),Ts(i),'PERRY');
    gCO(i) = delg(CO,delGo(6),delHo(6),Ts(i),'PERRY');
    gCO2(i) = delg(CO2,delGo(5),delHo(5),Ts(i),'PERRY');
    gw(i) = th.delg(Water,delGo(8),delHo(8),Ts(i),'PERRY');
    gO2(i) = delg(Oxygen,0,0,Ts(i),'NIST');
    gfCH4(i) = (gCH4(i) - gC(i) - 2*gH2(i));
    gfeth(i) = (geth(i) - gC(i)*2 - 3*gH2(i));
    gfact(i) = (gact(i) - gC(i)*2 - gH2(i));
    gfethy(i) = (gethy(i) - gC(i)*2 - gH2(i)*2);
    gfprop(i) = (gprop(i)-gC(i)*3-gH2(i)*4);
    gfmeoh(i) = (gmeoh(i)-gC(i)-gO2(i)/2-2*gH2(i));
    gfCO(i) = (gCO(i)-gC(i)-gO2(i)/2);
    gfCO2(i) = (gCO(i)-gC(i)-gO2(i));
    gfw(i) = (gw(i)-gH2(i)-gO2(i)/2);
end

figure(1)
plot(Ts,gfCH4/10^6,'g','DisplayName','methane');grid on;hold on;
legend('methane')
plot(Ts,gfeth/2/10^6,'k','DisplayName','ethane');grid on;hold on;
plot(Ts,gfact/2/10^6,'r','DisplayName','acetylene');grid on;hold on;
plot(Ts,gfethy/2/10^6,'b','DisplayName','ethylene');grid on;hold on;
plot(Ts,gfprop/3/10^6,'m','DisplayName','propane');grid on;hold on;
plot(Ts,zeros(1,length(Ts)),'DisplayName','g=0 line');hold off;
title('Francis diagram in kj/(mol*nC)')
xlabel('T in K')
ylabel('delg in kj/(mol*nC)')


figure(2)
title('value of gibbs free energy of formation in a range of temperatures')
xlabel('T in K')
ylabel('delg of formation in j/kmol')
plot(Ts,gfCH4);hold on;
plot(Ts,gfeth);hold on;
plot(Ts,gfprop);hold on;
plot(Ts,gfethy);hold on;
plot(Ts,gfact);hold on;
plot(Ts,gfmeoh);hold on;
plot(Ts,gfCO,'k*');hold on;legend('CO');
plot(Ts,gfCO2);hold on;
plot(Ts,gfw);hold off; grid on;



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
for i=1:length(Ys)
    Ys(i) = delh(species,delHof,Xs(i),fcorr)./(R*Xs(i).^2);
end
val_int = trapz(Xs,Ys);
% figure(4)
%plot(Xs,Ys)
out = (delGo/(R*Tref) - val_int)*R*T;
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
