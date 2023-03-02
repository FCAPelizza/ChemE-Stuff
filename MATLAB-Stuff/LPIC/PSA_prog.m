C_Co1 = [0,0,0.0007,0.0028,0.0063,0.0229,0.0556,0.0986,0.1653,0.2535,0.3458,];% []
C_Co2 = [0.4514,0.5611,0.6806,0.7743,0.8576,0.9236,0.9792,1.0000,1.0000]; % []
t1 = [0,8.9,9.0,9.2,9.4,9.6,9.8,10.0,10.2,10.4,10.6]; %h
t2 = [10.8,11,11.25,11.5,11.75,12.0,12.5,12.8,13]; %h
C_Co = [C_Co1,C_Co2];
t = [t1,t2]; % h
der1 = C_Co;
der1(1) = 0;

for i=1:(length(C_Co)-1)
    der1(i+1) = (C_Co(i+1)-C_Co(i))/(t(i+1)-t(i));
end
figure(2)
plot(t,der1)
der2 = der1;
for i=2:(length(C_Co)-1)
    der2(i+1) = (der1(i+1)-der1(i))/(t(i+1)-t(i));
end
figure(3) 
plot(t,der2)
t_star = t(13)+(t(14)-t(13))*1/3
inv_C_Co = 1-C_Co;
t_star_integral = trapz(t,inv_C_Co)
tb =t(6)+(t(7)-t(6))*0.82875
figure(1)
plot(t,C_Co,'k')
xline(t_star,'r')
xline(t_star_integral,'g')
xline(tb,'b')

A_PSA = 1175/300
LUB_over_L = (1-tb/t_star_integral)
% non uso il 10.63% del letto
rho_mol = 5.19972/13.3509*1000;
v = 300*10^3/rho_mol; %m/h
C0 = 273.86; %mol/m3
MM = 0.002; %kg/mol
capacity = v*C0*MM*t_star_integral
rhobed = 498; %kgbed/m^3
mueq = 0.186; %kgH2/kgbed
mu0 = 0.01; %kgH2/kgbed
L = capacity/(rhobed*(mueq-mu0))
LUB = L*LUB_over_L


