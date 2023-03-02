clear
clc
global F0 x_in mu delGo delHo cp_species P R
% vettore come [met et prop metoh CO2 CO H2 H2O]
% reazioni [SR metano; SR et; SR prop; WGS; form metoh]
methane = [0.3330e+05 0.7993e+05 2.0869e+03 0.4160e+05 991.96];
ethane = [0.4033e+05 1.3422e+05 1.6555e+03 0.7322e+05 752.87];
propane = [0.5192e+05 1.9245e+05 1.6265e+03 1.1680e+05 723.6];
methanol = [0.3925e+05 0.8790e+05 1.9165e+03 0.5365e+05 896.7];
hydrogen = [0.2762e+05 0.0956e+05 2.4660e+03 0.0376e+05 567.6];
CO = [0.2911e+05 0.0877e+05 3.0851e+03 0.0846e+05 1538.2];
CO2 = [0.2937e+05 0.3454e+05 1.4280e+03 0.2640e+05 588];
Water = [0.3336e+05 0.2679e+05 2.6105e+03 0.0890e+05 1169];
cp_species = [methane;ethane;propane;methanol;hydrogen;CO;CO2;Water]
x_in=[0.2022 0.0153 0.0049 0.0028 0.0053 0 0 0.7695];
mu = [-1,0,0,0,0,1,3,-1;0,-1,0,0,0,+2,+5,-2;0,0,-1,0,0,+3,+7,-3;0,0,0,0,+1,-1,+1,-1;0,0,0,+1,0,-1,-2,0];
delGo=1e+07*[-5.0496, -3.1926, -2.4390, -16.2320, -39.4370, -13.7150, 0, -22.8590];
delHo=1e+07.*[-7.5420, -8.3820, -10.4680, -20.0940, -39.3510, -11.0530, 0, -24.1814];
F0=1046,35; %kmol/h
P=20;%bar
primo_tent=[170 50 0.1 1000];

out=fsolve(@funtozero,primo_tent);
del = funtozero(out)
lambda = [out(1);F0*x_in(2)/-mu(2,2); F0*x_in(3)/-mu(3,3);out(2);out(3)];
Fout = F0 + sum(sum(lambda(1:5).*mu))
disp('[met et prop metoh CO2 CO H2 H2O]')
x_out=(F0.*x_in+sum(lambda(1:5).*mu))./Fout
Tout = out(4)


hi_in = zeros(1,8);
for i=1:8
    hi_in(i) = delH(cp_species(i,:),delHo(i),800);
end
hin = sum(hi_in.*x_in);
hi_out = zeros(1,8);
for i=1:8
    hi_out(i) = delH(cp_species(i,:),delHo(i),Tout);
end
hout = sum(hi_out.*x_out);

Q = Fout*hout-F0*hin;
Q_MW = Q*1e-06/3600

function out=funtozero(in)
global F0 x0 mu delGo delHo cp_species P x_in R
lambda = in(1:3)';
T = in(4);
lamSRet = F0*x_in(2)/-mu(2,2);
lamSRpr = F0*x_in(3)/-mu(3,3);
delGRT = zeros(1,8);
lambda_real = [lambda(1);lamSRet;lamSRpr;lambda(2:3)];
Fout = F0 + sum(sum(lambda_real.*mu));
x=(F0.*x_in+sum(lambda_real.*mu))./Fout;
for i = 1 : 8
    delGRT(i) = delgRT(cp_species(i,:),delGo(i),delHo(i),T);
end
delGrea_overRT = sum(delGRT.*mu,2);
Keq = exp(-delGrea_overRT);
Kp = prod((P*x).^mu,2);
out(1)=Keq(1)-Kp(1);
out(2)=Keq(4)-Kp(4);
out(3)=Keq(5)-Kp(5);
out(4)=((F0*x_in(1)-Fout*x(1))/(F0*x_in(1)))-0.85;
end
function output=delH(species,delHof,Tend)
C1 = species(1);
C2 = species(2);
C3 = species(3);
C4 = species(4);
C5 = species(5);
Tref = 298.15;
intcp = C1*Tend - C1*Tref - (2*C2*C3)/(exp(-(2*C3)/Tend) - 1) - (2*C4*C5)/(exp(-(2*C5)/Tend) + 1) + (2*C2*C3)/(exp(-(2*C3)/Tref) - 1) + (2*C4*C5)/(exp(-(2*C5)/Tref) + 1);
output = delHof + intcp;
end
function out=delgRT(species,delGo,delHof,T)
C1 = species(1);
C2 = species(2);
C3 = species(3);
C4 = species(4);
C5 = species(5);
Tref = 298.15;
R = 8314.5;
%intcp = @(Tend)C1*Tend - C1*Tref - (2*C2*C3)/(exp(-(2*C3)/Tend) - 1) - (2*C4*C5)/(exp(-(2*C5)/Tend) + 1) + (2*C2*C3)/(exp(-(2*C3)/Tref) - 1) + (2*C4*C5)/(exp(-(2*C5)/Tref) + 1);
h_over_RT2 = @(Tend) (delHof + C1*Tend - C1*Tref - (2*C2*C3)./(exp(-(2*C3)./Tend) - 1) - (2*C4*C5)./(exp(-(2*C5)./Tend) + 1) + (2*C2*C3)./(exp(-(2*C3)./Tref) - 1) + (2*C4*C5)./(exp(-(2*C5)./Tref) + 1))./(R*Tend.^2);
val_int = integral(h_over_RT2,Tref,T);
out = delGo/(R*Tref) - val_int;
end