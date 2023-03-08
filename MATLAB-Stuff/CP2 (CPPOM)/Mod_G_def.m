% ====================================================================
% == Course:     Chemical Plants II                                 ==
% == Title:      Spray dryer                                        ==
% == Solution:   Cocurrent configuration                            ==
% == Author:     Alessandro Di Pretoro (SuPER,POLIMI)               ==
% ==             alessandro.dipretoro@polimi.it                     ==
% ==             www.super.chem.polimi.it/ing-alessandro-di-pretoro ==
% == Date:       2019.12.02                                         ==
% ====================================================================

clc
clear all
close all

global P Gdry A B C rhoL Kpw nAvp h MWwater MWair dHev CpL CpG rhoG muG dp mp0 vg mdry Pw P0

% Physical Properties
Tp0=303;										% K
Tg0=403;										% K
P=1;											% atm
CpL=1;											% kcal/kg/K
CpG=0.25;										% kcal/kg/K
MWwater=18;                                     % kg/kmol
MWair=29;                                       % kg/kmol
R=0.0821;                                       % m^3*atm/K/kmol
dHev=540;										% kcal/kg
rhoL=1000;										% kg/m^3
rhoG=P/R/Tg0*MWair;                             % kg/m^3
muG=2.3*1e-5;                                   % kg/m/s
kG=8*1.0e-6;									% kcal/m/s/K
Diff=1.8*1e-5;                               	% m^2/s

% Antoine equation
A=18.3036;
B=3816.44;
C=-46.13;

% Process stream
Qmilk=1750/3600;                            	% kg/h -> kg/s
fat=4.76/100;                           		% kg/kg
Qfat=Qmilk*fat;                     			% kg/s
Ql=Qmilk-Qfat;                  				% kg/s
Win=Ql/Qfat;									% kg/kg
Wout=0.005; 									% kg/kg
dp=2*1e-4;     									% m
vp0=0.3;										% m/s
mp0=rhoL*3.14/6*(dp^3);         				% kg
mdry=mp0/(1+Win);                   			% kg
mout=mdry*(1+Wout);                     		% kg
nAvp=Qmilk/mp0;                             	% drops/s

% Gas stream
Gdry=20;										% kg/s
D=5.5;     										% m

% Momentum Balance
vg=Gdry/rhoG/3.14/D^2*4;                        % m/s
vs=(rhoL-rhoG)*dp^2*9.81/18/ muG;               % m/s
vp=vg+vs;                                       % m/s
vs0=vp0-vg;                                     % m/s

% Mass and heat transfer
Re=rhoG*vs*dp/muG;                              % Reynolds
Pr=muG*CpG /kG; 								% Prandtl
Sc= muG/rhoG/Diff;      						% Schmidt
Nu=2+0.4*Re^0.5*Pr^(1/3);               		% Nusselt
Sh=2+0.4*Re^0.5*Sc^(1/3);                   	% Sherwood
h=Nu*kG/dp;                             		% kcal/m^2/s/K
Kc=Sh*Diff/dp;                      			% m/s
Kpw=Kc/R/Tg0*MWwater;           				% kg/m^2/s/atm

% Integration
tspan=linspace(0,10);
y0=[mp0 0 Tp0 Tg0 0];

[t,y]=ode89(@solver,tspan,y0);

% Plots
figure(1)
plot(y(:,5),y(:,1))
xlabel('Length [m]')
ylabel('Particle mass [m]')

figure(2)
plot(y(:,5),y(:,3))
hold on
plot(y(:,5),y(:,4))
xlabel('Length [m]')
ylabel('Temperature [K]')
legend('Particle','Air')

figure(3)
plot(y(:,5),-(0.33333333333333333333333333333333.*9.81.*y(:,1).*(rhoG - 1.0*rhoL))./(dp.*muG.*pi.*rhoL))
hold on
plot(y(:,5),vg-(0.33333333333333333333333333333333.*9.81.*y(:,1).*(rhoG - 1.0*rhoL))./(dp.*muG.*pi.*rhoL))
plot(y(:,5),vg*ones(1*length(tspan)))
xlabel('Length [m]')
ylabel('Velocity [m/s]')
legend('vs','vp','vg')

figure(4)
plot(tspan,y(:,5))
xlabel('Time [s]')
ylabel('Axial coordinate [m]')


figure(5)
plot(t,y(:,1))
xlabel('time [s]')
ylabel('Particle mass [m]')

figure(6)
plot(t,y(:,3))
hold on
plot(y(:,5),y(:,4))
xlabel('time [s]')
ylabel('Temperature [K]')
legend('Particle','Air')

figure(7)
plot(t,-(0.33333333333333333333333333333333.*9.81.*y(:,1).*(rhoG - 1.0*rhoL))./(dp.*muG.*pi.*rhoL))
hold on
plot(t,vg-(0.33333333333333333333333333333333.*9.81.*y(:,1).*(rhoG - 1.0*rhoL))./(dp.*muG.*pi.*rhoL))
plot(t,vg*ones(1*length(tspan)))
xlabel('time [s]')
ylabel('Velocity [m/s]')
legend('vs','vp','vg')

figure(8)
plot(t,MWair/MWwater*y(:,2)*P/Gdry)
xlabel('time [s]')
ylabel('Pw [atm]')

figure(9)
plot(t,exp(A-B./(y(:,3)+C))./760)
xlabel('time [s]')
ylabel('P0 [atm]')

figure(10)
plot(t,exp(A-B./(y(:,3)+C))./760-MWair./MWwater.*y(:,2).*P./Gdry)
xlabel('time [s]')
ylabel('P0 [atm]')

function spraydryer=solver(t,x)

global P Gdry A B C rhoL Kpw nAvp h MWwater MWair dHev CpL CpG rhoG muG dp mp0 vg mdry Pw P0

mp=x(1);
Gvap=x(2);
Tp=x(3);
Tg=x(4);
g = 9.81;
pig = pi;
vs=(rhoL-rhoG)*dp^2*g/(18*muG);
%-(0.33333333333333333333333333333333*g*mp*(rhoG - 1.0*rhoL))/(dp*muG*pig*rhoL);
z=x(5);

Pw=MWair/MWwater*Gvap*P/Gdry;
P0=exp(A-B/(Tp+C))/760;                                                     % mmHg--->atm
Sp=3.14*((mp/rhoL*6/3.14)^(1/3))^2;

spraydryer(1)=heaviside(mp-mdry*1.00005)*Kpw*Sp*(Pw-P0);
spraydryer(2)=-heaviside(mp-mdry*1.00005)*Kpw*Sp*(Pw-P0)*nAvp;
spraydryer(3)=(h*Sp*(Tg-Tp)+heaviside(mp-mdry*1.00005)*Kpw*Sp*(Pw-P0)*dHev)/mp/CpL;
spraydryer(4)=-h*Sp*(Tg-Tp)*nAvp/(Gvap+Gdry)/CpG;
spraydryer(5)=(vs+vg);

spraydryer=spraydryer';
end