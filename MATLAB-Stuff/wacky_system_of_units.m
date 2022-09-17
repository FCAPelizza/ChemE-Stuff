clear all
% syms kg m K s
% mu = kg/(m*s);
% kt = kg*m/(s^3*K);
% D = m^2/s;
% P = kg/m/s^2;
% time = mu/P
% temperature = P*D/kt
% length = (D*mu/P)^0.5
% mass = mu*(D*mu/P)^0.5*mu/P
% speed = length/time % (D*P/mu)^0.5/(mu/P)
syms mu kt D P
time = mu/P;
temperature = P*D/kt;
length = (D*mu/P)^0.5;
mass = mu*(D*mu/P)^0.5*mu/P;
speed = length/time; % 
acceleration = speed/time;
force = acceleration * mass;
energy = force * length;
power = energy / time;
pretty(power)
%pretty(speed)