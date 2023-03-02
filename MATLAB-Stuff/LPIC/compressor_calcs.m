k = 1.30918;
Z = 0.995095;
MM = 39.3368;
R = 1544/MM;
beta = (20.42)^(1/4);
T = 286;
Ns = 300;
Ds = 1;
Q  = 24.85;
H = Z*R*T/((k-1)*MM)*k*(beta^(k/(k-1))-1)
N = Ns*H^(3/4)/Q
D = Ds*sqrt(Q)/(H^(1/4))
D_in_m = D*0.3048