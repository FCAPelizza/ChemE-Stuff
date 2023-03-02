m_wo_turbine= 1983; % kg/h
% 2110 kg/h esce se voglio avere come T out 300 F (circa 150 C)
m_w_turbine = 2416; % kg/h
deltam = m_w_turbine-m_wo_turbine
LHVmeth = 50.00; %Mj/kg
work_turb = 2.3774; % MW
eff_syst = work_turb/(deltam*LHVmeth/3600)
AGG = (m_w_turbine-m_wo_turbine)/m_wo_turbine