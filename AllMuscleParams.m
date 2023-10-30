clc
clear all
%% constants:

tstep = 0.001;
T_fin = 10;

%% Kinematic Boundaries
ThetaInit = 0*pi/180;
ThetaFin = 20*pi/180;
%% Skeletal Dynamics

% for QR problem, ThetaFin = 10*pi/180, T_fin = 10, Cload = -10, Kload =
% -10
% 

a0 = 0.286;
a1e = -0.014;
a2e = -3.96e-3;

Llarm = 0.5;
Marm = 2.5;
Iarm = Marm*(Llarm^2)/3;

%% Load constants

Mload = 5;
Cload = -10;
Kload = -10;
g = -9.8;

%%
c = 1.373e-4;
u = 52700;
k = 2.90;
q0 = 0.005;
m = 11.25;

width = 0.66;
arel = 0.41;
brel = 5.2;
qcrit = 0.03;
a = 1/(width^2);
brel = 5.2;


% CE force-velocity relationship
N    =   1.5; %[Fmax] eccentric force enhancement
K    =     5; %[] shape factor
%% Table 4.3B, pg 73, Kistemaker et al. 2006

% Selected muscle: MEF
lse_0 = 0.172; % Tendon natural length (m)
lce_opt = 0.092; % Optimal CE length (m)
lpe_0 = 1.4*lce_opt; % PE natural length (m)
%lpe_0 = 1*lce_opt; % PE natural length (m)
Fmax = 1420; % MVC value (N)

kse = Fmax/((0.04*lse_0)^2); % SE Stiffness (N/m)
kpe = (0.5*Fmax)/(1 + width - (lpe_0/lce_opt)); % PE Stiffness (N/m)

%% Muscle Kinematics & Kinetics initialization

lmtc_init = a0 + (a1e*ThetaInit) + (a2e*(ThetaInit^2));
%lmtc_init = lce_opt+lse_0;
lce_init = lmtc_init-lse_0;
lce_rel_init = lce_init/lce_opt;
