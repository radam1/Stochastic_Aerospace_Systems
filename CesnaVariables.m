% Filename: CesnaVariables.m
% 
% The purpose of this file is to load the state-space Cesna Variables into
% the main practical file. The values of this file are either drawn from
% the table in the practical instrucions or taken from tables 8.1, 8.2, and 8.3
% from the lecture notes. It does not add any stability augmentation to the
% system. 
%
% It is important to note that this code draws heavily from the "cit2a.m"
% file included in the lecture notes

% Center of Gravity
xcg = 0.30;

% Aircraft Weight and Velocity
W = 44675; % N
V = 51.4; % m/sec

% Mass and Altitude
m = 4556; % kg
h = 0; % m

% Wing and Air Density
S = 24.2; % m^2
rho = 1.225; % kg/m^3

% Mean Aerodynamic Chord and Other Lengths
c_bar = 2.022; % m
b = 13.36; % m
mu_c = 76;
mu_b = 11; 
lh = 5.5; % m
CL  = 1.1360; %From Flight Dynamics Lecture Notes

%Turbulence Parameters: 
Lg = 150; %From Instructions
B = b/(2 * Lg); %In order to use table 8.1 and 8.2 to find Iug, Iag and tau values
sigmaug = 3; %m/s
sigmaug_V = sigmaug/V;
sigmaag   = sigmaug/V;

sigmavg   = 2; %m/s
sigmabg   = sigmavg/V;

%Interpolating Table 8.1 to get moment values
Iug0 = ((0.0292262-0.0132418)/(0.05-0.03125) * (0.0445 - 0.03125) + 0.0132418)*sigmaug_V^2;
Iag0 =  (((0.0212969-0.0096887)/(0.05-0.03125) * (0.0445 - 0.03125)) + 0.0096887)*sigmaag^2;

%Interpolating Table 8.2 to get tau1,2,3 values
tau1 = ((0.106813-0.077782)/(0.05-0.03125) * (0.0445-0.03125) + 0.077782);
tau2 = ((0.569551-0.512936)/(0.05-0.03125) * (0.0445-0.03125) + 0.512936);
tau3 = ((0.427748-0.383390)/(0.05-0.03125) * (0.0445-0.03125) + 0.383390);

%Interpolating Table 8.3 to get tau4,5,6 values
tau4 = ((0.064521-0.047613)/(0.05-0.03125) * (0.0445-0.03125) + 0.047613);
tau5 = ((0.336211-0.310788)/(0.05-0.03125) * (0.0445-0.03125) + 0.310788);
tau6 = ((0.227862-0.214478)/(0.05-0.03125) * (0.0445-0.03125) + 0.214478);

% Stability Derivatives
KX2 = 0.012;
KZ2 = 0.037;
KXZ = 0.002;
KY2 = 0.980;

% Longitudinal Stability Derivatives
CX0 = 0;
CXu = -0.2173;
CXalpha = 0.4692;
CXq = 0;
CXdelta = 0;

CZ0 = -1.1360;
CZu = -2.2720;
CZalpha = -5.1300;
CZalphadot = -1.4050;
CZq = -3.8400;
CZdelta = -0.6238;

Cmu = 0;
Cmalpha = -0.4000;
Cmalphadot = -3.6150;
Cmq = -7.3500;
Cmdelta = -1.5530;

% Lateral Stability Derivatives
CYb = -0.9896;
CYp = -0.0870;
CYr = 0.4300;
CYd_a = 0;
CYd_r = 0.3037;

Clb = -0.0772;
Clp = -0.3415;
Clr = 0.2830;
Cld_a = -0.2349;
Cld_r = 0.0286;

Cnb = 0.1628;
Cnp = -0.0108;
Cnr = -0.1930;
Cnd_a = 0.0286;
Cnd_r = -0.1261;

CYfb = 0;
Clfb = 0;
Cnfb = 0;

%Stability Derivatives Derrived From Provided MATLAB Scripts(cit2a.m)
Clpw = 0.8*Clp;
Cnpw = 0.9*Cnp;
Clrw = 0.7*Clr;
Cnrw = 0.2*Cnr;

%From Equation Assumption in 6-1 of Flight Dynamics Notes
CYbetaDot = 0; 
CnbetaDot = 0; 

% Calculate Elements of The State Matrices
yb   = (V/b)*CYb/(2*mu_b);
yphi = (V/b)*CL/(2*mu_b);
yp   = (V/b)*CYp/(2*mu_b);
yr   = (V/b)*(CYr-4*mu_b)/(2*mu_b);
ybg  = yb;
ydr  = (V/b)*CYd_r/(2*mu_b);
den  = b*4*mu_b*(KX2*KZ2-KXZ^2)/V;
lb   = (Clb*KZ2+Cnb*KXZ)/den;
lp   = (Clp*KZ2+Cnp*KXZ)/den;
lr   = (Clr*KZ2+Cnr*KXZ)/den;
lda  = (Cld_a*KZ2+Cnd_a*KXZ)/den;
ldr  = (Cld_r*KZ2+Cnd_r*KXZ)/den;
lug  = (-Clrw*KZ2-Cnrw*KXZ)/den;
lbg  = lb;
lag  = (Clpw*KZ2+Cnpw*KXZ)/den;
nb   = (Clb*KXZ+Cnb*KX2)/den;
np   = (Clp*KXZ+Cnp*KX2)/den;
nr   = (Clr*KXZ+Cnr*KX2)/den;
nda  = (Cld_a*KXZ+Cnd_a*KX2)/den;
ndr  = (Cld_r*KXZ+Cnd_r*KX2)/den;
nug  = (-Clrw*KXZ-Cnrw*KX2)/den;
nbg  = nb;
nag  = (Clpw*KXZ+Cnpw*KX2)/den;
aug1 =-(V/Lg)^2*(1/(tau1*tau2));
aug2 =-(tau1+tau2)*(V/Lg)/(tau1*tau2);
aag1 =-(V/Lg)^2*(1/(tau4*tau5));
aag2 =-(tau4+tau5)*(V/Lg)/(tau4*tau5);
abg1 =-(V/Lg)^2;
abg2 =-2*(V/Lg);
bug1 = tau3*sqrt(Iug0*V/Lg)/(tau1*tau2);
bug2 = (1-tau3*(tau1+tau2)/(tau1*tau2))*sqrt(Iug0*(V/Lg)^3)/(tau1*tau2);
bag1 = tau6*sqrt(Iag0*V/Lg)/(tau4*tau5);
bag2 = (1-tau6*(tau4+tau5)/(tau4*tau5))*sqrt(Iag0*(V/Lg)^3)/(tau4*tau5);
bbg1 = sigmabg*sqrt(3*V/Lg);
bbg2 = (1-2*sqrt(3))*sigmabg*sqrt((V/Lg)^3);

% Create A and B(State-Space Matrices)
A = [yb yphi yp    yr 0    0    0    0    ybg  0;
     0  0    2*V/b 0  0    0    0    0    0    0;
     lb 0    lp    lr lug  0    lag  0    lbg  0;
     nb 0    np    nr nug  0    nag  0    nbg  0;
     0  0    0     0  0    1    0    0    0    0;
     0  0    0     0  aug1 aug2 0    0    0    0;
     0  0    0     0  0    0    0    1    0    0;
     0  0    0     0  0    0    aag1 aag2 0    0;
     0  0    0     0  0    0    0    0    0    1;
     0  0    0     0  0    0    0    0    abg1 abg2];

B = [0   ydr 0    0    0;
     0   0   0    0    0;
     lda ldr 0    0    0;
     nda ndr 0    0    0;
     0   0   bug1 0    0;
     0   0   bug2 0    0;
     0   0   0    bag1 0;
     0   0   0    bag2 0;
     0   0   0    0    bbg1;
     0   0   0    0    bbg2];
