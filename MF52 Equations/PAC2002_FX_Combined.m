function [FX_Combined] = PAC2002_FX_Combined(combined_coeff,long_coeff,X)

% PAC2002_FX Pure Longitudinal Slip Fit
% C is Input Coefficients
% X is Input Conditions [SL, SA, FZ, IA]

KAPPA = X(:,1); % Slip Ratio [-]
ALPHA = X(:,2)*pi/180; % degrees to [Rad]
FZ = abs(X(:,3)); % [N]
GAMMA = X(:,4).*pi/180; % degrees to [Rad]

global FZ0
global LFZO

FZ0PR = FZ0 * LFZO; %15, NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE
DFZ = (FZ-FZ0PR ) ./ FZ0PR ; %14, (%30)

LHX = 1;
LVX = 1;
LMUX = 1;
LGX = 1;
LCX = 1;
LEX = 1;
LKX = 1;
LXAL = 1;

% Pure Longitudinal Coefficients
% Taken from previous fitting process

PCX1=long_coeff(1);
PDX1=long_coeff(2);
PDX2=long_coeff(3);
PDX3=long_coeff(4);
PEX1=long_coeff(5);
PEX2=long_coeff(6);
PEX3=long_coeff(7);
PEX4=long_coeff(8);
PKX1=long_coeff(9);
PKX2=long_coeff(10);
PKX3=long_coeff(11);
PHX1=long_coeff(12);
PHX2=long_coeff(13);
PVX1=long_coeff(14);
PVX2=long_coeff(15);

% Longitudinal Force (Combined Slip) 

RBX1 = combined_coeff(1);
RBX2 = combined_coeff(2);
RCX1 = combined_coeff(3);
REX1 = combined_coeff(4);
REX2 = combined_coeff(5);
RHX1 = combined_coeff(6);

%% Longitudinal Code

% Pure Longitudinal Code
SHX = (PHX1 + PHX2.*DFZ).*LHX;
SVX = FZ.*(PVX1 + PVX2.*DFZ).*LVX.*LMUX;
KAPPAX = KAPPA + SHX;
GAMMAX = GAMMA.*LGX;
CX = PCX1.*LCX;
MUX = (PDX1 + PDX2.*DFZ).*(1 - PDX3.*(GAMMA.^2)).*LMUX;
DX = MUX.*FZ;
EX = (PEX1 + PEX2.*DFZ + PEX3.*(DFZ.^2)).*(1 - PEX4.*sign(KAPPAX)).*LEX;
KX = FZ.*(PKX1 + PKX2.*DFZ).*exp(PKX3.*DFZ).*LKX; % Longitudinal Slip Stiffness
BX = KX./(CX.*DX);
FX0 = DX.*sin( (CX.*atan(BX.*KAPPAX - EX.*(BX.*KAPPAX - atan(BX.*KAPPAX)))) + SVX);

% Combined Longitudinal Equations

SHXAL = RHX1;
CXAL = RCX1;
BXAL = RBX1.*cos( atan(RBX2.*KAPPA)).*LXAL; % cos term will always be positive regardless of slip ratio direction
ALPHAS = ALPHA + SHXAL;
EXAL = REX1 + REX2.*DFZ;
GXAL = ( cos(CXAL.*atan(BXAL.*ALPHAS - EXAL.*(BXAL.*ALPHAS - atan(BXAL.*ALPHAS))))) ./ ( cos(CXAL.*atan(BXAL.*SHXAL - EXAL.*(BXAL.*SHXAL - atan(BXAL.*SHXAL)))));
FX_C = GXAL.*FX0;
FX_Combined = FX_C;

end

