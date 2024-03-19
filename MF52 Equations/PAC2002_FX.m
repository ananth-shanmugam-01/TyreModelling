function [FX] = PAC2002_FX(B,X)
% PAC2002_FX Pure Longitudinal Slip Fit
% C is Input Coefficients
% X is Input Conditions [SL, SA, FZ, IA]

global FZ0
global LFZO

LHX = 1;
LVX = 1;
LMUX = 1;
LGX = 1;
LCX = 1;
LEX = 1;
LKX = 1;

PCX1 = B(1);
PDX1 = B(2);
PDX2 = B(3);
PDX3 = B(4);
PEX1 = B(5);
PEX2 = B(6);
PEX3 = B(7);
PEX4 = B(8);
PKX1 = B(9);
PKX2 = B(10);
PKX3 = B(11);
PHX1 = B(12);
PHX2 = B(13);
PVX1 = B(14);
PVX2 = B(15);

KAPPA = X(:,1); % Slip Ratio [-]
ALPHA = X(:,2)*pi/180; % degrees to [Rad]
FZ = abs(X(:,3)); % [N]
GAMMA = X(:,4).*pi/180; % degrees to [Rad]

FZ0PR = FZ0.* LFZO; %15, NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE
DFZ = (FZ - FZ0PR ) ./ FZ0PR ; %14, (%30)

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

FX = FX0;

end

