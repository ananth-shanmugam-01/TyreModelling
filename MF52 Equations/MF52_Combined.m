function [FY,FX] = MF52_Combined(A,X)
% Input - [SL, SA, FZ, IA] 
% Output - [FY, FX]
% Uses Coefficients Determined by Individual Lateral and Longitudinal
% Fitting Process

% Inputs

global FZ0

LFZO = 1.2; % Changed locally from 1.1 due to poor fit at loads above 1100N
% Poor extrapolation as expected, but safeguard during use case

KAPPA = X(:,1); % Slip Ratio [-]
ALPHA = X(:,2)*pi/180; % degrees to [Rad]
FZ = abs(X(:,3)); % [N]
GAMMA = X(:,4).*pi/180; % degrees to [Rad]

FZ0PR = FZ0 * LFZO; %15, NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE
DFZ = (FZ-FZ0PR ) ./ FZ0PR ; %14, (%30)

% Scaling Factors

LGAY    =   1;
LHY     =   1;
LVY     =   1;
LCY     =   1;
LEY     =   1;
LMUY    =   0.8; % Changed to Fit
LHX     =   1;
LVX     =   1;
LMUX    =   1;
LGX     =   1;
LCX     =   1;
LEX     =   1;
LKX     =   1;
LXAL    =   1;

% Pure Lateral Coefficients
PCY1	=	A(1);
PDY1	=	A(2);
PDY2	=	A(3);
PDY3	=	A(4);
PEY1	=	A(5);
PEY2	=	A(6);
PEY3	=	A(7);
PEY4	=	A(8);
PKY1	=	A(9);
PKY2	=	A(10);
PKY3	=	A(11);
PHY1	=	A(12);
PHY2	=	A(13);
PHY3	=	A(14);
PVY1	=	A(15);
PVY2	=	A(16);
PVY3	=	A(17);
PVY4	=	A(18);

% Pure Longitudinal Coefficients
PCX1	=	A(19);
PDX1	=	A(20);
PDX2	=	A(21);
PDX3	=	A(22);
PEX1	=	A(23);
PEX2	=	A(24);
PEX3	=	A(25);
PEX4	=	A(26);
PKX1	=	A(27);
PKX2	=	A(28);
PKX3	=	A(29);
PHX1	=	A(30);
PHX2	=	A(31);
PVX1	=	A(32);
PVX2	=	A(33);

% Combined Lateral Coefficients
RBY1	=	A(34);
RBY2	=	A(35);
RBY3	=	A(36);
RCY1	=	A(37);
REY1	=	A(38);
REY2	=	A(39);
RHY1	=	A(40);
RHY2	=	A(41);
RVY1	=	A(42);
RVY2	=	A(43);
RVY3	=	A(44);
RVY4	=	A(45);
RVY5	=	A(46);
RVY6	=	A(47);

% Combined Longitudinal Coefficients
RBX1	=	A(48);
RBX2	=	A(49);
RCX1	=	A(50);
REX1	=	A(51);
REX2	=	A(52);
RHX1	=	A(53);


% Lateral Force

GAMMAY = GAMMA.*LGAY;
SHY = (PHY1+ PHY2.*DFZ).*LHY + PHY3.*GAMMAY;
SVY = FZ.*((PVY1 + PVY2.*DFZ).*LVY + (PVY3 + PVY4.*DFZ)).*LMUY;
ALPHAY = ALPHA + SHY;
CY = PCY1.*LCY;
MUY = (PDY1 + PDY2.*DFZ).*(1 - PDY3.*(GAMMAY.^2)).*LMUY;
DY = MUY.*FZ;
EY = (PEY1 + PEY2.*DFZ).*(1 - (PEY3 + PEY4.*GAMMAY).*sign(ALPHAY)).*LEY;
KY0 = PKY1.*FZ0.*sin( 2.*atan(FZ./(PKY2.*FZ0PR))); % Cornering Stiffness
KVY0 = PHY3.*KY0 + FZ.*(PVY3 + PVY4.*DFZ); % Camber Stiffness
KY = KY0.*(1 - PKY3.*abs(GAMMAY));
BY = KY./(CY.*DY);
FY0 = DY.*sin(CY.*atan(BY.*ALPHAY - EY.*(BY.*ALPHAY - atan(BY.*ALPHAY)))) + SVY;
FY = FY0;

% Combined Lateral Force

BYK = RBY1.*cos(atan(RBY2.*(ALPHA - RBY3)));
CYK = RCY1;
EYK = REY1 + REY2.*DFZ;
SHYK = RHY1 + RHY2.*DFZ; 
KS = KAPPA + SHYK;
DVYK = MUY.*FZ.*(RVY1 + RVY2.*DFZ + RVY3.*GAMMA).*cos(atan(RVY4.*ALPHA));
SVYK = DVYK.*sin(RVY5.*atan(RVY6.*KAPPA));
DYK = FY0./( cos(CYK.*atan(BYK.*SHYK - EYK.*(BYK.*SHYK - atan(BYK.*SHYK)))));
GYK = ( cos(CYK.*atan( BYK.*KS - EYK.*(BYK.*KS - atan(BYK.*KS))))) ./ ( cos(CYK.*atan(BYK.*SHYK - EYK.*(BYK.*SHYK - atan(BYK.*SHYK)))));
FY_C = GYK.*FY + SVYK;
FY = FY_C;

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
FX0 = DX.*sin( (CX.*atan(BX.*KAPPAX - EX.*(BX.*KAPPAX - atan(BX.*KAPPAX)))) + SVX); % Pure Longitudinal Force

% Combined Longitudinal Equations

SHXAL = RHX1;
CXAL = RCX1;
BXAL = RBX1.*cos( atan(RBX2.*KAPPA)).*LXAL; 
ALPHAS = ALPHA + SHXAL;
EXAL = REX1 + REX2.*DFZ;
GXAL = ( cos(CXAL.*atan(BXAL.*ALPHAS - EXAL.*(BXAL.*ALPHAS - atan(BXAL.*ALPHAS))))) ./ ( cos(CXAL.*atan(BXAL.*SHXAL - EXAL.*(BXAL.*SHXAL - atan(BXAL.*SHXAL)))));
FX_C = GXAL.*FX0;
FX = FX_C;

end

