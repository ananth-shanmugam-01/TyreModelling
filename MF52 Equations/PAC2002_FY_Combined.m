function [FY_Combined, FY] = PAC2002_FY_Combined(combined_coeff,lat_coeff,X)
% Combined FY 

global FZ0
global LFZO

LGAY = 1;
LHY = 1;
LVY = 1;
LCY = 1;
LEY = 1;
LMUY = 0.8; % Changed to Fit
LKY = 1;

% Pure Lateral Slip Coefficients
% Taken for separate fitting process conducted for pure slip

PCY1=lat_coeff(1);
PDY1=lat_coeff(2);
PDY2=lat_coeff(3);
PDY3=lat_coeff(4);
PEY1=lat_coeff(5);
PEY2=lat_coeff(6);
PEY3=lat_coeff(7);
PEY4=lat_coeff(8);
PKY1=lat_coeff(9);
PKY2=lat_coeff(10);
PKY3=lat_coeff(11);
PHY1=lat_coeff(12);
PHY2=lat_coeff(13);
PHY3=lat_coeff(14);
PVY1=lat_coeff(15);
PVY2=lat_coeff(16);
PVY3=lat_coeff(17);
PVY4=lat_coeff(18);

% Lateral Force Coefficients (Combined Slip)
RBY1	=	combined_coeff(1);
RBY2	=	combined_coeff(2);
RBY3	=	combined_coeff(3);
RCY1	=	combined_coeff(4);
REY1	=	combined_coeff(5);
REY2	=	combined_coeff(6);
RHY1	=	combined_coeff(7);
RHY2	=	combined_coeff(8);
RVY1	=	combined_coeff(9);
RVY2	=	combined_coeff(10);
RVY3	=	combined_coeff(11);
RVY4	=	combined_coeff(12);
RVY5	=	combined_coeff(13);
RVY6	=	combined_coeff(14);

% C is Input Coefficients
% X is Input Conditions [SL, SA, FZ, IA]

% Pure Lateral Equations

KAPPA = X(:,1); % Slip Ratio [-]
ALPHA = X(:,2)*pi/180; % degrees to [Rad]
FZ = abs(X(:,3)); % [N]
GAMMA = X(:,4).*pi/180; % degrees to [Rad]

FZ0PR = FZ0 * LFZO; %15, NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE
DFZ = (FZ-FZ0PR ) ./ FZ0PR ; %14, (%30)

GAMMAY = GAMMA.*LGAY;
SHY = (PHY1+ PHY2.*DFZ).*LHY + PHY3.*GAMMAY;
SVY = FZ.*((PVY1 + PVY2.*DFZ).*LVY + (PVY3 + PVY4.*DFZ)).*LMUY;
ALPHAY = ALPHA + SHY;
CY = PCY1.*LCY;
MUY = (PDY1 + PDY2.*DFZ).*(1 - PDY3.*(GAMMAY.^2)).*LMUY;
DY = MUY.*FZ;
EY = (PEY1 + PEY2.*DFZ).*(1 - (PEY3 + PEY4.*GAMMAY).*sign(ALPHAY)).*LEY;

KY0 = PKY1.*FZ0.*sin( 2.*atan(FZ./(PKY2.*FZ0PR)))*LKY; % Cornering Stiffness
KVY0 = PHY3.*KY0 + FZ.*(PVY3 + PVY4.*DFZ); % Camber Stiffness
KY = KY0.*(1 - PKY3.*abs(GAMMAY));
BY = KY./(CY.*DY);
FY0 = DY.*sin(CY.*atan(BY.*ALPHAY - EY.*(BY.*ALPHAY - atan(BY.*ALPHAY)))) + SVY;
FY = FY0;

% Combined Lateral Force

BYK = RBY1.*cos(atan(RBY2.*(ALPHA - RBY3)));
CYK = RCY1;
EYK = REY1 + REY2.*DFZ;

SHYK = RHY1 + RHY2.*DFZ; % This kills 
KS = KAPPA + SHYK;
DVYK = MUY.*FZ.*(RVY1 + RVY2.*DFZ + RVY3.*GAMMA).*cos(atan(RVY4.*ALPHA));
SVYK = DVYK.*sin(RVY5.*atan(RVY6.*KAPPA));
DYK = FY0./( cos(CYK.*atan(BYK.*SHYK - EYK.*(BYK.*SHYK - atan(BYK.*SHYK)))));
GYK = ( cos(CYK.*atan( BYK.*KS - EYK.*(BYK.*KS - atan(BYK.*KS))))) ./ ( cos(CYK.*atan(BYK.*SHYK - EYK.*(BYK.*SHYK - atan(BYK.*SHYK)))));

FY_C = GYK.*FY + SVYK;

% FY = (DYK.*cos(CYK.*atan( BYK.*KS - EYK.*(BYK.*KS - atan(BYK.*KS))))) + SVYK;

FY_Combined = FY_C;

end