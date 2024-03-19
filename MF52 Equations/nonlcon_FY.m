function [c,ceq] = nonlcon_FY(A,X)
% Non Linear Constraints for PAC2002 FY
% Pure Lateral Slip Fitting Code
%  Input C - Pure Lateral Coefficients
%  Input X - Input Load Conditions (SL, SA, FZ, IA)
% SL not used for pure lateral slip, kept for ease of use.

global FZ0
global LFZO

LGAY = 1;
LHY = 1;
LVY = 1;
LCY = 1;
LEY = 1;
LMUY = 1;
LKY = 1;

PCY1 = A(1);
PDY1 = A(2);
PDY2 = A(3);
PDY3 = A(4);
PEY1 = A(5);
PEY2 = A(6);
PEY3 = A(7);
PEY4 = A(8);
PKY1 = A(9);
PKY2 = A(10);
PKY3 = A(11);
PHY1 = A(12);
PHY2 = A(13);
PHY3 = A(14);
PVY1 = A(15);
PVY2 = A(16);
PVY3 = A(17);
PVY4 = A(18);

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

c(1) = double(any(EY > 1));
c(2) = abs(PEY3) - 1;
c(3) = abs(PEY4) - 1;
c(4) = abs(PVY2) - 1;
c(5) = abs(PVY4) - 1;
c(6) = abs(PVY1) - 1;
c(7) = abs(PVY2) - 1;

ceq = [];

end