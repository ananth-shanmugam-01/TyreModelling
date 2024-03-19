function [lat_coefficients, fy_table] = MF52_lateral_fitting(lateral_data)

load(lateral_data);

disp('Cornering (Free Rolling)....')
disp(['Selected Tire - ', tireid])
disp(['Camber Range - ', num2str(unique(round(IA))'), ' degrees'])

% Identify zero points to sort data

m=1:length(SA);     % point counter
sa_spline = spline(m,SA); 

z = fnzeros(sa_spline);  % location of zero crossings
z=round(z(1,:));   % save only whole numbers
z(1:4:length(z))=[];  % Every 4th point pertains to the end of negative slip sweep (return to zero SA)

% Optional Plot to Check Location of Zeros

% figure
% hold on
% plot(SA)
% plot(z,zeros(length(z),1),'o')
% hold off

% Sorting condition data into values divided by zero points
% This is useful to normalise CP force and moments

sl = [];
sa = [];
fz = [];
fy = [];
fy_sp = []; % Init array for csaps fit of fy-sa
fx = [];
mz = [];
mz_sp = []; % Init array for csaps fit of mz-sa
rl = [];
ia = [];
p = [];
temp_centre = []; % Init temperature array for middle surface temperature (TSTC)
temp_inner = []; % Init temperature array for inner surface temperature (TSTI)
temp_outer = []; % Init temperature array for output surface temperature (TSTO)

for n=1:3:length(z) % Vertical Concatenation of the data points between zero slip angle points
 
 sl_1 = SL(z(n):z(n+2)); % reference slip ratio used for csaps functions
 sl = [sl ; SL(z(n):z(n+2))];
    
 sa_1 = SA(z(n):z(n+2)); % reference slip angle used for csaps functions
 sa = [sa ; SA(z(n):z(n+2))];

 fz_1= FZ(z(n):z(n+2)); % Collect set of FZ values within this range of zero points
 fz_2 = mean(fz_1).*ones(length(fz_1),1); % collect set of mean values for this range of FZ
 fz = [fz ; fz_2]; % Append to array with clean fz values for fitting

 fx_1= FX(z(n):z(n+2)); % Collect set of FX values within this range of zero points
 fx_1 = mean(fz_2).*fx_1./fz_1; % each fx value is scaled slightly by deviation from mean fz
 fx = [fx; fx_1];

 fy_1= FY(z(n):z(n+2)); % Collect set of FY values within this range of zero points
 fy_1 = mean(fz_2).*fy_1./fz_1; % each fy value is scaled slightly by deviation from mean fz
 fy = [fy; fy_1];
 fy_sp = [fy_sp ; fnval(csaps(sa_1,fy_1,0.1),sa_1)]; % create csaps function and paste evaluated results

 mz_1= MZ(z(n):z(n+2)); % Repeat fy procedures for Mz
 mz_1 = mean(fz_2).*mz_1./fz_1; % each fy value is scaled slightly by deviation from mean fz
 mz = [mz; mz_1];
 mz_sp = [mz_sp ; fnval(csaps(sa_1,mz_1,0.1),sa_1)]; % create csaps function and paste evaluated results

 rl = [rl ; RL(z(n):z(n+2))]; % Loaded Radius

 ia = [ia ; round(IA(z(n):z(n+2)))]; % round values to obtain whole numbers

 p_1 = mean(P(z(n):z(n+2))).*ones(length(z(n):z(n+2)),1);
 p = [p ; round(p_1*0.145038)]; % convert to psi, round values to obtain whole numbers

 temp_centre = [temp_centre; TSTC(z(n):z(n+2))]; % Middle Temperature in degrees
 temp_inner = [temp_inner; TSTI(z(n):z(n+2))]; % Inner Temperature in degrees
 temp_outer = [temp_outer; TSTO(z(n):z(n+2))]; % Outer Temperature in degrees

end 

disp(['Pressure Range - ', num2str(unique(round(p))'), ' Psi'])

clear t

% INPUT = [SL ,SA, FZ, IA]

t(:,1) = sl;
t(:,2) = sa;
t(:,3) = fz;
t(:,4) = ia;
t(:,5) = p;
t(:,6) = fy;
t(:,7) = fx;
t(:,8) = mz;
t(:,9) = fy_sp;
t(:,10) = temp_inner;
t(:,11) = temp_centre;
t(:,12) = temp_outer;

t = sortrows(t,[4,5,3]); % Sort the table based on camber, pressure and fz

disp('Data Normalised and Sorted....')

% Create Bins for Inclination Angle and Pressure
% Camber Angles of 0, 2, 4 degrees
% Pressures of 10, 12, 14 psi

ia0_p10 = t(t(:,4) == 0 & t(:,5) == 10,:);
ia0_p12 = t(t(:,4) == 0 & t(:,5) == 12,:);
ia0_p14 = t(t(:,4) == 0 & t(:,5) == 14,:);

ia2_p10 = t(t(:,4) == 2 & t(:,5) == 10,:);
ia2_p12 = t(t(:,4) == 2 & t(:,5) == 12,:);
ia2_p14 = t(t(:,4) == 2 & t(:,5) == 14,:);

ia4_p10 = t(t(:,4) == 4 & t(:,5) == 10,:);
ia4_p12 = t(t(:,4) == 4 & t(:,5) == 12,:);
ia4_p14 = t(t(:,4) == 5 & t(:,5) == 14,:);

disp('Camber and Pressure Groups Created....')

% Pressure-based bin for MF5.2
lat_p10 = [ia0_p10 ; ia2_p10 ; ia4_p10];
lat_p12 = [ia0_p12 ; ia2_p12 ; ia4_p12];
lat_p14 = [ia0_p14; ia2_p14; ia4_p14]; 

save('pure_lat_p10.mat','lat_p10');
save('pure_lat_p12.mat','lat_p12');
save('pure_lat_p14.mat','lat_p14');

% Yet to create the UI and option to select pressure value
% Current work is for 10psi, although any other pressure can be analysed
% with a bit of inconvenience

disp('Evaluating 10 Psi...')

INPUT = [zeros(height(lat_p10),1),lat_p10(:,2),lat_p10(:,3),lat_p10(:,4)]; % [SL ,SA, FZ, IA]

global FZ0 R0;
FZ0     =  abs(mean(lat_p10(:,3))); % = FNOMIN = 'nominal wheel load'
R0      =   0.2032; % Loaded radius - 0.2032m for 16" & 0.2286m for 18" % Does not Affect Any Calculation for Lateral

% Perturbation Inputs Lateral Coefficients

PCY1 = 1.193;
PDY1 = -0.990;
PDY2 = 0.145;
PDY3 = -11.23;
PEY1 = -1.003;
PEY2 = -0.537;
PEY3 = -0.083;
PEY4 = -4.787;
PKY1 = -14.95;
PKY2 = 2.130;
PKY3 = -0.028;
PHY1 = 0.003;
PHY2 = -0.001;
PHY3 = 0.075;
PVY1 = 0.053;
PVY2 = -0.073;
PVY3 = -0.532;
PVY4 = 0.039;

% Starting Coefficients
A_str ={'PCY1' 'PDY1' 'PDY2' 'PDY3' 'PEY1' 'PEY2' 'PEY3' 'PEY4' 'PKY1' 'PKY2' 'PKY3' 'PHY1' 'PHY2' 'PHY3' 'PVY1' 'PVY2' 'PVY3' 'PVY4'};
A_old =[PCY1 PDY1 PDY2 PDY3 PEY1 PEY2 PEY3 PEY4 PKY1 PKY2 PKY3 PHY1 PHY2 PHY3 PVY1 PVY2 PVY3 PVY4];

options = optimset('MaxFunEvals',2000000,'MaxIter',2000000,'Display','final','Algorithm','interior-point','FinDiffType','central');
fun = @(x) sum((lat_p10(:,9) - PAC2002_FY(x,INPUT)).^2); % MSE
x0 = A_old;
A = [];
b = [];
[x, ~] = fmincon(fun,x0,A,b,[],[],[],[],@(x) nonlcon_FY(x,INPUT),options);

fy_coeff = [{'pressure'},{mean(lat_p10(:,5))};A_str',num2cell(x')];

% Outputs

lat_coefficients = x;
fy_table = [A_str',num2cell(x)'];

% Testing

test_0 = PAC2002_FY(x,INPUT);

figure
plot3(lat_p10(:,2),lat_p10(:,4),lat_p10(:,9),'r.','DisplayName','Measured')
hold on
plot3(lat_p10(:,2),lat_p10(:,4),test_0,'b.','DisplayName','MF52 Fit')
hold off
legend
xlabel('SA (deg)')
ylabel('IA (deg)')
zlabel('FY (N)')
title('Lateral Force Fitting Comparison')
subtitle([tireid,' ',testid])

disp('Cornering Tyre Modelling Complete...')

end
