function [comb_lat_coefficients, comb_long_coefficients, comb_lat_coefficients_table, comb_long_coefficients_table] = MF52_combined_fitting(lat_coefficients,long_coefficients, longitudinal_data)

lat_coeff = lat_coefficients; % Pure Lateral Coefficients from Prev Fit
long_coeff = long_coefficients; % Pure Longitudinal Coefficients from Prev Fit

load(longitudinal_data);

disp('Combined Slip Fitting....')
disp(tireid)
disp(testid)

P = round(P .*0.145038);
IA = round(IA);
SA = round(SA);
V = round(V);

m=1:length(SL);     % point counter
sl_spline = spline(m,SL-0.02); 

% Remember to use actual SL graph, not SR as it is used for CALSPAN machine
% control

z = fnzeros(sl_spline);  % location of zero crossings
z=round(z(1,:));   % save only whole numbers
z(1:4:length(z))=[];  % Every 4th point is removed

% Subset the data
% create empty arrays, vertcat each new data set to it. 
sa = [];
sl = [];
fz = [];
fx = [];
fy = [];
ia = [];
p = [];
v = [];
temp_centre = [];
temp_inner = [];
temp_outer =[];

for n=1:3:length(z)-2 % Vertical Concatenation of the data points between zero slip angle points
 
 sa_1 = SA(z(n):z(n+2)); % reference slip angle used for csaps functions
 sa = [sa ; sa_1];

 sl_1 = SL(z(n):z(n+2)); % reference slip ratios used for csaps functions, change to percentage
 sl = [sl ; sl_1];

 fz_1= FZ(z(n):z(n+2)); % Collect set of FZ values within this range of zero points
 fz_2 = mean(fz_1).*ones(length(fz_1),1); % collect set of mean values for this range of FZ
 fz = [fz ; fz_2]; % Append to array with clean fz values for fitting

 fx_1= FX(z(n):z(n+2)); % Collect set of FX values within this range of zero points
 fx_1 = mean(fz_2).*fx_1./fz_1; % each fy value is scaled slightly by deviation from mean fz
 fx = [fx; fx_1];

 fy_1= FY(z(n):z(n+2)); % Collect set of FY values within this range of zero points
 fy_1 = mean(fz_2).*fy_1./fz_1; % each fy value is scaled slightly by deviation from mean fz
 fy = [fy; fy_1];

 ia = [ia ; round(IA(z(n):z(n+2)))]; % round values to obtain whole numbers

 p_1 = mean(P(z(n):z(n+2))).*ones(length(z(n):z(n+2)),1);
 p = [p ; p_1]; % convert to psi, round values to obtain whole numbers

 v_1 = mean(V(z(n):z(n+2))).*ones(length(z(n):z(n+2)),1); % Create groups for Velocities
 v = [v ; v_1]; % convert to psi, round values to obtain whole numbers

 temp_centre = [temp_centre; TSTC(z(n):z(n+2))]; % Middle Temperature in degrees
 temp_inner = [temp_inner; TSTI(z(n):z(n+2))]; % Inner Temperature in degrees
 temp_outer = [temp_outer; TSTO(z(n):z(n+2))]; % Outer Temperature in degrees

end 

% Horcat arrays

t(:,1) = sl;
t(:,2) = round(sa);
t(:,3) = fz;
t(:,4) = ia;
t(:,5) = p;
t(:,6) = fy;
t(:,7) = fx;
t(:,10) = temp_inner;
t(:,11) = temp_centre;
t(:,12) = temp_outer;
t(:,13) = round(v,-1);

t = sortrows(t,[4,5,3]); % Sort the table based on camber, pressure and fz

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
ia4_p14 = t(t(:,4) == 4 & t(:,5) == 14,:);

disp('Camber and Pressure Groups Created....')

% Pressure-based bin for MF5.2
comb_p10 = [ia0_p10 ; ia2_p10 ; ia4_p10];
comb_p12 = [ia0_p12 ; ia2_p12 ; ia4_p12];
comb_p14 = [ia0_p14; ia2_p14; ia4_p14]; 

save('comb_p10.mat','comb_p10');
save('comb_p12.mat','comb_p12');
save('comb_p14.mat','comb_p14');

%% Combined Lateral 

% Initial guess for combined lateral coefficients

RBY1=15.5;
RBY2=15.13;
RBY3=-0.001377;
RCY1=0.9231;
REY1=-0.561;
REY2=-0.06803;
RHY1=0.00928;
RHY2=0.01527;
RVY1=0;
RVY2=0;
RVY3=0;
RVY4=0;
RVY5=0;
RVY6=0;

INPUT = [comb_p10(:,1),comb_p10(:,2),comb_p10(:,3),comb_p10(:,4)]; % [SL, SA, FZ, IA]

global FZ0
FZ0 = abs(mean(comb_p10(:,3)));

C_str = ["RBY1"	"RBY2"	"RBY3"	"RCY1"	"REY1"	"REY2"	"RHY1"	"RHY2"	"RVY1"	"RVY2"	"RVY3"	"RVY4"	"RVY5"	"RVY6"];
C_old = [RBY1	RBY2	RBY3	RCY1	REY1	REY2	RHY1	RHY2	RVY1	RVY2	RVY3	RVY4	RVY5	RVY6]';

options = optimoptions("fmincon",'EnableFeasibilityMode',true,'Algorithm','sqp','MaxFunctionEvaluations',300000,'MaxIterations',300000);

fy_fun = @(x) sum((comb_p10(:,6) - PAC2002_FY_Combined(x,lat_coeff,INPUT)).^2); % Combined Lateral Optimisation Function

y0 = C_old;

y = fmincon(fy_fun,y0,[],[],[],[],[],[],@(x) nonlcon_FY_Combined(x,lat_coeff,INPUT),options);

disp('Combined Lateral Tyre Modelling Complete...')

%% Combined Longitudinal

% Initial Guess for Combined Longitudinal Coefficients

RBX1    =   12.35; % Pacejka Representative Parameters
RBX2    =   -10.77;
RCX1    =   1.092;
REX1    =   1.188;
REX2    =   0.006;
RHX1    =   0.007;

D_str = ["RBX1"	"RBX2"	"RCX1"	"REX1"	"REX2"	"RHX1"];
D_old = [RBX1	RBX2	RCX1	REX1	REX2	RHX1]';

A = [];
b = [];

fun = @(x) sum((comb_p10(:,7) - PAC2002_FX_Combined(x,long_coeff,INPUT)).^2); % Optimisation Function

z0 = D_old;

z = fmincon(fun,z0,A,b,[],[],[],[],@(x) nonlcon_FX_Combined(x,long_coeff,INPUT),options);

disp('Combined Longitudinal Tyre Modelling Complete...')

%% Results

comb_lat_coefficients = y;
comb_long_coefficients = z;

comb_lat_coefficients_table = [C_str',num2cell(y)];
comb_long_coefficients_table = [D_str',num2cell(z)];

test_fy = PAC2002_FY_Combined(y,lat_coeff,INPUT);
test_fx = PAC2002_FX_Combined(z,long_coeff,INPUT);

figure
plot3(comb_p10(:,1),comb_p10(:,2),comb_p10(:,6),'r.','DisplayName','Measured')
hold on
plot3(INPUT(:,1),INPUT(:,2),test_fy,'b.','DisplayName','MF52 Fit')
hold off
legend Location northeast
xlabel('SL (-)')
ylabel('SA (deg)')
zlabel('FY (N)')
title('Combined Lateral Force Fitting Comparison')
subtitle([tireid,' ',testid])

figure
plot3(comb_p10(:,1),comb_p10(:,2),comb_p10(:,7),'r.','DisplayName','Measured')
hold on
plot3(INPUT(:,1),INPUT(:,2),test_fx,'b.','DisplayName','MF52 Fit')
hold off
legend Location northeast
xlabel('SL (-)')
ylabel('SA (deg)')
zlabel('FX (N)')
title('Combined Longitudinal Force Fitting Comparison')
subtitle([tireid,' ',testid])

end