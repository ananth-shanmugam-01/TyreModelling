function [long_coefficients, fx_table] = MF52_longitudinal_fitting(longitudinal_data)

load(longitudinal_data);

disp('Drive/Brake (Pure Longitudinal Slip)....')

P = round(P .*0.145038);
IA = round(IA);
SA = round(SA);
V = round(V);

disp(['Selected Tire - ', tireid])
disp(['Camber Range - ', num2str(unique(round(IA))'), ' degrees'])
disp(['Pressure Range - ', num2str(unique(round(P))'), ' Psi'])

m=1:length(SL);     % point counter
sl_spline = spline(m,SL-0.02); 

% Remember to use actual SL graph, not SR as it is used for CALSPAN machine
% control

z = fnzeros(sl_spline);  % location of zero crossings
z=round(z(1,:));   % save only whole numbers
z(1:4:length(z))=[];  % Every 4th point is removed

% figure
% hold on
% fnplt(sl_spline)
% plot(z,zeros(length(z)),'o')
% hold off

% Sorting condition data into values divided by zero points

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
% Pure Longitudinal Slip, 0 Slip Angle
% Camber Angles of 0, 2, 4 degrees
% Pressures of 10, 12, 14 psi

ia0_p10 = t(t(:,2) == 0 & t(:,4) == 0 & t(:,5) == 10,:);
ia0_p12 = t(t(:,2) == 0 & t(:,4) == 0 & t(:,5) == 12,:);
ia0_p14 = t(t(:,2) == 0 & t(:,4) == 0 & t(:,5) == 14,:);

ia2_p10 = t(t(:,2) == 0 & t(:,4) == 2 & t(:,5) == 10,:);
ia2_p12 = t(t(:,2) == 0 & t(:,4) == 2 & t(:,5) == 12,:);
ia2_p14 = t(t(:,2) == 0 & t(:,4) == 2 & t(:,5) == 14,:);

ia4_p10 = t(t(:,2) == 0 & t(:,4) == 4 & t(:,5) == 10,:);
ia4_p12 = t(t(:,2) == 0 & t(:,4) == 4 & t(:,5) == 12,:);
ia4_p14 = t(t(:,2) == 0 & t(:,4) == 4 & t(:,5) == 14,:);

disp('Camber and Pressure Groups Created....')

% Pressure-based bin for MF5.2
long_p10 = [ia0_p10 ; ia2_p10 ; ia4_p10];
long_p12 = [ia0_p12 ; ia2_p12 ; ia4_p12];
long_p14 = [ia0_p14; ia2_p14; ia4_p14]; 

save('pure_long_p10.mat','long_p10');
save('pure_long_p12.mat','long_p12');
save('pure_long_p14.mat','long_p14');

% Yet to create the UI and option to select pressure value
% Current work is for 10psi, although any other pressure can be analysed
% with a bit of inconvenience

disp('Evaluating 10 Psi...')

INPUT = [long_p10(:,1),zeros(height(long_p10),1),long_p10(:,3),long_p10(:,4)]; % [SL, SA, FZ, IA]

global FZ0 R0
FZ0     =  abs(mean(INPUT(:,3))); % = FNOMIN = 'nominal wheel load' % Trial with ABS value
R0      =   0.2032; % Loaded radius - 0.2032m

% Initial Guess for Longitudinal 

PCX1 = 1.685; %Shape factor Cfx for longitudinal force
PDX1 = 1.210; %Longitudinal friction Mux at Fznom
PDX2 = -0.037 ;%Variation of friction Mux with load
PDX3 = 0; %Variation of friction Mux with camber
PEX1 = 0.344; %Longitudinal curvature Efx at Fznom
PEX2 = 0.095; %Variation of curvature Efx with load
PEX3 = -0.002; %Variation of curvature Efx with load squared
PEX4 = 0; %Factor in curvature Efx while driving
PKX1 = 21.51; %Longitudinal slip stiffness Kfx/Fz at Fznom
PKX2 = -0.163 ;%Variation of slip stiffness Kfx/Fz with load
PKX3 = 0.245; %Exponent in slip stiffness Kfx/Fz with load
PHX1 = -0.002; %Horizontal shift Shx at Fznom
PHX2 = 0.002; %Variation of shift Shx with load
PVX1 = 0;%Vertical shift Svx/Fz at Fznom
PVX2 = 0 ;%Variation of shift Svx/Fz with load
  
clear x B
% List of Starting Longitudinal Force Coefficients
B_str ={'PCX1' 'PDX1' 'PDX2' 'PDX3' 'PEX1' 'PEX2' 'PEX3' 'PEX4' 'PKX1' 'PKX2' 'PKX3' 'PHX1' 'PHX2' 'PVX1' 'PVX2'};
B_old =[PCX1 PDX1 PDX2 PDX3 PEX1 PEX2 PEX3 PEX4 PKX1 PKX2 PKX3 PHX1 PHX2 PVX1 PVX2];

options =optimset('MaxFunEvals',2000000,'MaxIter',2000000,'Display','final','TolX',1e-7,'TolFun',1e-7);

fun = @(x) sum((long_p10(:,7) - PAC2002_FX(x,INPUT)).^2); % MSE
x0 = B_old;

tic
x = fmincon(fun,x0,[],[],[],[],[],[],@(x) nonlcon_FX(x,INPUT),options);
toc

% Outputs

long_coefficients = x;
fx_table = [B_str',num2cell(x)'];

test_fx = PAC2002_FX(x,INPUT); % Pure FX Test

figure
plot3(long_p10(:,1),long_p10(:,4),long_p10(:,7),'r.','DisplayName','Measured')
hold on
plot3(INPUT(:,1),INPUT(:,4),test_fx,'b.','DisplayName','MF52')
hold off
legend
xlabel('SL')
ylabel('IA (deg)')
zlabel('FX')
title('Longitudinal Force Fitting')
subtitle([tireid,' ',testid])

end