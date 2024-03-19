
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code serves to interact with the individual (Lateral & Longitudinal)
% Pacejka MF52 Equations. Uses a function minimiser solver to obtain
% coefficients that can achieve a close fit between measured and simulated
% forces. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

addpath(genpath(pwd));

% Load Lateral Raw Data 
clc;
[file,folder]=uigetfile; % Select lateral data .mat file
lateral_data=fullfile(folder,file);

% Pure Lateral Modelling
global LFZO
LFZO = 1.1; % Empirical testing led to this value being selected for a fair fit
[lat_coefficients, fy_table] = MF52_lateral_fitting(lateral_data); % Results coefficients

% Load Longitudinal Raw Data 

[file,folder]=uigetfile;
longitudinal_data=fullfile(folder,file);

% Pure Longitudinal Modelling
[long_coefficients, fx_table] = MF52_longitudinal_fitting(longitudinal_data);

% Combined Slip Modelling
% Uses preceding results lat_coefficients, long_coefficients and combined
% slip data from longitudinal data set

[lat_comb_coefficients, long_comb_coefficients, comb_lat_coefficients_table, comb_long_coefficients_table] = MF52_combined_fitting(lat_coefficients, long_coefficients, longitudinal_data);

%% Final Coefficients

final_coefficients = [fy_table;fx_table;comb_lat_coefficients_table;comb_long_coefficients_table];

%% Load Segmented Raw Data

% From prior analysis, used to verify the results

load('comb_p10.mat')
load('pure_lat_p10.mat')
load('pure_long_p10.mat')

%% Coefficient Testing Plots

% Directory for arrays and their values
% t(:,1) = sl; % Slip Ratio
% t(:,2) = sa; % Slip Angle (deg)
% t(:,3) = fz; % Vertical Load (N)
% t(:,4) = ia; % Inclination Angle (deg)
% t(:,5) = p;  % Inflation Pressure (Psi)
% t(:,6) = fy; % Lateral Force (N)
% t(:,7) = fx; % Longitudinal Force (N)
% t(:,8) = mz; % Aligning Moment (Nm)
% t(:,9) = fy_sp; % Spline fit Lateral Force (N), use to improve fitting
% t(:,10) = temp_inner; % TSTI (Celcius)
% t(:,11) = temp_centre; % TSTC (Celcius)
% t(:,12) = temp_outer; % TSTO (Celcius)

sa_range = linspace(-12,12,50);
sr_range = linspace(-0.2,0.2,50);
[SA, SL] = meshgrid(sa_range,sr_range);
X = SA(:);
Y = SL(:);

y = round(lat_p10(:,3),-2)==-1100 & lat_p10(:,4) == 0;
z = round(long_p10(:,3),-2)==-1100 & long_p10(:,4) == 0;

test_input = [Y, X, -1100.*ones(length(X),1),zeros(length(X),1);...
    Y, X, -700.*ones(length(X),1),zeros(length(X),1);...
    Y, X, -400.*ones(length(X),1),zeros(length(X),1);];

[fy_test, fx_test] = MF52_Combined(str2double(final_coefficients(:,2)),test_input);

figure
plot3(test_input(:,2),test_input(:,1),fy_test,'.','DisplayName','MF52 Fit') % SL, SA, FY
hold on
plot3(comb_p10(:,2),comb_p10(:,1),comb_p10(:,6),'ro','DisplayName','Measured')
plot3(lat_p10(:,2),zeros(height(lat_p10),1),lat_p10(:,6),'bo','DisplayName','Pure Slip FY')
hold off
xlabel('SA')
ylabel('SL')
zlabel('FY')
title('Complete Lateral Modelling Analysis')
legend location northeast

figure
plot3(test_input(:,2),test_input(:,1),fx_test,'.','DisplayName','MF52 Fit')
hold on 
plot3(comb_p10(:,2),comb_p10(:,1),comb_p10(:,7),'ro','DisplayName','Combined Measured')
plot3(long_p10(:,2),long_p10(:,1),long_p10(:,7),'bo','DisplayName','Pure Slip FX')
hold off
xlabel('SA')
ylabel('SL')
zlabel('FX')
title('Complete Longitudinal Modelling Analysis')
legend location northeast

%% Friction Circle

[fy_sweep_1300,fx_sweep_1300] = MF52_Combined(str2double(final_coefficients(:,2)),[Y,X,-1300.*ones(height(X),1),zeros(height(X),1)]);

[fy_sweep_1100,fx_sweep_1100] = MF52_Combined(str2double(final_coefficients(:,2)),[Y,X,-1100.*ones(height(X),1),zeros(height(X),1)]);

[fy_sweep_900, fx_sweep_900] = MF52_Combined(str2double(final_coefficients(:,2)),[Y,X,-900.*ones(height(X),1),zeros(height(X),1)]);

[fy_sweep_700, fx_sweep_700] = MF52_Combined(str2double(final_coefficients(:,2)),[Y,X,-700.*ones(height(X),1),zeros(height(X),1)]);

[fy_sweep_500, fx_sweep_500] = MF52_Combined(str2double(final_coefficients(:,2)),[Y,X,-400.*ones(height(X),1),zeros(height(X),1)]);

figure
hold on
plot(fx_sweep_1300,fy_sweep_1300,'.','DisplayName','Fz = 1300')
plot(fx_sweep_1100,fy_sweep_1100,'.','DisplayName','Fz = 1100')
plot(fx_sweep_900,fy_sweep_900,'.','DisplayName','Fz = 900')
plot(fx_sweep_700,fy_sweep_700,'.','DisplayName','Fz = 700')
plot(fx_sweep_500,fy_sweep_500,'.','DisplayName','Fz = 400')
hold off
ylabel('FY (N)')
xlabel('FX (N)')
legend location northeast
title('Force Ellipse')
subtitle('R20 18"x6" 7" Rim @ 0 IA')