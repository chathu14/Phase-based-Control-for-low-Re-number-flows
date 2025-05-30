clear;clc;close;

%% Load baseline data

dt = 0.001;
load('Data/Aero_Baseline/time.mat');
time= time*dt-time(1)*dt;
load('Data/Aero_Baseline/theta.mat');
alpha_org = theta;clear theta;
load('Data/Aero_Baseline/moment.mat');
moment_org = moment; clear moment;
load('Data/Aero_Baseline/CL.mat');
CL_org = CL; clear CL;
load('Data/Aero_Baseline/CD.mat');
CD_org = CD; clear CD;
load('Data/Aero_Baseline/theta_dot.mat');
alpha_dot_org = theta_dot; clear theta_dot;

%% Calculate FFT

[Freq,yamp]       = Func_sptrm(alpha_org,1/dt,2);
[Freq_CL,yamp_CL] = Func_sptrm(CL_org,1/dt,2);
[Freq_m,yamp_m]   = Func_sptrm(moment_org,1/dt,2);
[a,b]             = max(yamp);
T_org             = 1/Freq(b);
[a_CL,b_CL]       = max(yamp_CL(2:end));
T_org_CL          = 1/Freq_CL(b_CL+1);
[a_m,b_m]         = max(yamp_m(2:end));
T_org_m           = 1/Freq_m(b_m+1);

%% Pert simulations parameters

T_pert_period = 5.000000000023873;   % needed to be T_org above
Initial_Time  = 0.14983357145777632; % initial kicktime
N             = 32;                  %
KickTime      = linspace(Initial_Time, Initial_Time+T_pert_period, N+2) ;

%% Determine input impulse as a times-series

u_stiff = zeros(length(time),length(KickTime));
u_heave = zeros(length(time),length(KickTime));
u_moment = zeros(length(time),length(KickTime));
u_pitch = zeros(length(time),length(KickTime));
for i = 1:length(KickTime)
    u_stiff(:,i) = input_impulse(time,10,KickTime(i));
    u_heave(:,i) = input_impulse(time,0.01,KickTime(i));
    u_moment(:,i) = input_impulse(time,0.05,KickTime(i));
    u_pitch(:,i) = input_impulse(time,0.0001,KickTime(i));
end

%% Determine phase of oscillation
phase0    = time(1);%time(193);
phase     = 2*pi*(KickTime(2:end-1)-phase0)/T_org;

%% Load perturbation data

load Data/Aero_Stiffness/theta.mat;KM_stiff = -10;
alpha_stiff_pert = theta; clear theta;
load Data/Aero_Stiffness/moment.mat;
moment_stiff_pert = moment;clear moment;
load Data/Aero_Heave/theta.mat;KM_heave = -0.01;
alpha_heave_pert = theta;clear theta;
load Data/Aero_Heave/moment.mat;
moment_heave_pert = moment;clear moment;

Phase_Sensitivity_Stiff = calculatePhaseSensitivity(time, 94, 99.5,...
    KM_stiff, alpha_org, alpha_stiff_pert);

Phase_Sensitivity_Heave = calculatePhaseSensitivity(time, 95, 96.5,...
    KM_heave, alpha_org, alpha_heave_pert);

[~,ind2_stiff]  = find(max(alpha_stiff_pert)>0.1);
[~,ind2_heave]  = find(max(alpha_heave_pert)>0.1);
Phase_Sensitivity_Stiff_temp = Phase_Sensitivity_Stiff;Phase_Sensitivity_Stiff_temp(ind2_stiff) = nan;
Phase_Sensitivity_Heave_temp = Phase_Sensitivity_Heave;Phase_Sensitivity_Heave_temp(24:26) = nan;
Phase_Sensitivity_Heave1 = [Phase_Sensitivity_Heave(end-1); Phase_Sensitivity_Heave(end); Phase_Sensitivity_Heave_temp; Phase_Sensitivity_Heave(1); Phase_Sensitivity_Heave(2)];
Phase_Sensitivity_Stiff1 = [Phase_Sensitivity_Stiff(end-1); Phase_Sensitivity_Stiff(end); Phase_Sensitivity_Stiff_temp; Phase_Sensitivity_Stiff(1); Phase_Sensitivity_Stiff(2)];
phase1 = [phase(end-1)-2*pi phase(end)-2*pi phase 2*pi+phase(1) 2*pi+phase(2)];

phase_refined = -0.2:0.001:2*pi+0.2;
PS_Heave = spline(phase1,Phase_Sensitivity_Heave1,phase_refined);
PS_Stiff = spline(phase1,Phase_Sensitivity_Stiff1,phase_refined);
PS_Grad_Heave = gradient(PS_Heave,phase_refined);
PS_Grad_Stiff = gradient(PS_Stiff,phase_refined);

save('PRC_aeroelastic_Heave.mat','phase_refined','PS_Heave','PS_Grad_Heave');
save('PRC_aeroelastic_Stiff.mat','phase_refined','PS_Stiff','PS_Grad_Stiff');

% Load the custom colormap
colormap_values = load('colormap_values.txt');

% Ensure colormap_values has the correct dimensions
if size(colormap_values, 2) == 4
    colormap_values = colormap_values(:, 1:3); % Remove alpha channel if present
end

% Normalize phase data to the range [0, 1]
phase_normalized = (phase - min(phase)) / (max(phase) - min(phase));
n_bins = size(colormap_values, 1); % Number of colors in the colormap

% Map normalized phase to colormap indices
color_indices = round(phase_normalized * (n_bins - 1)) + 1;

figure;

% Stiffness Perturbation Plot
subplot(221);
hold on;

% Plot all points with scatter using the colormap
for i = 1:length(phase)
    scatter(phase(i), Phase_Sensitivity_Stiff_temp(i), 50, colormap_values(color_indices(i), :), 'filled');
end

% Highlight specific points (ind2_stiff) with the colormap
for i = 1:length(ind2_stiff)
    scatter(phase(ind2_stiff(i)), Phase_Sensitivity_Stiff(ind2_stiff(i)), 80, ...
        colormap_values(color_indices(ind2_stiff(i)), :), 'filled'); % No 's' for square markers
end

% Add the refined curve
plot(phase_refined, PS_Stiff, '-', 'Color', [0.8, 0, 0]); % Line with a reddish color

% Vertical lines
xline(phase(ind2_stiff), 'k--');

% Labels and formatting
xlabel('\phi'); ylabel('Z(\phi)');
title('Stiffness Perturbation');
set(gca, 'FontSize', 14);
xlim([0 2*pi]);
grid off;

% Heave Perturbation Plot
subplot(222);
hold on;
for i = 1:length(phase)
    scatter(phase(i), Phase_Sensitivity_Heave(i), 50, colormap_values(color_indices(i), :), 'filled');
end
plot(phase_refined, PS_Heave, '-', 'Color', [0, 0, 0.8]); % Line with a bluish color
plot(phase(ind2_heave), Phase_Sensitivity_Heave(ind2_heave), 'bs', 'MarkerFaceColor', 'b');
xline(phase(ind2_stiff), 'k--');
xlabel('\phi'); ylabel('Z(\phi)');
title('Heave Perturbation');
set(gca, 'FontSize', 14);
xlim([0 2*pi]);
grid off;

% Add a colormap to the figure
colormap(colormap_values);

% Remove colorbar by commenting out or deleting the following line:
% colorbar;

% Save the figure
print('-depsc', 'Sensitivity_1.eps');

figure;
subplot(211);
plot(time,alpha_stiff_pert(:,5),'-');hold on;
plot(time,alpha_stiff_pert(:,6),'-');hold on;
plot(time,alpha_org,'k-');hold on;
grid off;%grid minor;
xlabel('time');ylabel('\alpha');
title('Stiffness pert');
ylim([-0.4,0.4]);
set(gca,'Fontsize',14);
subplot(212);
plot(time,moment_stiff_pert(:,5),'-');hold on;
plot(time,moment_stiff_pert(:,6),'-');hold on;
plot(time,moment_org,'k-');hold on;
grid off; %grid minor;
xlabel('time');ylabel('C_M');
title('Stiffness pert');
ylim([-0.1,0.3]);
set(gca,'Fontsize',14);
print('-depsc', 'Response.eps');



function Phase_Sensitivity = calculatePhaseSensitivity(Time, Tstart, Tend,...
    KickMag, data_org, data_pert)

% Define finer time intervals for precise interpolation
TI           = Tstart:0.00005:Tend;

% Interpolate the base response and find peak time
data_I       = interp1(Time, data_org, TI, 'spline');
[~, Index]   = max(data_I);
PT_Base      = TI(Index);

% Preallocate for kick peak times
PT_Kick = nan(size(data_pert, 2), 1);

% Compute peak times for perturbed data
for i = 1:size(data_pert, 2)
    data_I       = interp1(Time, data_pert(:, i), TI, 'spline');
    %[~, Index]   = findpeaks(data_I);
    [~, Index] = max(data_I);
    PT_Kick(i) = TI(Index);%TI(Index(13));%
end

% Wrap around for cyclic data
%PT_Kick(N + 1) = PT_Kick(1);

% Calculate phase sensitivity
Phase_Sensitivity = (PT_Base - PT_Kick) / KickMag;
end


function [Freq, YSptrm] = Func_sptrm(st, fs, flag)
% Gives freq vector and amp vector bu FFT
L_Y    = length(st) ;
avg    = mean(st) ;
Y      = st - avg ;

if flag == 1
    if mod(L_Y, 2) ~= 0
        warning('Vector length must be even, last element is truncated.')
        Y(end) = [] ;
    end
    NFFT = L_Y ;

elseif flag == 2
    NFFT = 2^nextpow2(L_Y) ;
end

Freq   = transpose(fs/2*linspace(0, 1, NFFT/2+1)) ;
Yfft   = fft(Y, NFFT)/L_Y ;
YAmp   = Yfft(1:(0.5*NFFT+1)) ;
YSptrm = 2 *abs(YAmp) ;
YSptrm = transpose(YSptrm) ;
YSptrm(1) = abs(avg) ;

% disp(['FFT length: ' num2str(length(YAmp))])
% Phase = phase(Yshift((0.5*L_Y+1):end)) ;
% Phase = transpose(Phase) ;

end

function u = input_impulse(time,I,KickTime)

% Parameters
sigma = 0.01;    % Standard deviation
% Gaussian impulse function
u = (I / (sqrt(2*pi) * sigma)) * exp(-0.5 * ((time - KickTime) / sigma).^2);
end
