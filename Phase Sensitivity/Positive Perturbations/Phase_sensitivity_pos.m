clear;clc;close;

%% Load baseline data

dt = 0.001;
load('Baseline/time.mat');
time= time*dt-time(1)*dt;
load('Baseline/theta.mat');
alpha_org = theta;clear theta;
load('Baseline/moment.mat');
moment_org = moment; clear moment;
load('Baseline/CL.mat');
CL_org = CL; clear CL;
load('Baseline/CD.mat');
CD_org = CD; clear CD;
load('Baseline/theta_dot.mat');
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
    u_stiff(:,i) = input_impulse(time,-10,KickTime(i));
    u_heave(:,i) = input_impulse(time,0.01,KickTime(i));
    u_moment(:,i) = input_impulse(time,0.05,KickTime(i));
    u_pitch(:,i) = input_impulse(time,0.0001,KickTime(i));
end

%% Determine phase of oscillation
phase0    = time(1);%time(193);
phase     = 2*pi*(KickTime(2:end-1)-phase0)/T_org;

%% Load perturbation data

load Stiffness/theta.mat;KM_stiff = 10;
alpha_stiff_pert = theta_data; clear theta_data;
load Stiffness/moment.mat;
moment_stiff_pert = moment_data; clear moment_data;

Phase_Sensitivity_Stiff = calculatePhaseSensitivity(time, 94, 99.5,...
    KM_stiff, alpha_org, alpha_stiff_pert);

[~,ind2_stiff]  = find(max(alpha_stiff_pert)>0.1);
Phase_Sensitivity_Stiff_temp = Phase_Sensitivity_Stiff;Phase_Sensitivity_Stiff_temp(ind2_stiff) = nan;
Phase_Sensitivity_Stiff1 = [Phase_Sensitivity_Stiff(end-1); Phase_Sensitivity_Stiff(end); Phase_Sensitivity_Stiff_temp; Phase_Sensitivity_Stiff(1); Phase_Sensitivity_Stiff(2)];
phase1 = [phase(end-1)-2*pi phase(end)-2*pi phase 2*pi+phase(1) 2*pi+phase(2)];

phase_refined = -0.2:0.001:2*pi+0.2;
PS_Stiff = spline(phase1,Phase_Sensitivity_Stiff1,phase_refined);
PS_Grad_Stiff = gradient(PS_Stiff,phase_refined);

save('PRC_aeroelastic_Stiff.mat','phase_refined','PS_Stiff','PS_Grad_Stiff');

colormap_values = load('colormap_values.txt');

if size(colormap_values, 2) == 4
    colormap_values = colormap_values(:, 1:3); % Remove alpha channel if present
end

phase_normalized = (phase - min(phase)) / (max(phase) - min(phase));
n_bins = size(colormap_values, 1); % Number of colors in the colormap

color_indices = round(phase_normalized * (n_bins - 1)) + 1;

figure;
subplot(221);
hold on;

for i = 1:length(phase)
    scatter(phase(i), Phase_Sensitivity_Stiff_temp(i), 50, colormap_values(color_indices(i), :), 'filled');
    disp(phase(i));

end

for i = 1:length(ind2_stiff)
    scatter(phase(ind2_stiff(i)), Phase_Sensitivity_Stiff(ind2_stiff(i)), 80, ...
        colormap_values(color_indices(ind2_stiff(i)), :), 'filled'); 
end

plot(phase_refined, PS_Stiff, '-', 'Color', [0.5, 0.5, 0.5],'LineWidth', 1.5); % Line with gray color

% Vertical lines
xline(phase(ind2_stiff), 'k--');

% Labels and formatting
xlabel('\phi'); ylabel('Z(\phi)');
set(gca, 'FontSize', 18);
xlim([0 2*pi]);
ylim([-0.1 0.32]);
grid off;
box on

% Set x-ticks at multiples of pi
xticks(0:pi/2:2*pi)
xticklabels({'0', 'π/2', 'π', '3π/2', '2π'})
xtickangle(0)

% Axis appearance
set(gca, 'TickLength', [0.02 0.02], 'LineWidth', 1.5)
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')

ax = gca;
ax.XAxis.MinorTickValues = 0:pi/6:2*pi;

yticks(-0.1:0.1:0.4);

yt = yticks;
minorY = [];
for i = 1:length(yt)-1
    delta = (yt(i+1) - yt(i)) / 3;
    minorY = [minorY, yt(i) + delta, yt(i) + 2*delta];
end

ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = minorY;
grid off;
yline(0, '-', 'Color', [0 0 0], 'LineWidth', 1.0);
colormap(colormap_values);

width = 10.5;   % inches
 height = 8.7;  % inches
fig = gcf;
set(fig, 'Units', 'inches');
set(fig, 'Position', [1 1 width height]);  % [left bottom width height]
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [width height]);
set(fig, 'PaperPosition', [0 0 width height]);  % exact match
print('-depsc', 'pos_sensitivity.eps');


growth_rate = zeros(length(phase),1);
for  i = 1:length(phase)
    growth_rate(i) = growth_rate_func(time(1:end),moment_stiff_pert(1:end,i));
end

figure;
hold on;
for i = 1:length(phase)
    %scatter(phase(i), growth_rate(i), 50, colormap_values(color_indices(i), :), 'filled');
   % scatter(phase(i), growth_rate(i), 50, [0 0 0], 'filled');
    scatter(phase(i), growth_rate(i), 50, '^', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2);

    % Save phase values to a file (one value per row)
   writematrix(phase(:), 'pert_phase.dat', 'Delimiter', ' ');

end
% Vertical lines
xline(phase(ind2_stiff), 'k--');

% Labels and formatting
xlabel('\phi'); ylabel('\lambda');
%title('Stiffness Perturbation');
set(gca, 'FontSize', 18);
box on

xlim([0 2*pi]);
ylim([-0.005 0.0155]);

% Set x-ticks at multiples of pi
xticks(0:pi/2:2*pi)
xticklabels({'0', 'π/2', 'π', '3π/2', '2π'})

% Set major y-ticks at 0.1 interval
ylimVals = ylim;  % get current y-axis limits
yticks(ylimVals(1):0.1:ylimVals(2));

% Axis appearance
set(gca, 'TickLength', [0.02 0.02], 'LineWidth', 1.5)
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca,'Ytick',[-0.005:0.005:0.015])
% Add two minor ticks between each major x-tick (pi/6 spacing)
ax = gca;
ax.XAxis.MinorTickValues = 0:pi/6:2*pi;

% Add two minor ticks between each major y-tick
yline(0, 'k--');
% Optional: turn off grid
grid off;
% Save the figure
width = 4.75;   % inches
height = 4.1;  % inches
fig = gcf;
set(fig, 'Units', 'inches');
set(fig, 'Position', [1 1 width height]);  % [left bottom width height]
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [width height]);
set(fig, 'PaperPosition', [0 0 width height]);  % exact match
print('-depsc', 'growth_rate_pos.eps');

function Phase_Sensitivity = calculatePhaseSensitivity(Time, Tstart, Tend,...
    KickMag, data_org, data_pert)

% Define finer time intervals for precise interpolation
TI           = Tstart:0.00005:Tend;

% Interpolate the base response and find peak time
data_I       = interp1(Time, data_org, TI, 'spline');
[~, Index]   = max(-data_I);
PT_Base      = TI(Index);

% Preallocate for kick peak times
PT_Kick = nan(size(data_pert, 2), 1);

% Compute peak times for perturbed data
for i = 1:size(data_pert, 2)
    data_I       = interp1(Time, data_pert(:, i), TI, 'spline');
    %[~, Index]   = findpeaks(data_I);
    [~, Index] = max(-data_I);
    PT_Kick(i) = TI(Index);%TI(Index(13));%
end

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

