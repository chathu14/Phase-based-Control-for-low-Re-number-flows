clear all;    % Clear all variables from workspace
close all;    % Close all figure windows
clc;  
%% === Common Parameters ===
start_idx = 5001;
dt = 0.001;

%% Load Data
load('Pos_Stiffness/theta.mat');
load('Pos_Stiffness/time.mat');
load('Pos_Stiffness/moment.mat');
load('Pos_Stiffness/theta_dot.mat');

theta_data_1 = theta_data;
time_data_1 = time_data;
moment_data_1 = moment_data;
theta_dot_data_1 = theta_dot_data;

[~, num_columns_1] = size(time_data_1);
energy_result_1 = zeros(1, num_columns_1);

for col_1 = 1:num_columns_1
    theta_col_1     = theta_data_1(start_idx:end, col_1);
    moment_col_1    = moment_data_1(start_idx:end, col_1);
    theta_dot_col_1 = theta_dot_data_1(start_idx:end, col_1);
    time_col_1      = time_data_1(start_idx:end, 1) * dt;

    integral_val_1 = 0;
    for j_1 = 1:(length(time_col_1) - 1)
        dt_local_1 = time_col_1(j_1+1) - time_col_1(j_1);
        f1_1 = moment_col_1(j_1) * theta_dot_col_1(j_1);
        f2_1 = moment_col_1(j_1+1) * theta_dot_col_1(j_1+1);
        integral_val_1 = integral_val_1 + 0.5 * (f1_1 + f2_1) * dt_local_1;
    end

    energy_result_1(col_1) = integral_val_1;
end

%% Load Data (Neg stiff)
load('Neg_Stiffness/time.mat');       % loads variable 'time'
load('Neg_Aero_Stiffness/moment.mat');     % loads variable 'moment'
load('Neg_Aero_Stiffness/theta_dot.mat');  % loads variable 'theta_dot'
 start_idx = 5001;
 dt = 0.001;

time_data_2 = time;
moment_data_2 = moment;
theta_dot_data_2 = theta_dot;

[~, num_columns_2] = size(time_data_2);
energy_result_2 = zeros(1, num_columns_2);

for col_2 = 1:num_columns_2
    % No theta_data for part 2, so skip theta_col_2

    moment_col_2    = moment_data_2(start_idx:end, col_2);
    theta_dot_col_2 = theta_dot_data_2(start_idx:end, col_2);
    time_col_2      = time_data_2(start_idx:end, 1) * dt;

    integral_val_2 = 0;
    for j_2 = 1:(length(time_col_2) - 1)
        dt_local_2 = time_col_2(j_2+1) - time_col_2(j_2);
        f1_2 = moment_col_2(j_2) * theta_dot_col_2(j_2);
        f2_2 = moment_col_2(j_2+1) * theta_dot_col_2(j_2+1);
        integral_val_2 = integral_val_2 + 0.5 * (f1_2 + f2_2) * dt_local_2;
    end

    energy_result_2(col_2) = integral_val_2;
end

%% Load Phase and Plot 
phase = load('pert_phase.dat');

energy_plot_1 = energy_result_1;
energy_plot_2 = energy_result_2(2:end-1);
phase_plot_2  = phase(2:end-1);

figure;
hold on;
p1 = plot(phase, energy_plot_2, 'o', ...
    'Color', [0 0 1], ...    'MarkerFaceColor', 'none', ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'MarkerSize', 7);  % empty gray circle

p2 = plot(phase, energy_plot_1, '^', ...
    'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'MarkerSize', 7);  % filled black triangl

yline(0, 'k--', 'LineWidth', 1);  % Dashed black line at y = 0


% Labels and formatting
xlabel('\phi'); ylabel('E');

% xlabel('$\phi$', 'Interpreter', 'latex');
ylabel('\itE');
set(gca, 'FontSize', 18);

% Legend with same line+marker styles
% Add legend with labels
%legend([p1, p2], {'+ve perturbation', '-ve perturbation'}, 'Location', 'west');
legend([p1, p2], {'-\beta', '+\beta'}, 'Location', 'west', 'FontSize', 20);legend boxoff;
% Set x and y limits
xlim([0 2*pi]);

% ylimVals = ylim;
% if ylimVals(1) > -0.1 || ylimVals(2) < 0.2
%     ylim([-0.1 0.2]);
% else
%     ylim(ylimVals);
% end

% Set x-ticks and labels
xticks(0:pi/2:2*pi)
xticklabels({'0', 'π/2', 'π', '3π/2', '2π'})


% Major y-ticks at 0.1 interval
xlim([0 2*pi]);
ylim([-0.2 0.6]);

% Axis appearance
set(gca, 'TickLength', [0.02 0.02], 'LineWidth', 1.5)
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')

% Minor ticks
ax = gca;
ax.XAxis.MinorTickValues = 0:pi/6:2*pi;
xline(1.269, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
xline(4.537, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
xline(4.719, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
xline(4.719, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
xline(4.900, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
xline(5.263, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
xline(1.269, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
xline(5.445, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
xline(5.626, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');% yt = yticks;
% minorY = [];
% for i = 1:length(yt)-1
%     delta = (yt(i+1) - yt(i)) / 3;
%     minorY = [minorY, yt(i) + delta, yt(i) + 2*delta];
% end


box on;
grid off;

% Set figure size
%width = 4.955;   % inches
width = 6.7;   % inches

height = 4;      % inches
fig = gcf;
set(fig, 'Units', 'inches');
set(fig, 'Position', [1 1 width height]);
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [width height]);
set(fig, 'PaperPosition', [0 0 width height]);

%drawnow;
print('-depsc', 'energy_fig2.eps');
hold off;

