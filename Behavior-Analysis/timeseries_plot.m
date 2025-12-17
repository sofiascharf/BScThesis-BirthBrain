%% --- USER SETTINGS ---
refTime = datetime('21-Oct-2025 17:56:55','InputFormat','dd-MMM-yyyy HH:mm:ss'); % defines time = 0
startMinutes = -180;   % same window logic
endMinutes   = 60;     % optional upper limit (adjust as needed)

%% --- LOAD DATA ---
dataTS = readtable('all_calls_timeseries.csv');

% Expected columns: "EventDateTime" and "Frequency_kHz"
dataTS.EventDateTime = datetime(dataTS.EventDateTime, 'InputFormat','dd-MMM-yyyy HH:mm:ss');

% Compute minutes relative to refTime
dataTS.Minutes = minutes(dataTS.EventDateTime - refTime);

% Filter by chosen window
dataTS = dataTS(dataTS.Minutes >= startMinutes & dataTS.Minutes <= endMinutes, :);

%% --- PLOT SCATTER ---
figure('Color','w');
s = scatter(dataTS.Minutes, dataTS.Frequency_kHz, 30, 'black', 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',0.8);
s.MarkerFaceAlpha = 0.5;

% Add shading 
% Parturition
x_points_par = [0,0,15,15];
y_points = [0,12,12,0];
color_par = [0.75,0,0.75];
% Pre-parturition 
x_points_pre = [-60,-60,0,0];
color_pre = [0.98, 0.31, 0.73];
% Post-parturition
x_ponts_post = [15,15,60,60];
color_post = [0.494, 0.1840, 0.5560];

hold on;

area_par = fill(x_points_par, y_points, color_par);
area_par.FaceAlpha = 0.5;

area_pre = fill(x_points_pre, y_points, color_pre);
area_pre.FaceAlpha = 0.5;

area_post = fill(x_ponts_post, y_points, color_post);
area_post.FaceAlpha = 0.5;

xlabel('Time relative to Parturition (min)');
ylabel('Mean frequency (kHz)');
title('Call Frequency Over Time');
grid on;

% Prettify
set(gca,'Box','off','Color','none','LineWidth',1,'FontSize',11);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.GridColor = [0.85 0.85 0.85];

% Optional: set limits consistent with bar plots
xlim([startMinutes, endMinutes]);

%% --- LEGEND FOR SHADING REGIONS ---
% Create invisible patch handles for legend
patch_pre  = patch(nan, nan, color_pre,  'EdgeColor','none', 'FaceAlpha',0.5);
patch_par  = patch(nan, nan, color_par,  'EdgeColor','none', 'FaceAlpha',0.5);
patch_post = patch(nan, nan, color_post, 'EdgeColor','none', 'FaceAlpha',0.5);

legend([patch_pre, patch_par, patch_post], ...
       {'Pre-parturition', 'Parturition', 'Post-parturition'}, ...
       'Location','northwest', 'Box','off', 'FontSize',10);
