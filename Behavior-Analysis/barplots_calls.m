%% --- USER SETTINGS ---
refTime      = datetime('21-Oct-2025 17:56:55','InputFormat','dd-MMM-yyyy HH:mm:ss'); % time = 0
startMinutes = -120;   % lower limit for x-axis
binWidth10   = 10;     % bin width for 10-min file
binWidth20   = 20;     % bin width for 20-min file

% Phase boundaries (in minutes)
pre_parturition_start = -60;
parturition_start = 0;
parturition_end   = 15;

% Define colors (matching your scatter shades)
color_baseline = [1, 1, 1];
color_pre  = [0.98, 0.31, 0.73];        % pre-parturition (light purple)
color_par  = [0.75, 0, 0.75];           % parturition (same hue)
color_post = [0.494, 0.184, 0.556];     % post-parturition (darker purple)

% --- OUTPUT FOLDER ---
outdir = fullfile('.', 'Plots');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% --- LOAD DATA ---
data10 = readtable(fullfile('..','Data','Detections_Audible','calls_10min_summary.csv'));
data20 = readtable(fullfile('..','Data','Detections_Audible','calls_20min_summary.csv'));

% Convert time columns
data10.BinCenter = datetime(data10.BinCenter,'InputFormat','dd-MMM-yyyy HH:mm:ss');
data20.BinCenter = datetime(data20.BinCenter,'InputFormat','dd-MMM-yyyy HH:mm:ss');

% Compute bin start times relative to refTime
data10.StartMin = minutes(data10.BinCenter - refTime) - binWidth10/2;
data20.StartMin = minutes(data20.BinCenter - refTime) - binWidth20/2;

% Filter window
data10 = data10(data10.StartMin >= startMinutes, :);
data20 = data20(data20.StartMin >= startMinutes, :);

% Ensure a bin starting at time 0 exists
if ~ismember(0, data10.StartMin)
    newRow = array2table(zeros(1,width(data10)),'VariableNames',data10.Properties.VariableNames);
    newRow.BinCenter = refTime + minutes(binWidth10/2);
    newRow.CallCount = 0;
    newRow.StartMin  = 0;
    data10 = [data10; newRow];
end
if ~ismember(0, data20.StartMin)
    newRow = array2table(zeros(1,width(data20)),'VariableNames',data20.Properties.VariableNames);
    newRow.BinCenter = refTime + minutes(binWidth20/2);
    newRow.CallCount = 0;
    newRow.StartMin  = 0;
    data20 = [data20; newRow];
end

% Sort after insertion
data10 = sortrows(data10,'StartMin');
data20 = sortrows(data20,'StartMin');

%% --- ASSIGN COLORS + PHASE LABELS + RATES (CALLS/MIN) ---
barColors10   = zeros(height(data10),3);
phaseLabels10 = strings(height(data10),1);

for i = 1:height(data10)
    t = data10.StartMin(i);
    if t < pre_parturition_start
        barColors10(i, :) = color_baseline;
        phaseLabels10(i)  = "baseline";
    elseif t < parturition_start && t >= pre_parturition_start
        barColors10(i,:)  = color_pre;
        phaseLabels10(i)  = "pre";
    elseif t >= parturition_start && t < parturition_end
        barColors10(i,:)  = color_par;
        phaseLabels10(i)  = "par";
    else
        barColors10(i,:)  = color_post;
        phaseLabels10(i)  = "post";
    end
end

barColors20   = zeros(height(data20),3);
phaseLabels20 = strings(height(data20),1);

for i = 1:height(data20)
    t = data20.StartMin(i);
    if t < pre_parturition_start
        barColors20(i, :) = color_baseline;
        phaseLabels20(i)  = "baseline";
    elseif t < parturition_start && t >= pre_parturition_start
        barColors20(i,:)  = color_pre;
        phaseLabels20(i)  = "pre";
    elseif t >= parturition_start && t < parturition_end
        barColors20(i,:)  = color_par;
        phaseLabels20(i)  = "par";
    else
        barColors20(i,:)  = color_post;
        phaseLabels20(i)  = "post";
    end
end

data10.Phase = categorical(phaseLabels10, {'baseline','pre','par','post'});
data20.Phase = categorical(phaseLabels20, {'baseline','pre','par','post'});

data10.Rate = data10.CallCount ./ binWidth10;  % calls per minute
data20.Rate = data20.CallCount ./ binWidth20;  % calls per minute

%% --- PLOT 10-MIN BINS ---
figure('Color','w');
barWidth = 0.8;
b = bar(data10.StartMin, data10.CallCount, barWidth, ...
        'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b.CData = barColors10;

xticks(data10.StartMin);
xlabels10 = strings(height(data10),1);
for i = 1:height(data10)
    xlabels10(i) = sprintf('%d to %d', round(data10.StartMin(i)), ...
                           round(data10.StartMin(i)+binWidth10));
end
set(gca,'XTickLabel',xlabels10,'XTickLabelRotation',45);

xlabel('Time window relative to Parturition (min)');
ylabel('Number of calls');
title('10-min Binned Calls');
grid on;
set(gca,'Box','off','Color','none','LineWidth',1,'FontSize',11);
ax = gca; ax.GridColor = [0.8 0.8 0.8];
xlim([startMinutes, max(data10.StartMin)+binWidth10]);

hold on;
bar_pre  = bar(nan, nan, 'FaceColor', color_pre,  'EdgeColor','k', 'LineWidth',1.2);
bar_par  = bar(nan, nan, 'FaceColor', color_par,  'EdgeColor','k', 'LineWidth',1.2);
bar_post = bar(nan, nan, 'FaceColor', color_post, 'EdgeColor','k', 'LineWidth',1.2);
legend([bar_pre, bar_par, bar_post], ...
       {'Pre-parturition', 'Parturition', 'Post-parturition'}, ...
       'Location', 'northwest', 'Box', 'off', 'FontSize', 10);

fig10 = gcf;
save_if_missing(fig10, fullfile(outdir, 'AudSqueaks10'));

%% --- PLOT 20-MIN BINS ---
figure('Color','w');
b2 = bar(data20.StartMin, data20.CallCount, barWidth, ...
         'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b2.CData = barColors20;

xticks(data20.StartMin);
xlabels20 = strings(height(data20),1);
for i = 1:height(data20)
    xlabels20(i) = sprintf('%d to %d', round(data20.StartMin(i)), ...
                           round(data20.StartMin(i)+binWidth20));
end
set(gca,'XTickLabel',xlabels20,'XTickLabelRotation',45);

xlabel('Time window relative to Parturition (min)');
ylabel('Number of calls');
title('20-min Binned Calls');
grid on;
set(gca,'Box','off','Color','none','LineWidth',1,'FontSize',11);
ax = gca; ax.GridColor = [0.8 0.8 0.8];
xlim([startMinutes, max(data20.StartMin)+binWidth20]);

hold on;
bar_pre  = bar(nan, nan, 'FaceColor', color_pre,  'EdgeColor','k', 'LineWidth',1.2);
bar_par  = bar(nan, nan, 'FaceColor', color_par,  'EdgeColor','k', 'LineWidth',1.2);
bar_post = bar(nan, nan, 'FaceColor', color_post, 'EdgeColor','k', 'LineWidth',1.2);
legend([bar_pre, bar_par, bar_post], ...
       {'Pre-parturition', 'Parturition', 'Post-parturition'}, ...
       'Location', 'northwest', 'Box', 'off', 'FontSize', 10);

fig20 = gcf;
save_if_missing(fig20, fullfile(outdir, 'AudSqueaks20'));

%% --- STATS ON 10-MIN RATES (CALLS/MIN) ACROSS PHASES ---
y = data10.Rate;
g = data10.Phase;
phaseCats = categories(g);  % {'baseline','pre','par','post'}

% One-way ANOVA
[p_anova10, tbl_anova10, stats_anova10] = anova1(y, g, 'off');

% Residuals (group-mean predicted)
yhat = zeros(size(y));
for i = 1:numel(phaseCats)
    idx = g == phaseCats{i};
    yhat(idx) = mean(y(idx), 'omitnan');
end
resid10 = y - yhat;
resid10 = resid10(~isnan(resid10));

% Normality test
if numel(resid10) >= 5
    [h_norm10, p_norm10] = lillietest(resid10);
else
    h_norm10 = 1;
    p_norm10 = NaN;
end

use_nonparam = (h_norm10 ~= 0);

if ~use_nonparam
    fprintf('\nResiduals normal-ish (Lilliefors p = %.3g). Using ANOVA.\n', p_norm10);
    fprintf('ANOVA p = %.3g\n', p_anova10);
    c = multcompare(stats_anova10, 'CType','tukey', 'Display','off');
    gnames = cellstr(stats_anova10.gnames);
else
    fprintf('\nResiduals non-normal (Lilliefors p = %.3g). Using Kruskal–Wallis.\n', p_norm10);
    [p_kw10, tbl_kw10, stats_kw10] = kruskalwallis(y, g, 'off');
    fprintf('KW p = %.3g\n', p_kw10);
    c = multcompare(stats_kw10, 'Display','off');
    gnames = cellstr(stats_kw10.gnames);
end

% Build pairwise p-value matrix IN PHASE CATEGORY ORDER
nPhase = numel(phaseCats);
pairP  = nan(nPhase);
for k = 1:size(c,1)
    % indices in stats order
    i1_stats = c(k,1);
    i2_stats = c(k,2);
    pval     = c(k,6);

    name1 = strtrim(gnames{i1_stats});
    name2 = strtrim(gnames{i2_stats});

    % indices in phaseCats order (plot order)
    i1_cat = find(strcmp(phaseCats, name1));
    i2_cat = find(strcmp(phaseCats, name2));

    if ~isempty(i1_cat) && ~isempty(i2_cat)
        pairP(i1_cat, i2_cat) = pval;
        pairP(i2_cat, i1_cat) = pval;
    end
end

%% --- PHASE MEAN BARPLOT + STARS FOR SIGNIFICANT PAIRS ---
meanRate = zeros(nPhase,1);
semRate  = zeros(nPhase,1);

for i = 1:nPhase
    idx = g == phaseCats{i};
    r   = y(idx);
    meanRate(i) = mean(r, 'omitnan');
    semRate(i)  = std(r, 'omitnan') / sqrt(sum(~isnan(r)));
end

figure('Color','w');
bh = bar(1:nPhase, meanRate, 'FaceColor','flat', ...
    'EdgeColor','k', 'LineWidth',1.2);

% Color bars by phase
cPhase = zeros(nPhase,3);
for i = 1:nPhase
    switch char(phaseCats{i})
        case 'baseline'
            cPhase(i,:) = color_baseline;
        case 'pre'
            cPhase(i,:) = color_pre;
        case 'par'
            cPhase(i,:) = color_par;
        case 'post'
            cPhase(i,:) = color_post;
        otherwise
            cPhase(i,:) = [0.7 0.7 0.7];
    end
end
bh.CData = cPhase;

hold on;
errorbar(1:nPhase, meanRate, semRate, semRate, ...
    'k', 'LineStyle','none', 'LineWidth',1.2);

set(gca,'XTick',1:nPhase,'XTickLabel',phaseCats,'FontSize',11);
xlabel('Phase');
ylabel('Calls / min');
title('Call rate per phase (10-min bins)');
box off;
grid on;
ax = gca; ax.GridColor = [0.8 0.8 0.8];

% --- Pairwise significance bars: ONLY significant pairs get stars ---
alpha_sig = 0.05;

yBase = max(meanRate + semRate);
yStep = range(meanRate + semRate);
if yStep == 0
    yStep = max(meanRate + semRate) * 0.2 + 1;
end
yCurrent = yBase + 0.15*yStep;

pairs  = nchoosek(1:nPhase,2);
pairPs = arrayfun(@(k) pairP(pairs(k,1), pairs(k,2)), 1:size(pairs,1));
[~, order] = sort(pairPs, 'ascend', 'MissingPlacement','last');
pairs  = pairs(order,:);
pairPs = pairPs(order);

for kk = 1:size(pairs,1)
    i1   = pairs(kk,1);
    i2   = pairs(kk,2);
    pval = pairPs(kk);

    if isnan(pval) || pval >= alpha_sig
        continue;   % only draw significant ones
    end

    x1 = i1; x2 = i2;
    yLine = yCurrent;
    tickHeight = 0.015*yStep;

    plot([x1 x2], [yLine yLine], 'k-', 'LineWidth',1.3);
    plot([x1 x1], [yLine-tickHeight yLine], 'k-', 'LineWidth',1.3);
    plot([x2 x2], [yLine-tickHeight yLine], 'k-', 'LineWidth',1.3);

    text((x1+x2)/2, yLine + 0.01*yStep, p2stars(pval), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'Color','k', ...
        'FontSize',12, 'FontWeight','bold');

    yCurrent = yCurrent + 0.12*yStep;
end

figPhase = gcf;
save_if_missing(figPhase, fullfile(outdir, 'AudSqueaks_PhaseComparison'));

%% --- HELPER FUNCTIONS ---------------------------------------------------
function save_if_missing(figHandle, basepath)
    if ~endsWith(basepath, '.fig', 'IgnoreCase', true)
        basepath = basepath + ".fig";
    end
    if ~isfile(basepath)
        savefig(figHandle, basepath);
        fprintf('Saved: %s\n', basepath);
    else
        fprintf('Exists, not overwriting: %s\n', basepath);
    end
end

function s = p2stars(p)
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = '';
    end
end
