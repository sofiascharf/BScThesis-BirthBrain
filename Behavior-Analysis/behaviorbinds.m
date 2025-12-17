%% --- USER SETTINGS (match your call plots) ---
refTime      = datetime('21-Oct-2025 17:56:55','InputFormat','dd-MMM-yyyy HH:mm:ss'); % time = 0
startMinutes = -120;   % lower limit for x-axis
binWidth10   = 10;     % minutes (must match calls_10min_summary.csv)
binWidth20   = 20;     % minutes (must match calls_20min_summary.csv)

% Phase boundaries (in minutes)
baseline_start    = -120;   % forced white
baseline_end      = -60;
parturition_start = 0;
parturition_end   = 15;

% Colors (same palette; baseline must be RGB to work with CData)
color_baseline = [1, 1, 1];          % white
color_pre      = [0.98, 0.31, 0.73]; % pre-parturition (−60..0)
color_par      = [0.75, 0, 0.75];    % 0..15
color_post     = [0.494, 0.184, 0.556];  % >=15

% Output folder for summary plots (optional)
outdir = fullfile('.', 'Plots');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% --- LOAD EXISTING CALL BINS TO MIRROR WINDOWS EXACTLY ---
data10 = readtable(fullfile('..','Data','Detections_Audible','calls_10min_summary.csv'));
data20 = readtable(fullfile('..','Data','Detections_Audible','calls_20min_summary.csv'));

data10.BinCenter = datetime(data10.BinCenter,'InputFormat','dd-MMM-yyyy HH:mm:ss');
data20.BinCenter = datetime(data20.BinCenter,'InputFormat','dd-MMM-yyyy HH:mm:ss');

data10.StartMin = minutes(data10.BinCenter - refTime) - binWidth10/2;
data20.StartMin = minutes(data20.BinCenter - refTime) - binWidth20/2;

data10 = data10(data10.StartMin >= startMinutes, :);
data20 = data20(data20.StartMin >= startMinutes, :);

% Ensure a bin starting at t=0 exists (so colors/labels line up)
if ~ismember(0, data10.StartMin)
    newRow = array2table(zeros(1,width(data10)),'VariableNames',data10.Properties.VariableNames);
    newRow.BinCenter = refTime + minutes(binWidth10/2);
    newRow.StartMin  = 0;
    data10 = [data10; newRow];
end
if ~ismember(0, data20.StartMin)
    newRow = array2table(zeros(1,width(data20)),'VariableNames',data20.Properties.VariableNames);
    newRow.BinCenter = refTime + minutes(binWidth20/2);
    newRow.StartMin  = 0;
    data20 = [data20; newRow];
end

data10 = sortrows(data10,'StartMin');
data20 = sortrows(data20,'StartMin');

%% --- Assign bin colors with BASELINE override ----------------------------
makeBarColors = @(starts, w) arrayfun(@(x) ...
    pickColorWithBaseline(x, baseline_start, baseline_end, ...
                              parturition_start, parturition_end, ...
                              color_baseline, color_pre, color_par, color_post), ...
    starts, 'UniformOutput', false);
colorMat = @(Ccell) vertcat(Ccell{:});

barColors10 = colorMat(makeBarColors(data10.StartMin, binWidth10));
barColors20 = colorMat(makeBarColors(data20.StartMin, binWidth20));

% --- Assign Phase labels (baseline, pre, par, post) for stats ---
phaseLabels10 = strings(height(data10),1);
for i = 1:height(data10)
    t = data10.StartMin(i);
    if t >= baseline_start && t < baseline_end
        phaseLabels10(i) = "baseline";
    elseif t >= baseline_end && t < parturition_start
        phaseLabels10(i) = "pre";
    elseif t >= parturition_start && t < parturition_end
        phaseLabels10(i) = "par";
    else
        phaseLabels10(i) = "post";
    end
end
data10.Phase = categorical(phaseLabels10, {'baseline','pre','par','post'});

% (optional) also label 20-min bins if needed later
phaseLabels20 = strings(height(data20),1);
for i = 1:height(data20)
    t = data20.StartMin(i);
    if t >= baseline_start && t < baseline_end
        phaseLabels20(i) = "baseline";
    elseif t >= baseline_end && t < parturition_start
        phaseLabels20(i) = "pre";
    elseif t >= parturition_start && t < parturition_end
        phaseLabels20(i) = "par";
    else
        phaseLabels20(i) = "post";
    end
end
data20.Phase = categorical(phaseLabels20, {'baseline','pre','par','post'});

%% --- LOAD BEHAVIOR TIME SERIES ---
TS = readtable(fullfile('..','Data','Annotations','timeseries_annotations.csv'));

TS.Minutes = minutes(TS.EventDateTime - refTime);

% Keep only the window we care about (a bit beyond last bin end for safety)
lastEnd10 = max(data10.StartMin) + binWidth10;
lastEnd20 = max(data20.StartMin) + binWidth20;
lastEnd   = max(lastEnd10, lastEnd20);
TS = TS(TS.Minutes >= startMinutes & TS.Minutes < lastEnd, :);

% ACTIVE-anywhere: active_in_nest OR outside_nest
TS.ActiveAnywhere = TS.active_in_nest | TS.outside_nest;

%% --- Helper to compute % in-bin for a logical vector at 1 s resolution ---
pctInBin = @(mins, logicalVec, bStart, bWidth) ...
    mean( logicalVec( mins >= bStart & mins < (bStart + bWidth) ) ) * 100;

%% --- 1) % TIME IN NEST ---------------------------------------------------
pctInNest10 = arrayfun(@(s) pctInBin(TS.Minutes, TS.inside_nest, s, binWidth10), data10.StartMin);
pctInNest20 = arrayfun(@(s) pctInBin(TS.Minutes, TS.inside_nest, s, binWidth20), data20.StartMin);

% Plot 10-min time course
figure('Color','w');
barWidth = 0.8;
b = bar(data10.StartMin, pctInNest10, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b.CData = barColors10;
xticks(data10.StartMin);
xlabels10 = strings(height(data10),1);
for i = 1:height(data10)
    xlabels10(i) = sprintf('%d to %d', round(data10.StartMin(i)), round(data10.StartMin(i)+binWidth10));
end
set(gca,'XTickLabel',xlabels10,'XTickLabelRotation',45);
xlabel('Time window relative to Parturition (min)');
ylabel('% of time in bin');
title('% Time IN Nest (10-min bins)');
grid on; box off; set(gca,'Color','none','LineWidth',1,'FontSize',11);
ax=gca; ax.GridColor=[0.8 0.8 0.8];
xlim([startMinutes, max(data10.StartMin)+binWidth10]);
ylim([0 100]);

% Plot 20-min time course
figure('Color','w');
b2 = bar(data20.StartMin, pctInNest20, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b2.CData = barColors20;
xticks(data20.StartMin);
xlabels20 = strings(height(data20),1);
for i = 1:height(data20)
    xlabels20(i) = sprintf('%d to %d', round(data20.StartMin(i)), round(data20.StartMin(i)+binWidth20));
end
set(gca,'XTickLabel',xlabels20,'XTickLabelRotation',45);
xlabel('Time window relative to Parturition (min)');
ylabel('% of time in bin');
title('% Time IN Nest (20-min bins)');
grid on; box off; set(gca,'Color','none','LineWidth',1,'FontSize',11);
ax=gca; ax.GridColor=[0.8 0.8 0.8];
xlim([startMinutes, max(data20.StartMin)+binWidth20]);
ylim([0 100]);

% Phase summary + stars (10-min)
phaseSummaryPlotWithStats(pctInNest10(:), data10.Phase, ...
    color_baseline, color_pre, color_par, color_post, ...
    '% Time IN Nest (10-min, per phase)', ...
    outdir, 'PhaseSummary_InNest');

%% --- 2) % TIME ACTIVE (active_in_nest OR outside_nest) -------------------
pctActive10 = arrayfun(@(s) pctInBin(TS.Minutes, TS.ActiveAnywhere, s, binWidth10), data10.StartMin);
pctActive20 = arrayfun(@(s) pctInBin(TS.Minutes, TS.ActiveAnywhere, s, binWidth20), data20.StartMin);

% Plot 10-min
figure('Color','w');
b = bar(data10.StartMin, pctActive10, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b.CData = barColors10;
xticks(data10.StartMin);
set(gca,'XTickLabel',xlabels10,'XTickLabelRotation',45);
xlabel('Time window relative to Parturition (min)');
ylabel('% of time in bin');
title('% Time ACTIVE (10-min bins)');
grid on; box off; set(gca,'Color','none','LineWidth',1,'FontSize',11);
ax=gca; ax.GridColor=[0.8 0.8 0.8];
xlim([startMinutes, max(data10.StartMin)+binWidth10]);
ylim([0 100]);

% Plot 20-min
figure('Color','w');
b2 = bar(data20.StartMin, pctActive20, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b2.CData = barColors20;
xticks(data20.StartMin);
set(gca,'XTickLabel',xlabels20,'XTickLabelRotation',45);
xlabel('Time window relative to Parturition (min)');
ylabel('% of time in bin');
title('% Time ACTIVE (20-min bins)');
grid on; box off; set(gca,'Color','none','LineWidth',1,'FontSize',11);
ax=gca; ax.GridColor=[0.8 0.8 0.8];
xlim([startMinutes, max(data20.StartMin)+binWidth20]);
ylim([0 100]);

% Phase summary + stars
phaseSummaryPlotWithStats(pctActive10(:), data10.Phase, ...
    color_baseline, color_pre, color_par, color_post, ...
    '% Time ACTIVE (10-min, per phase)', ...
    outdir, 'PhaseSummary_Active');

%% --- 3) % TIME CIRCLING (anywhere) --------------------------------------
vars = string(TS.Properties.VariableNames);
circCols = contains(lower(vars), 'circling');
if any(circCols)
    TS.CirclingAnywhere = false(height(TS),1);
    for v = vars(circCols)
        TS.CirclingAnywhere = TS.CirclingAnywhere | TS.(v);
    end
else
    warning('No columns containing "circling" found in timeseries_annotations.csv. Using zeros.');
    TS.CirclingAnywhere = false(height(TS),1);
end

pctCircling10 = arrayfun(@(s) pctInBin(TS.Minutes, TS.CirclingAnywhere, s, binWidth10), data10.StartMin);
pctCircling20 = arrayfun(@(s) pctInBin(TS.Minutes, TS.CirclingAnywhere, s, binWidth20), data20.StartMin);

% Plot 10-min
figure('Color','w');
b = bar(data10.StartMin, pctCircling10, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b.CData = barColors10;
xticks(data10.StartMin);
set(gca,'XTickLabel',xlabels10,'XTickLabelRotation',45);
xlabel('Time window relative to Parturition (min)');
ylabel('% of time in bin');
title('% Time CIRCLING (10-min bins)');
grid on; box off; set(gca,'Color','none','LineWidth',1,'FontSize',11);
xlim([startMinutes, max(data10.StartMin)+binWidth10]); ylim([0 100]);

% Plot 20-min
figure('Color','w');
b2 = bar(data20.StartMin, pctCircling20, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b2.CData = barColors20;
xticks(data20.StartMin);
set(gca,'XTickLabel',xlabels20,'XTickLabelRotation',45);
xlabel('Time window relative to Parturition (min)');
ylabel('% of time in bin');
title('% Time CIRCLING (20-min bins)');
grid on; box off; set(gca,'Color','none','LineWidth',1,'FontSize',11);
xlim([startMinutes, max(data20.StartMin)+binWidth20]); ylim([0 100]);

% Phase summary + stars
phaseSummaryPlotWithStats(pctCircling10(:), data10.Phase, ...
    color_baseline, color_pre, color_par, color_post, ...
    '% Time CIRCLING (10-min, per phase)', ...
    outdir, 'PhaseSummary_Circling');

%% --- 4) % TIME NESTING (anywhere) ---------------------------------------
nestCols = contains(lower(vars), 'nesting');
if any(nestCols)
    TS.NestingAnywhere = false(height(TS),1);
    for v = vars(nestCols)
        TS.NestingAnywhere = TS.NestingAnywhere | TS.(v);
    end
else
    warning('No columns containing "nesting" found in timeseries_annotations.csv. Using zeros.');
    TS.NestingAnywhere = false(height(TS),1);
end

pctNesting10 = arrayfun(@(s) pctInBin(TS.Minutes, TS.NestingAnywhere, s, binWidth10), data10.StartMin);
pctNesting20 = arrayfun(@(s) pctInBin(TS.Minutes, TS.NestingAnywhere, s, binWidth20), data20.StartMin);

% Plot 10-min
figure('Color','w');
b = bar(data10.StartMin, pctNesting10, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b.CData = barColors10;
xticks(data10.StartMin);
set(gca,'XTickLabel',xlabels10,'XTickLabelRotation',45);
xlabel('Time window relative to Parturition (min)');
ylabel('% of time in bin');
title('% Time NESTING (10-min bins)');
grid on; box off; set(gca,'Color','none','LineWidth',1,'FontSize',11);
xlim([startMinutes, max(data10.StartMin)+binWidth10]); ylim([0 100]);

% Plot 20-min
figure('Color','w');
b2 = bar(data20.StartMin, pctNesting20, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
b2.CData = barColors20;
xticks(data20.StartMin);
set(gca,'XTickLabel',xlabels20,'XTickLabelRotation',45);
xlabel('Time window relative to Parturition (min)');
ylabel('% of time in bin');
title('% Time NESTING (20-min bins)');
grid on; box off; set(gca,'Color','none','LineWidth',1,'FontSize',11);
xlim([startMinutes, max(data20.StartMin)+binWidth20]); ylim([0 100]);

% Phase summary + stars
phaseSummaryPlotWithStats(pctNesting10(:), data10.Phase, ...
    color_baseline, color_pre, color_par, color_post, ...
    '% Time NESTING (10-min, per phase)', ...
    outdir, 'PhaseSummary_Nesting');

%% --- HELPERS -------------------------------------------------------------
function c = pickColorWithBaseline(startMin, bStart, bEnd, pStart, pEnd, cBase, cPre, cPar, cPost)
    % Baseline override first
    if startMin >= bStart && startMin < bEnd
        c = cBase;
    elseif startMin < pStart
        % Pre is now effectively [-60..0) because baseline is white
        c = cPre;
    elseif startMin >= pStart && startMin < pEnd
        c = cPar;
    else
        c = cPost;
    end
end

function phaseSummaryPlotWithStats(y, Phase, ...
    color_baseline, color_pre, color_par, color_post, ...
    titleStr, outdir, baseFilename)

    % Remove NaNs / undefined
    valid = ~isnan(y) & ~isundefined(Phase);
    y = y(valid);
    g = Phase(valid);

    if numel(y) < 2
        warning('[%s] Not enough data points for stats.', titleStr);
        return;
    end

    % --- RUN OMNIBUS TEST (ANOVA or KW) ----------------------------------
    % First ANOVA
    [p_anova, ~, stats_anova] = anova1(y, g, 'off');
    gnames = cellstr(stats_anova.gnames);  % group labels in stats order

    % Rebuild grouping index so that indices 1..nGroup match gnames
    gcat = categorical(cellstr(g), gnames);
    groupIdx = grp2idx(gcat);
    nGroup   = numel(gnames);

    % Residuals for normality
    yhat = zeros(size(y));
    for i = 1:nGroup
        vals = y(groupIdx == i);
        if isempty(vals)
            continue;
        end
        m = mean(vals, 'omitnan');
        yhat(groupIdx == i) = m;
    end
    resid = y - yhat;
    resid = resid(~isnan(resid));

    if numel(resid) >= 5
        [h_norm, p_norm] = lillietest(resid);
    else
        h_norm = 1;    % force nonparametric if too few
        p_norm = NaN;
    end

    use_nonparam = (h_norm ~= 0);

    if ~use_nonparam
        fprintf('\n[%s] Residuals ~normal (p = %.3g). Using ANOVA.\n', titleStr, p_norm);
        fprintf('[%s] ANOVA p = %.3g\n', titleStr, p_anova);
        c = multcompare(stats_anova, 'CType','tukey', 'Display','off');
    else
        fprintf('\n[%s] Residuals non-normal (p = %.3g). Using Kruskal–Wallis.\n', titleStr, p_norm);
        [p_kw, ~, stats_kw] = kruskalwallis(y, g, 'off');
        fprintf('[%s] Kruskal–Wallis p = %.3g\n', titleStr, p_kw);

        c = multcompare(stats_kw, 'Display','off');

        % Use KW group names & enforce that order on g
        gnames = cellstr(stats_kw.gnames);
        gcat   = categorical(cellstr(g), gnames);
        groupIdx = grp2idx(gcat);
        nGroup   = numel(gnames);
    end

    % --- GROUP Ns ---------------------------------------------------------
    groupNs = zeros(nGroup,1);
    for i = 1:nGroup
        groupNs(i) = sum(groupIdx == i);
    end
    fprintf('[%s] Group Ns:\n', titleStr);
    for i = 1:nGroup
        fprintf('  %s: n = %d\n', gnames{i}, groupNs(i));
    end

    % --- Pairwise p-value matrix in STATS order (1..nGroup) --------------
    pairP = nan(nGroup);
    for k = 1:size(c,1)
        i1   = c(k,1);
        i2   = c(k,2);
        pval = c(k,6);
        pairP(i1,i2) = pval;
        pairP(i2,i1) = pval;
    end

    % Also print as a proper table so you can inspect easily
    pairTable = array2table(pairP, 'VariableNames', gnames, 'RowNames', gnames);
    fprintf('\n[%s] Pairwise p-values:\n', titleStr);
    disp(pairTable);

    % --- Means/SEMs in STATS order ---------------------------------------
    meanY = zeros(nGroup,1);
    semY  = zeros(nGroup,1);
    for i = 1:nGroup
        vals = y(groupIdx == i);
        meanY(i) = mean(vals, 'omitnan');
        semY(i)  = std(vals, 'omitnan') / max(1, sqrt(sum(~isnan(vals))));
    end

    % --- Plot bar summary in STATS order ---------------------------------
    figure('Color','w');
    bh = bar(1:nGroup, meanY, 'FaceColor','flat', ...
        'EdgeColor','k', 'LineWidth',1.2);

    % Color bars by phase name
    cPhase = zeros(nGroup,3);
    for i = 1:nGroup
        switch strtrim(gnames{i})
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
    errorbar(1:nGroup, meanY, semY, semY, ...
        'k', 'LineStyle','none', 'LineWidth',1.2);

    set(gca,'XTick',1:nGroup,'XTickLabel',gnames,'FontSize',11);
    xlabel('Phase');
    ylabel('% of time in bin');
    title(titleStr);
    box off;
    grid on;
    ax = gca; ax.GridColor = [0.8 0.8 0.8];

    % --- Significant pair connectors only --------------------------------
    alpha_sig = 0.05;

    yBase = max(meanY + semY);
    yStep = range(meanY + semY);
    if yStep == 0
        yStep = max(meanY + semY) * 0.2 + 1;
    end
    yCurrent = yBase + 0.15*yStep;

    pairs  = nchoosek(1:nGroup,2);
    pairPs = arrayfun(@(k) pairP(pairs(k,1), pairs(k,2)), 1:size(pairs,1));
    [~, order] = sort(pairPs, 'ascend', 'MissingPlacement','last');
    pairs  = pairs(order,:);
    pairPs = pairPs(order);

    nSig = 0;
    for kk = 1:size(pairs,1)
        i1   = pairs(kk,1);
        i2   = pairs(kk,2);
        pval = pairPs(kk);

        if isnan(pval) || pval >= alpha_sig
            continue;   % only draw significant ones
        end
        nSig = nSig + 1;

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

    if nSig == 0
        fprintf('[%s] No pairwise comparisons with p < %.3f => no stars drawn.\n', ...
                titleStr, alpha_sig);
    end

    if nargin >= 8 && ~isempty(outdir) && ~isempty(baseFilename)
        save_if_missing(gcf, fullfile(outdir, baseFilename));
    end
end


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
