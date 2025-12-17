%% MULTIMODAL TIMELINE: FILM STRIP + PROBABILISTIC SHADING

%% --- USER SETTINGS ---
refTime = datetime('21-Oct-2025 17:56:55','InputFormat','dd-MMM-yyyy HH:mm:ss');

% Time window of interest (relative to refTime, in minutes)
startMinutes = -120;     % baseline start
endMinutes   = 0;        % up to parturition (0)

% Phase boundaries (purely for visual guides)
baseline_start        = -120;
baseline_end          = -60;
pre_parturition_start = -60;
parturition_time      = 0;

% Probabilistic shading bin width (minutes)
probBinWidth = 5;   % e.g. 5-min bins

% File names
callsFile = fullfile('..','Data','Detections_Audible','all_calls_timeseries.csv');
annotsFile = fullfile('..', 'Data', 'Annotations', 'timeseries_annotations.csv');

% Column name assumptions (change if your CSV headers differ)
timeCol_calls   = 'EventDateTime';   % timestamp in calls table
freqCol_calls   = 'Frequency_kHz';   % call frequency
timeCol_annots  = 'EventDateTime';   % timestamp in annotation table
activeVarName   = 'active';          % generic activity flag (excluded from behaviors)

% Call duration settings
durCol_calls        = 'Duration_s';   % duration column in all_calls_timeseries
durationThreshold_s = 0.1;           % seconds; <= short, > long

% Visual width/height of call strips
stripWidth_min   = 0.3;              % width of each call strip in minutes (purely visual)
stripHeight_kHz  = 4;                % height in kHz for each strip (purely visual)

% Colors for short vs long calls
shortColor = [0.3 0.3 0.3];          % dark gray for short calls
longColor  = [0.9 0.3 0.3];          % reddish for long calls


%% --- LOAD DATA ---
T_calls  = readtable(callsFile);
T_annots = readtable(annotsFile);

% Ensure datetime
T_calls.(timeCol_calls)   = datetime(T_calls.(timeCol_calls),   'InputFormat','dd-MMM-yyyy HH:mm:ss');
T_annots.(timeCol_annots) = datetime(T_annots.(timeCol_annots), 'InputFormat','dd-MMM-yyyy HH:mm:ss');

% Compute Minutes relative to refTime
T_calls.Minutes  = minutes(T_calls.(timeCol_calls)   - refTime);
T_annots.Minutes = minutes(T_annots.(timeCol_annots) - refTime);

% Restrict to analysis window
callsMask  = (T_calls.Minutes  >= startMinutes) & (T_calls.Minutes  <= endMinutes);
annotsMask = (T_annots.Minutes >= startMinutes) & (T_annots.Minutes <= endMinutes);

T_calls_win  = T_calls(callsMask, :);
T_annots_win = T_annots(annotsMask, :);

%% --- IDENTIFY BEHAVIOR VARIABLES (WITH CUSTOM EXCEPTIONS) ---
varNamesAnn = T_annots_win.Properties.VariableNames;

% Explicit sets
forcedBehavior      = {'nest_entry','nest_exit','nesting','circling'};
forcedLocationOnly  = {'inside_nest','outside_nest'};
forcedNonBehavior   = {'parturition'};  % treat as phase marker, not behavior

% Core non-behavior columns
coreNonBehavior = ismember(varNamesAnn, {timeCol_annots, 'Minutes'});

% Generic nest-like detection
isNestLike = contains(varNamesAnn, 'nest', 'IgnoreCase', true);

% Active is never a behavior in this context
isActiveVar = strcmp(varNamesAnn, activeVarName);

% Columns that are explicitly not behaviors
isForcedNonBehavior = ismember(varNamesAnn, [forcedLocationOnly, forcedNonBehavior]);

% First pass: behavior candidates are those that are not
%   - time/Minutes
%   - active
%   - generic nest-like
%   - explicitly forced non-behavior
isBehaviorCandidate = ~(coreNonBehavior | isActiveVar | isNestLike | isForcedNonBehavior);

% Start behaviorVars from these candidates
behaviorVars = varNamesAnn(isBehaviorCandidate);

% Now *force in* specific behaviors even if they contain "nest"
for b = 1:numel(forcedBehavior)
    vname = forcedBehavior{b};
    if ismember(vname, varNamesAnn) && ~ismember(vname, behaviorVars)
        behaviorVars{end+1} = vname; %#ok<AGROW>
    end
end

% Remove any duplicates just in case
behaviorVars = unique(behaviorVars, 'stable');

nBeh = numel(behaviorVars);
if nBeh == 0
    error('No behavior variables detected. Check the forcedBehavior lists / column names.');
end


%% ========================================================================
%  FIGURE 1: TEMPORAL FILM STRIP
%  - Top: squeaks (time vs frequency)
%  - Bottom: behavioral bouts, 1 lane per behavior
% ========================================================================

figure('Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

%% --- TOP PANEL: CALL STRIPS (SHORT vs LONG) ---
nexttile;
hold on;

if ~isempty(T_calls_win)
    % Extract call info
    tCall   = T_calls_win.Minutes;
    fCall   = T_calls_win.(freqCol_calls);
    durCall = T_calls_win.(durCol_calls);

    % Short vs long classification
    isShort = durCall <= durationThreshold_s;
    isLong  = durCall >  durationThreshold_s;

    % Set frequency axis limits based on data
    freqMin = min(fCall) - 1;
    freqMax = max(fCall) + 1;
    ylim([freqMin freqMax]);
else
    % Fallback if somehow no calls in window
    freqMin = 0;
    freqMax = 80;
    ylim([freqMin freqMax]);
end

% Phase shading (using current y-limits)
yl = ylim;

% Baseline shading
if baseline_start < endMinutes && baseline_end > startMinutes
    x1 = max(baseline_start, startMinutes);
    x2 = min(baseline_end,   endMinutes);
    patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], ...
          'EdgeColor','none','FaceAlpha',0.3);
end

% Pre-parturient shading
if pre_parturition_start < endMinutes && parturition_time > startMinutes
    x1 = max(pre_parturition_start, startMinutes);
    x2 = min(parturition_time,      endMinutes);
    patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.9 1.0], ...
          'EdgeColor','none','FaceAlpha',0.3);
end

% Draw call strips on top of shading
if ~isempty(T_calls_win)
    % SHORT calls
    shortIdx = find(isShort);
    for ii = shortIdx'
        x1 = tCall(ii) - stripWidth_min/2;
        x2 = tCall(ii) + stripWidth_min/2;
        y1 = fCall(ii) - stripHeight_kHz/2;
        y2 = fCall(ii) + stripHeight_kHz/2;
        rectangle('Position',[x1 y1 x2-x1 y2-y1], ...
                  'FaceColor', [shortColor 0.9], ...
                  'EdgeColor','none');
    end

    % LONG calls
    longIdx = find(isLong);
    for ii = longIdx'
        x1 = tCall(ii) - stripWidth_min/2;
        x2 = tCall(ii) + stripWidth_min/2;
        y1 = fCall(ii) - stripHeight_kHz/2;
        y2 = fCall(ii) + stripHeight_kHz/2;
        rectangle('Position',[x1 y1 x2-x1 y2-y1], ...
                  'FaceColor', [longColor 0.9], ...
                  'EdgeColor','none');
    end

    % Dummy patches for legend
    hShort = patch(NaN, NaN, shortColor, 'FaceAlpha',0.9, 'EdgeColor','none');
    hLong  = patch(NaN, NaN, longColor,  'FaceAlpha',0.9, 'EdgeColor','none');
    legend([hShort hLong], {'Short calls', 'Long calls'}, 'Location','best');
end

xline(0,'r--','LineWidth',1.5);   % parturition

xlabel('Time (min relative to parturition)');
ylabel('Call frequency (kHz)');
title('Vocalizations (call strips: short vs long)');
xlim([startMinutes endMinutes]);
box on;

%% --- BOTTOM PANEL: BEHAVIOR LANES ---
nexttile;
hold on;

cmap = lines(nBeh);   % one color per behavior
yMargin = 0.35;

for b = 1:nBeh
    vname = behaviorVars{b};
    vals  = T_annots_win.(vname);
    
    % Treat numeric/logical nonzero as "behavior present"
    isOn = false(height(T_annots_win),1);
    if islogical(vals)
        isOn = vals;
    elseif isnumeric(vals)
        isOn = vals ~= 0;
    end
    
    if ~any(isOn)
        continue;  % no bouts for this behavior
    end
    
    % Find contiguous segments where isOn == true
    onIdx = find(isOn);
    
    % We work with the full mask to detect transitions
    m = isOn;
    dm = diff([false; m; false]);   % pad for edge transitions
    startIdx = find(dm == 1);
    endIdx   = find(dm == -1) - 1;
    
    for k = 1:numel(startIdx)
        sIdx = startIdx(k);
        eIdx = endIdx(k);
        tStart = T_annots_win.Minutes(sIdx);
        tEnd   = T_annots_win.Minutes(eIdx);
        
        % Clip to our plotting window
        tStart = max(tStart, startMinutes);
        tEnd   = min(tEnd,   endMinutes);
        if tEnd <= tStart
            continue;
        end
        
        yCenter = b;
        rectPos = [tStart, yCenter - yMargin,  tEnd - tStart,  2*yMargin];
        rectangle('Position', rectPos, ...
                  'FaceColor', [cmap(b,:) 0.6], ...
                  'EdgeColor', 'none');
    end
end

% Phase shading (behind rectangles)
yl_b = [0.5 nBeh+0.5];
uistack(gca,'bottom');  % keep axes but we’ll redraw patches then behavior

% Draw phase shading AFTER resetting hold
hold off;
hold on;
% Baseline shading
if baseline_start < endMinutes && baseline_end > startMinutes
    x1 = max(baseline_start, startMinutes);
    x2 = min(baseline_end, endMinutes);
    patch([x1 x2 x2 x1], [yl_b(1) yl_b(1) yl_b(2) yl_b(2)], [0.9 0.9 0.9], ...
          'EdgeColor','none','FaceAlpha',0.3);
end
% Pre-parturient shading
if pre_parturition_start < endMinutes && parturition_time > startMinutes
    x1 = max(pre_parturition_start, startMinutes);
    x2 = min(parturition_time, endMinutes);
    patch([x1 x2 x2 x1], [yl_b(1) yl_b(1) yl_b(2) yl_b(2)], [0.8 0.9 1.0], ...
          'EdgeColor','none','FaceAlpha',0.3);
end

% Re-draw behavior rectangles on top of shading
for b = 1:nBeh
    vname = behaviorVars{b};
    vals  = T_annots_win.(vname);
    
    isOn = false(height(T_annots_win),1);
    if islogical(vals)
        isOn = vals;
    elseif isnumeric(vals)
        isOn = vals ~= 0;
    end
    
    if ~any(isOn)
        continue;
    end
    
    m = isOn;
    dm = diff([false; m; false]);
    startIdx = find(dm == 1);
    endIdx   = find(dm == -1) - 1;
    
    for k = 1:numel(startIdx)
        sIdx = startIdx(k);
        eIdx = endIdx(k);
        tStart = T_annots_win.Minutes(sIdx);
        tEnd   = T_annots_win.Minutes(eIdx);
        tStart = max(tStart, startMinutes);
        tEnd   = min(tEnd,   endMinutes);
        if tEnd <= tStart, continue; end
        
        yCenter = b;
        rectPos = [tStart, yCenter - yMargin,  tEnd - tStart,  2*yMargin];
        rectangle('Position', rectPos, ...
                  'FaceColor', [cmap(b,:) 0.8], ...
                  'EdgeColor', 'none');
    end
end

xline(0,'r--','LineWidth',1.5);
xlim([startMinutes endMinutes]);
ylim(yl_b);
set(gca,'YTick',1:nBeh,'YTickLabel',behaviorVars);
xlabel('Time (min relative to parturition)');
ylabel('Behavior');
title('Behavior bouts (film strip)');
box on;

%% ========================================================================
%  FIGURE 2: PROBABILISTIC SHADING
%  - Heatmap of behavior occupancy over time
%  - Calls overlaid as markers
% ========================================================================

% Define time bins
binEdges   = startMinutes:probBinWidth:endMinutes;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
nBins      = numel(binCenters);

% Compute occupancy: fraction of annotation samples in each bin where behavior == 1
occupancy = nan(nBeh, nBins);

for b = 1:nBeh
    vname = behaviorVars{b};
    vals  = T_annots_win.(vname);
    
    isOn = false(height(T_annots_win),1);
    if islogical(vals)
        isOn = vals;
    elseif isnumeric(vals)
        isOn = vals ~= 0;
    end
    
    t = T_annots_win.Minutes;
    for k = 1:nBins
        idx = (t >= binEdges(k)) & (t < binEdges(k+1));
        if any(idx)
            occupancy(b,k) = sum(isOn(idx)) / sum(idx);
        else
            occupancy(b,k) = NaN;  % no data for that bin
        end
    end
end

figure('Color','w');
imagesc(binCenters, 1:nBeh, occupancy);
set(gca,'YDir','normal');
set(gca,'YTick',1:nBeh,'YTickLabel',behaviorVars);
xlabel('Time (min relative to parturition)');
ylabel('Behavior');
title(sprintf('Behavior occupancy (%.1f-min bins)', probBinWidth));
colormap(parula);
c = colorbar;
ylabel(c,'Fraction of annotation samples with behavior on','Rotation',90);

xlim([startMinutes endMinutes]);

hold on;
% Overlay calls as markers above the top behavior row
if ~isempty(T_calls_win)
    y_calls = ones(height(T_calls_win),1)*(nBeh + 0.5);
    plot(T_calls_win.Minutes, y_calls, 'kv', 'MarkerFaceColor','k', 'MarkerSize',4);
    ylim([0.5 nBeh + 1]);
else
    ylim([0.5 nBeh + 0.5]);
end

xline(0,'r--','LineWidth',1.5);  % parturition
box on;
