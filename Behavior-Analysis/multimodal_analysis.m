


%% --- USER SETTINGS ---
refTime = datetime('21-Oct-2025 17:56:55','InputFormat','dd-MMM-yyyy HH:mm:ss');

% Pre-parturition window (relative to refTime, in minutes)
pre_parturition_start = -60;
parturition_start     = 0;

% File names
callsFile = fullfile('..','Data','Detections_Audible','all_calls_timeseries.csv');
annotsFile = fullfile('..', 'Data', 'Annotations', 'timeseries_annotations.csv');

% Column name assumptions (edit if needed)
timeCol_calls  = 'EventDateTime';  % timestamp column in all_calls_timeseries
timeCol_annots = 'EventDateTime';  % timestamp column in timeseries_annotations

activeVarName = 'active';          % activity flag column

% Explicit location variables if they exist (we also auto-detect "nest" names)
locationVars_manual = {'in_nest','nest_entrance','off_nest','nest_inside','nest_edge'};

% Max allowed time difference between call and annotation (sec)
% [] = always-nearest
maxTimeDiffSec = [];  % keep [] unless you want a hard cutoff

%% --- LOAD DATA ---
T_calls  = readtable(callsFile);
T_annots = readtable(annotsFile);

% Ensure datetime
T_calls.(timeCol_calls)   = datetime(T_calls.(timeCol_calls),   'InputFormat','dd-MMM-yyyy HH:mm:ss');
T_annots.(timeCol_annots) = datetime(T_annots.(timeCol_annots), 'InputFormat','dd-MMM-yyyy HH:mm:ss');

% Compute Minutes relative to refTime
T_calls.Minutes  = minutes(T_calls.(timeCol_calls)   - refTime);
T_annots.Minutes = minutes(T_annots.(timeCol_annots) - refTime);

%% --- FOCUS ON PRE-PARTURITION CALLS ONLY ---
preMask = (T_calls.Minutes >= pre_parturition_start) & (T_calls.Minutes < parturition_start);
T_calls_pre = T_calls(preMask, :);

if isempty(T_calls_pre)
    warning('No calls found in the pre-parturition window.');
    return;
end

%% --- CONVERT TO TIMETABLES & ALIGN ANNOTATIONS TO CALL TIMES ---
ttCalls = table2timetable(T_calls_pre, 'RowTimes', T_calls_pre.(timeCol_calls));
ttAnn   = table2timetable(T_annots,    'RowTimes', T_annots.(timeCol_annots));

% Retime annotations onto call times using nearest neighbor
ttAnnOnCalls = retime(ttAnn, ttCalls.Time, 'nearest');

%% --- MERGE CALLS + ANNOTATIONS (safe version) ---
ttSync = ttCalls;

for v = annVars
    vname = v{1};
    if ismember(vname, ttAnnOnCalls.Properties.VariableNames)
        ttSync.(vname) = ttAnnOnCalls.(vname);
    else
        warning('Annotation variable %s not found in ttAnnOnCalls, skipping.', vname);
    end
end


%% --- BEHAVIOR LABEL ASSIGNMENT PER CALL ---
varNamesAnn = T_annots.Properties.VariableNames;

% Identify core non-behavior columns
coreNonBehavior = ismember(varNamesAnn, {timeCol_annots, 'Minutes'});

% Active column?
hasActive = ismember(activeVarName, varNamesAnn);

% Identify location-like columns
isManualLocation = ismember(varNamesAnn, locationVars_manual);
isNestLike       = contains(varNamesAnn, 'nest', 'IgnoreCase', true);  % anything with "nest"

isLocation = isManualLocation | isNestLike;

% Behavior candidates = not time, not location
isBehaviorCandidate = ~(coreNonBehavior | isLocation);
if hasActive
    % active is treated separately and never as a behavior label
    isBehaviorCandidate = isBehaviorCandidate & ~strcmp(varNamesAnn, activeVarName);
end

behaviorVars = varNamesAnn(isBehaviorCandidate);

% Convenience flag for circling if it exists
hasCircling = ismember('circling', behaviorVars);

nCalls = height(ttSync);
behaviorLabels = strings(nCalls,1);

for i = 1:nCalls
    % --- get the annotation values for this call time ---
    rowSync = ttSync(i, :);
    
    % Detect if we have *any* annotation info
    annotHasData = false;
    for v = behaviorVars
        if ismember(v{1}, ttSync.Properties.VariableNames)
            val = rowSync.(v{1});
            if ~all(ismissing(val))
                annotHasData = true;
                break;
            end
        end
    end
    if hasActive && ismember(activeVarName, ttSync.Properties.VariableNames)
        valA = rowSync.(activeVarName);
        if ~all(ismissing(valA))
            annotHasData = true;
        end
    end
    
    if ~annotHasData
        behaviorLabels(i) = "no annotation";
        continue;
    end
    
    % (Optional) max time difference cutoff – skipped here because we don't
    % have the original annotation time after retime. If you want this,
    % you’ll need to store the offset before retime.
    
    % --- read active flag ---
    isActive = false;
    if hasActive && ismember(activeVarName, ttSync.Properties.VariableNames)
        valA = rowSync.(activeVarName);
        if islogical(valA) || isnumeric(valA)
            isActive = logical(valA);
        end
    end
    
    % --- which behavior variables (non-location) are true? ---
    behTrue = false(1, numel(behaviorVars));
    for b = 1:numel(behaviorVars)
        vname = behaviorVars{b};
        if ~ismember(vname, ttSync.Properties.VariableNames)
            behTrue(b) = false;
            continue;
        end
        val = rowSync.(vname);
        if islogical(val) || isnumeric(val)
            behTrue(b) = logical(val);
        else
            behTrue(b) = false;
        end
    end
    
    % --- APPLY YOUR RULE ---
    % 1) If active = 1 and NO behavior vars are true -> "undefined activity"
    % 2) Else if circling is true -> "circling"
    % 3) Else if any other behavior true -> name of that behavior
    % 4) Else if nothing true but annotations exist -> "no behavior/other"
    
    if isActive && ~any(behTrue)
        behaviorLabels(i) = "undefined activity";
    else
        % circling gets priority if present
        if hasCircling
            idxC = find(strcmp(behaviorVars, 'circling'));
            if ~isempty(idxC) && behTrue(idxC)
                behaviorLabels(i) = "circling";
                continue;
            end
        end
        
        idxTrue = find(behTrue);
        if ~isempty(idxTrue)
            % take the first true behavior
            vname = behaviorVars{idxTrue(1)};
            behaviorLabels(i) = string(vname);
        else
            behaviorLabels(i) = "no behavior/other";
        end
    end
end

%% --- SUMMARY TABLE & PIE CHART (BEHAVIORS ONLY) ---
% We keep "undefined activity" and real behaviors,
% exclude pure "no annotation" from the pie.
pieMask = behaviorLabels ~= "no annotation";
labelsPie = behaviorLabels(pieMask);

if isempty(labelsPie)
    warning('No annotated behaviors found for pre-parturition calls.');
    return;
end

[uniqueLabels, ~, labelIdx] = unique(labelsPie);
counts = accumarray(labelIdx, 1);

% Sort for nicer ordering
[countsSorted, sortIdx] = sort(counts, 'descend');
labelsSorted = uniqueLabels(sortIdx);

% Summary table on full set (including no annotation)
[allLabels, ~, allIdx] = unique(behaviorLabels);
allCounts = accumarray(allIdx, 1);
summaryTable = table(allLabels, allCounts, ...
    allCounts / sum(allCounts) * 100, ...
    'VariableNames', {'Behavior','Count','Percent'});
disp('--- Call-aligned behaviors in PRE-PARTURITION ---');
disp(summaryTable);

% Pie chart
figure('Color','w');
pie(countsSorted, labelsSorted);
title('What is the mouse doing when it squeaks? (Pre-parturition – behaviors only)');
axis equal;
