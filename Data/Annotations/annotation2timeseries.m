function annotation2timeseries()
%% --- USER SETTINGS REMEMBER TO UPDATE EACH FIELD ---
baseStart   = datetime(2025,12,12,09,11,34,483);  % File 1 starts here
sessionDurS = 3600;                           % Each file is 1 hour
dt          = seconds(1);                     % Output sampling
outCSV      = 'timeseries_m5_{parturition}_annotations.csv';

files = { ...
 'Raw data-EphysBirthManualScoring-Trial     3.xlsx'
 'Raw data-EphysBirthManualScoring-Trial     4.xlsx'
 'Raw data-EphysBirthManualScoring-Trial     5.xlsx'
 'Raw data-EphysBirthManualScoring-Trial     6.xlsx'
 'Raw data-EphysBirthManualScoring-Trial     7.xlsx'};

% Column hints (case-insensitive “contains” matching)
col.time    = ["TrialTime","Trial Time","onset","start","date","datetime"];
col.label   = ["Behavior","behaviour","label","code","event"];
col.state   = "Event";  % values like 'state start' / 'state stop'

% Keys
KEY_ENTRY   = "nest entry";
KEY_EXIT    = "nest exit";
KEY_STILL   = "still";             % inactive when inside
KEY_ACTIVE  = "active";            % explicit 'Active (in nest)'
IN_NEST_TAG = "(in nest)";         % any label containing this is an in-nest behavior

%% --- LOAD ALL, CLEAN PAST HEADER, RECONSTRUCT INTERVALS ---
allIntervals = table(); % columns: start_dt, end_dt, label
entryStarts  = datetime.empty(0,1);
exitStarts   = datetime.empty(0,1);

for i = 1:numel(files)
    file = files{i};
    fileT0 = baseStart + seconds((i-1)*sessionDurS);

    Traw = readtable(file,'TextType','string','VariableNamingRule','preserve');
    v = string(Traw.Properties.VariableNames);

    cTime  = pickCol(v, col.time);
    cLabel = pickCol(v, col.label);
    cState = pickCol(v, col.state);
    if cTime=="" || cLabel=="" || cState==""
        error('Could not detect Time/Label/State columns in %s', file);
    end

    % --- strip header junk: keep only rows with 'state start/stop'
    st = lower(strtrim(Traw.(cState)));
    keepRows = contains(st,"state start") | contains(st,"state stop");
    T = Traw(keepRows, :);

    tcol = T.(cTime);                           % numeric seconds or datetimes
    lab  = lower(strtrim(T.(cLabel)));
    st   = lower(strtrim(T.(cState)));

    % absolute time
    if isdatetime(tcol)
        tAbs = tcol;
    else
        tAbs = fileT0 + seconds(double(tcol));
    end

    % sort
    [tAbs, idx] = sort(tAbs);
    lab = lab(idx); st = st(idx);

    % collect entry/exit (state start only)
    entryStarts = [entryStarts; tAbs(contains(lab,KEY_ENTRY) & contains(st,"start"))]; %#ok<AGROW>
    exitStarts  = [exitStarts;  tAbs(contains(lab,KEY_EXIT)  & contains(st,"start"))]; %#ok<AGROW>

    % rebuild intervals for every distinct label
    uLabs = unique(lab);
    for L = uLabs'
        rows = lab==L;
        tL = tAbs(rows); sL = st(rows);

        iStart = find(contains(sL,"start"));
        iStop  = find(contains(sL,"stop"));

        Z = table();
        for k = 1:numel(iStart)
            sIdx = iStart(k);
            eIdx = iStop(find(iStop > sIdx, 1, 'first'));
            if isempty(eIdx)
                stopTime = fileT0 + seconds(sessionDurS); % cap to file end
            else
                stopTime = tL(eIdx);
            end
            startTime = tL(sIdx);
            if startTime < stopTime
                Z = [Z; table(startTime, stopTime, L, ...
                    'VariableNames',{'start_dt','end_dt','label'})]; %#ok<AGROW>
            end
        end
        allIntervals = [allIntervals; Z]; %#ok<AGROW>
    end
end

% Clean + sort
if ~isempty(allIntervals)
    ok = ~isnat(allIntervals.start_dt) & ~isnat(allIntervals.end_dt) & ...
         (allIntervals.start_dt < allIntervals.end_dt);
    allIntervals = sortrows(allIntervals(ok,:), 'start_dt');
end
entryStarts = sort(entryStarts); exitStarts = sort(exitStarts);

%% --- UNIFIED TIMEBASE ---
t0 = baseStart;
t1 = baseStart + seconds(numel(files)*sessionDurS);
timevec = (t0:dt:t1)';

S = struct();
S.EventDateTime = timevec;

%% --- INSIDE/OUTSIDE with correct INITIAL STATE inference ---
% Determine initial state at t0
inside0 = false;
if ~isempty(entryStarts) || ~isempty(exitStarts)
    firstEntry = iif(~isempty(entryStarts), entryStarts(1), datetime(Inf,1,1));
    firstExit  = iif(~isempty(exitStarts),  exitStarts(1),  datetime(Inf,1,1));
    if firstExit < firstEntry
        % First thing that happens is an EXIT -> must have been INSIDE before it
        inside0 = true;
    else
        % First is an ENTRY -> was OUTSIDE before it
        inside0 = false;
    end
else
    % No entry/exit markers at all: infer from behaviors active at t0
    inside0 = anyIntervalActiveAt(allIntervals, t0, IN_NEST_TAG);
end

% Build inside timeline by toggling on entry/exit starts
inside = false(size(timevec));
state  = inside0;
ePtr = 1; xPtr = 1;
for k = 1:numel(timevec)
    tk = timevec(k);
    % apply any entries/exits at or before tk
    while ePtr<=numel(entryStarts) && entryStarts(ePtr) <= tk
        state = true;  ePtr = ePtr + 1;
    end
    while xPtr<=numel(exitStarts)  && exitStarts(xPtr)  <= tk
        state = false; xPtr = xPtr + 1;
    end
    inside(k) = state;
end
S.inside_nest  = inside;
S.outside_nest = ~inside;

%% --- SPECIFIC BEHAVIOR COLUMNS (preserved) ---
behLabels = string([]);
if ~isempty(allIntervals), behLabels = unique(allIntervals.label); end
for L = behLabels'
    S.(mk(L)) = false(size(timevec));
end
for r = 1:height(allIntervals)
    m = (timevec >= allIntervals.start_dt(r)) & (timevec < allIntervals.end_dt(r));
    S.(mk(allIntervals.label(r)))(m) = true;
end

%% --- ACTIVE vs INACTIVE in nest ---
isStillCols = behColsMatching(behLabels, KEY_STILL) & behColsMatching(behLabels, IN_NEST_TAG);
isActCols   = behColsMatching(behLabels, KEY_ACTIVE) & behColsMatching(behLabels, IN_NEST_TAG);
isInNestSpec= behColsMatching(behLabels, IN_NEST_TAG);

stillTrack = false(size(timevec));
for L = behLabels(isStillCols)'
    stillTrack = stillTrack | S.(mk(L));
end

explicitActive = false(size(timevec));
for L = behLabels(isActCols)'
    explicitActive = explicitActive | S.(mk(L));
end

% any in-nest specific behavior except 'still'
anyInNestSpecific = false(size(timevec));
for L = behLabels(isInNestSpec & ~isStillCols)'
    anyInNestSpecific = anyInNestSpecific | S.(mk(L));
end

S.inactive_in_nest = S.inside_nest & stillTrack;
S.active_in_nest   = S.inside_nest & (explicitActive | anyInNestSpecific) & ~S.inactive_in_nest;

%% --- GAP FILL (inside only) for active/inactive by midpoint rule ---
gapMask = S.inside_nest & ~(S.inactive_in_nest | S.active_in_nest);
runs = logicalRuns(gapMask);
for rr = 1:size(runs,1)
    i1 = runs(rr,1); i2 = runs(rr,2);
    prevIdx = find(~gapMask(1:i1-1) & S.inside_nest(1:i1-1), 1, 'last');
    nextIdx = i2 + find(~gapMask(i2+1:end) & S.inside_nest(i2+1:end), 1, 'first');
    if ~isempty(prevIdx) && ~isempty(nextIdx)
        mid = floor((i1 + nextIdx)/2);
        if S.inactive_in_nest(prevIdx), S.inactive_in_nest(i1:mid) = true;
        else,                           S.active_in_nest(i1:mid)   = true; end
        if S.inactive_in_nest(nextIdx), S.inactive_in_nest(mid+1:i2) = true;
        else,                           S.active_in_nest(mid+1:i2)   = true; end
    end
end

%% --- WRITE OUTPUT ---
core = ["EventDateTime","inside_nest","outside_nest","inactive_in_nest","active_in_nest"];
Tout = table(S.EventDateTime, S.inside_nest, S.outside_nest, S.inactive_in_nest, S.active_in_nest, ...
             'VariableNames', core);

for f = setdiff(fieldnames(S), {'EventDateTime','inside_nest','outside_nest','inactive_in_nest','active_in_nest'})'
    Tout.(f{1}) = S.(f{1});
end

writetable(Tout, outCSV);
fprintf('Wrote %s (%d rows, %d cols).\n', outCSV, height(Tout), width(Tout));

%% --- HELPERS ---
function name = pickCol(vnames, opts)
    vn = lower(strrep(vnames,'_',' '));
    for o = lower(opts)
        hit = find(contains(vn, o), 1, 'first');
        if ~isempty(hit), name = vnames(hit); return; end
    end
    name = "";
end

function nm = mk(s)
    nm = matlab.lang.makeValidName(lower(regexprep(string(s),'\s+','_')));
end

function mask = behColsMatching(labels, needle)
    labels = lower(string(labels));
    mask = contains(labels, lower(string(needle)));
end

function runs = logicalRuns(m)
    m = m(:)'; dm = diff([false m false]);
    starts = find(dm==1); stops = find(dm==-1)-1;
    runs = [starts(:) stops(:)];
end

function tf = anyIntervalActiveAt(tbl, t, containsStr)
    if isempty(tbl), tf = false; return; end
    tf = any(tbl.start_dt <= t & t < tbl.end_dt & contains(lower(tbl.label), lower(containsStr)));
end

function y = iif(cond, a, b), if cond, y = a; else, y = b; end
end
end