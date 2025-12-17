function concat_calls()
% CONCAT_USV_CALLS_V6
% Uses midpoint frequency = minFreq + freqRange/2 for all frequency measures.
% Otherwise identical to v5 (robust to parentheses filenames and time alignment).

    %---------------------- USER SETTINGS ----------------------%
    rootDir  = fileparts(fileparts(mfilename('fullpath')));
    detectionsFolder = fullfile(rootDir, 'Detections_Audible');

    outCSV_calls   = fullfile(rootDir, 'all_calls_timeseries.csv');
    outCSV_10sum   = fullfile(rootDir, 'calls_10min_summary.csv');
    outCSV_20sum   = fullfile(rootDir, 'calls_20min_summary.csv');
    outCSV_10minTS = fullfile(rootDir, 'calls_10min_timeseries.csv');

    filePrefix = '3rd_mouse';
    slotMinutes = 10;

    windowDate   = [2025 10 21];
    windowStartT = [14 56 55];
    windowEndT   = [18 56 55];
    %-----------------------------------------------------------%

    if ~isfolder(detectionsFolder)
        error('Detections folder not found: %s', detectionsFolder);
    end

    files = dir(fullfile(detectionsFolder, [filePrefix '20*.mat']));
    if isempty(files)
        error('No matching .mat files found in %s', detectionsFolder);
    end

    % Parse datetime and (optional) parentheses interval
    baseRx = '(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-(\d{2})';
    parenRx = '\((\d{2}):(\d{2}):(\d{2})-(\d{2}):(\d{2}):(\d{2})\)';
    fileInfo = struct('name',{files.name},'folder',{files.folder},'startDT',datetime.empty);

    for k = 1:numel(files)
        f = files(k).name;
        t1 = regexp(f, baseRx, 'tokens', 'once');
        if isempty(t1), continue; end
        yyyy=str2double(t1{1}); mm=str2double(t1{2}); dd=str2double(t1{3});
        hh=str2double(t1{4}); mi=str2double(t1{5}); ss=str2double(t1{6});
        baseDT = datetime(yyyy,mm,dd,hh,mi,ss);
        t2 = regexp(f, parenRx, 'tokens', 'once');
        if ~isempty(t2)
            hh2=str2double(t2{1}); mi2=str2double(t2{2}); ss2=str2double(t2{3});
            baseDT = datetime(yyyy,mm,dd,hh2,mi2,ss2);
        end
        fileInfo(k).name = f;
        fileInfo(k).folder = files(k).folder;
        fileInfo(k).startDT = baseDT;
    end
    fileInfo = fileInfo(~arrayfun(@isempty,{fileInfo.startDT}));
    [~,ord] = sort([fileInfo.startDT]);
    fileInfo = fileInfo(ord);

    ws = datetime(windowDate(1),windowDate(2),windowDate(3), ...
                  windowStartT(1),windowStartT(2),windowStartT(3));
    we = datetime(windowDate(1),windowDate(2),windowDate(3), ...
                  windowEndT(1),windowEndT(2),windowEndT(3));
    expectedStarts = ws:minutes(slotMinutes):we;
    nSlots = numel(expectedStarts);
    tol = seconds(2);

    % Assign each file to nearest expected slot
    slotFile = cell(nSlots,1);
    for k = 1:numel(fileInfo)
        [d,idx] = min(abs(expectedStarts - fileInfo(k).startDT));
        if d <= tol
            slotFile{idx} = fileInfo(k);
        end
    end

    % Accumulators
    allDT = datetime.empty(0,1);
    allFreq = []; allDur = []; allFileStart = datetime.empty(0,1);
    slotCount = zeros(nSlots,1); slotMeanF = nan(nSlots,1);
    slotSDF = nan(nSlots,1); slotMissing = false(nSlots,1);

    fprintf('Processing %d slots...\n', nSlots);

    for i = 1:nSlots
        st = expectedStarts(i);
        if isempty(slotFile{i})
            slotMissing(i) = true;
            continue;
        end
        info = slotFile{i};
        S = load(fullfile(info.folder,info.name));
        [T,~] = find_detection_table_with_vector_first_col(S,info.name);
        if isempty(T) || height(T)==0
            slotCount(i)=0; continue;
        end
        [beginSec,minFreq,duration,freqRange,~] = expand_first_col_to_four(T,info.name);
        good = isfinite(beginSec)&beginSec>=0&isfinite(minFreq)&isfinite(duration)&isfinite(freqRange);
        beginSec = beginSec(good); minFreq=minFreq(good);
        duration=duration(good); freqRange=freqRange(good);
        if isempty(beginSec), slotCount(i)=0; continue; end

        % Compute midpoint frequency
        centerFreq = minFreq + freqRange/2;

        % Real absolute call times
        beginDT = info.startDT + seconds(beginSec);
        inwin = beginDT>=ws & beginDT<=we;
        beginDT = beginDT(inwin);
        centerFreq = centerFreq(inwin);
        duration = duration(inwin);

        n = numel(beginDT);
        if n>0
            allDT=[allDT; beginDT(:)];
            allFreq=[allFreq; centerFreq(:)];
            allDur=[allDur; duration(:)];
            allFileStart=[allFileStart; repmat(info.startDT,n,1)];
        end
        slotCount(i)=n;
        if n>0
            slotMeanF(i)=mean(centerFreq);
            slotSDF(i)=std(centerFreq,0);
        end
    end

    %----- CSV 1: per-call -----
    PerCall = table(allFileStart, allDT, allFreq, allDur, ...
        'VariableNames', {'FileStartTime','EventDateTime','Frequency_kHz','Duration_s'});
    [~,ord] = sort(PerCall.EventDateTime);
    PerCall = PerCall(ord,:);
    writetable(PerCall, outCSV_calls);
    fprintf('Wrote per-call CSV: %s\n', outCSV_calls);

    %----- CSV 2: 10-min timeseries -----
    TS10 = table(expectedStarts(:),slotCount,slotMeanF,slotSDF,slotMissing, ...
        'VariableNames',{'SlotStart','CallCount','MeanFreq_kHz','SDFreq_kHz','MissingFile'});
    writetable(TS10,outCSV_10minTS);

    %----- CSV 3–4: summaries -----
    write_summary_bins(PerCall, ws, we, 10, outCSV_10sum);
    write_summary_bins(PerCall, ws, we, 20, outCSV_20sum);
end


%% ---------- Helper functions ----------
function write_summary_bins(PerCall, ws, we, mins, outPath)
    if isempty(PerCall)
        edges = ws:minutes(mins):we+seconds(1);
        centers = edges(1:end-1)+minutes(mins/2);
        Cnt = zeros(numel(centers),1);
        Mf = nan(size(Cnt)); Sf = nan(size(Cnt));
    else
        edges = ws:minutes(mins):we+seconds(1);
        [~,~,bin] = histcounts(PerCall.EventDateTime,edges);
        nb = numel(edges)-1;
        Cnt = accumarray(max(bin,1),1,[nb 1],@sum,0);
        Mf = accumarray(max(bin,1),PerCall.Frequency_kHz,[nb 1],@mean,NaN);
        Sf = nan(nb,1);
        for b=1:nb
            idx = (bin==b);
            if any(idx), Sf(b)=std(PerCall.Frequency_kHz(idx),0); end
        end
        centers = edges(1:end-1)+minutes(mins/2);
    end
    T = table(centers(:),Cnt,Mf,Sf, ...
        'VariableNames',{'BinCenter','CallCount','MeanFreq_kHz','SDFreq_kHz'});
    writetable(T,outPath);
    fprintf('Wrote %d-min summary: %s\n', mins, outPath);
end

function [T, scoreCol] = find_detection_table_with_vector_first_col(S, filename)
    T=[]; scoreCol=0; fn=fieldnames(S);
    candNames=[{'calls'}, fn(:)']; candNames=unique(candNames,'stable');
    for i=1:numel(candNames)
        nm=candNames{i};
        if ~isfield(S,nm), continue; end
        val=S.(nm);
        if istable(val)&&width(val)>=1
            firstCol=val{:,1};
            if iscell(firstCol)
                ok=all(cellfun(@(x)isnumeric(x)&&isvector(x)&&numel(x)>=4,firstCol));
            elseif isnumeric(firstCol)
                ok=ismatrix(firstCol)&&size(firstCol,2)>=4;
            else, ok=false;
            end
            if ok, T=val; return; end
        end
    end
end

function [beginSec,minFreq,duration,freqRange,goodMask] = expand_first_col_to_four(T, filename)
    firstCol=T{:,1};
    n=height(T);
    beginSec=NaN(n,1); minFreq=NaN(n,1); duration=NaN(n,1); freqRange=NaN(n,1);
    if iscell(firstCol)
        for i=1:n
            v=firstCol{i};
            if isnumeric(v)&&isvector(v)&&numel(v)>=4
                v=double(v(:).');
                beginSec(i)=v(1); minFreq(i)=v(2); duration(i)=v(3); freqRange(i)=v(4);
            end
        end
    elseif isnumeric(firstCol)
        beginSec=double(firstCol(:,1));
        minFreq=double(firstCol(:,2));
        duration=double(firstCol(:,3));
        freqRange=double(firstCol(:,4));
    end
    goodMask=isfinite(beginSec)&isfinite(minFreq)&isfinite(duration)&isfinite(freqRange);
    beginSec=beginSec(goodMask); minFreq=minFreq(goodMask);
    duration=duration(goodMask); freqRange=freqRange(goodMask);
end
