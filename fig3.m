% Figure 3: track processing pipeline (raw echogram -> detections ->
% tracks -> world-frame projection) and example trajectories, for a
% window spanning the first few shots of firing period F1.

close all; clear; clc;

here = fileparts(mfilename('fullpath'));
data_dir = fullfile(here, 'data');
addpath(fullfile(here, 'utils'));

%% ===================== load =====================
TRACKS = load(fullfile(data_dir, 'auv', 'tracks.mat'), 'allTracks');
SHOTS = load(fullfile(data_dir, 'airgun', 'shots.mat'), 'shots');

% raw echogram + detections for the shown window 
ECHO = load(fullfile(data_dir, 'auv', 'sp_194843.mat'));
DETS = load(fullfile(data_dir, 'auv', 'detections_194843.mat'), 'detections');

%% ===================== compute =====================
allTracks = TRACKS.allTracks;
shotTimes = sort(SHOTS.shots.time);
rawDets = DETS.detections;

% restrict shot times to the track-data time window (tracks.mat is already cropped to end of Firing 2)
allHitTimes = [];
for i = 1:length(allTracks)
    allHitTimes = [allHitTimes; vertcat(allTracks(i).posLog{:, 5})]; %#ok<AGROW>
end
tData_start = min(allHitTimes);
tData_end = max(allHitTimes);
shotTimes = shotTimes(shotTimes >= tData_start & shotTimes <= tData_end);

t0 = tData_start;

% shown-shot window 
nShow = 5;      % number of shots to show
padPre = 44;    % seconds of context before the first shown shot
padPost = 10;   % seconds of context after the last shown shot
shotSec = seconds(shotTimes - t0);
nShowActual = min(nShow, numel(shotSec));
winStartMin = (shotSec(1) - padPre) / 60;
winEndMin = (shotSec(nShowActual) + padPost) / 60;
shotsInWinMin = shotSec(shotSec >= winStartMin*60 & shotSec <= winEndMin*60) / 60;
fprintf('Shown window: %.2f - %.2f min (%d of %d shots)\n', ...
    winStartMin, winEndMin, nShowActual, numel(shotSec));

rawRlim = [5 60];
rawEcho.Sp = ECHO.Sp;
rawEcho.range = ECHO.range;
rawEcho.time_abs = ECHO.time;
if isempty(rawEcho.time_abs.TimeZone), rawEcho.time_abs.TimeZone = 'UTC'; end

% the 4 example tracks shown in the 2x2 3D grid
exampleTrackIdx = [5866 6171 5750 5744];

%% ===================== plot: processing-chain row + echogram/timeseries + 2x2 3D grid =====================
fig = figure('Units', 'centimeters', 'Position', [2 2 34.5 24.34], 'Color', 'w');
divCmap = buildDivColormap(256);
cmaxDown = 20;   % cm/s, hit-level v_Down colour map limits

% --- top row: 5 processing-chain panels ---
% (1) raw echogram (Sp), (2) raw detections in 2D (3) those detections in 
% 3D, (4) RTS-smoothed tracks in 3D (coloured by track ID), (5) those 
% tracks projected into the world frame (coloured by track ID) with the 
% AUV's own path overlaid in black. Panels 1-3 share one Sp colorbar
cbGapFrac1 = 0.02;
cbFrac1 = 0.07;
imgFrac1 = 1 - cbGapFrac1 - cbFrac1;   % panel 1's image/colorbar split
procBottom = 0.7330;
procHeight = 0.2465;
procGaps = [0.050 0.021 0.012 0.03];   % gaps
wts = [1.4/imgFrac1 1.4 0.75 0.75 2.1];   % panel widths 
totalSpan = 0.9713 - 0.0722;
panelUnit = (totalSpan - sum(procGaps)) / sum(wts);
procW = wts * panelUnit;
procLeftEdges = 0.0722 + [0, cumsum(procW(1:end-1) + procGaps)];
procLeftEdges(5) = procLeftEdges(4) + procW(4) + 0.022;   % panel 5 pulled a bit closer to panel 4
procPos = arrayfun(@(i) [procLeftEdges(i) procBottom procW(i) procHeight], 1:5, 'UniformOutput', false);

% panel 1's own image+colorbar split, so the colorbar stays inside its box
pp1 = procPos{1};
posRawImg = [pp1(1), pp1(2), pp1(3)*imgFrac1, pp1(4)];
posRawCb = [pp1(1)+pp1(3)*(imgFrac1+cbGapFrac1), pp1(2), pp1(3)*cbFrac1, pp1(4)];

axProc = gobjects(5, 1);
for k = 1:5
    if k == 1
        axProc(k) = axes('Position', posRawImg); %#ok<LAXES>
        plotRawEchoPanel(axProc(k), rawEcho, t0, winStartMin, winEndMin, rawRlim, posRawImg, posRawCb);
        continue
    end
    if k == 2
        axProc(k) = axes('Position', procPos{k}); %#ok<LAXES>
        plotDetectionsPanel2D(axProc(k), rawDets, t0, winStartMin, winEndMin, rawRlim, procPos{k});
        continue
    end
    if k == 3
        axProc(k) = axes('Position', procPos{k}); %#ok<LAXES>
        plotDetectionsSp3D(axProc(k), rawDets, t0, winStartMin, winEndMin, procPos{k});
        % capture panel 3's limits/ticks/aspect so panel 4 can match exactly
        panel3Lims = struct('XLim', xlim(axProc(k)), 'YLim', ylim(axProc(k)), 'ZLim', zlim(axProc(k)), ...
            'XTick', get(axProc(k), 'XTick'), 'YTick', get(axProc(k), 'YTick'), 'ZTick', get(axProc(k), 'ZTick'), ...
            'DataAspectRatio', get(axProc(k), 'DataAspectRatio'), ...
            'PlotBoxAspectRatio', get(axProc(k), 'PlotBoxAspectRatio'));
        continue
    end
    if k == 4
        axProc(k) = axes('Position', procPos{k}); %#ok<LAXES>
        plotRTSTracks3D(axProc(k), allTracks, t0, winStartMin, winEndMin, procPos{k}, panel3Lims);
        continue
    end
    axProc(k) = axes('Position', procPos{k}); %#ok<LAXES>
    plotWorldTracks3D(axProc(k), allTracks, t0, winStartMin, winEndMin, procPos{k});
end

% --- left block: echogram (range vs. time) + hit-level vDown time series ---
posE = [0.0722 0.1906 0.4008 0.4601];    % echogram
posV = [0.0722 0.0411 0.4008 0.1328];    % time series, same left/width as echogram
posC = [0.4816 0.1906 0.01243 0.4601];   % colorbar, spans the echogram height

axE = axes('Position', posE);
axT = axes('Position', posV);
plotEchoTimeSeries(axE, axT, allTracks, t0, winStartMin, winEndMin, shotsInWinMin, cmaxDown, divCmap, posE, posC, exampleTrackIdx, rawRlim);

% --- right block: 4 example tracks in 3D, 2x2 grid, square panels ---
gridLeft = 0.5607;
gridWidth = 0.42;
gridBottom = 0.0711;
gridHeight = 0.56;
gw = gridWidth / 2;
gh = gridHeight / 2;
gapX = 0.07;
gapY = 0.10;
figWcm = fig.Position(3);
figHcm = fig.Position(4);
availPanelWnorm = gw - gapX/2;
availPanelHnorm = gh - gapY/2;
panelWnorm = availPanelWnorm;
panelHnorm = panelWnorm * figWcm/figHcm;   % square, if width-limited
if panelHnorm > availPanelHnorm   % height-limited instead
    panelHnorm = availPanelHnorm;
    panelWnorm = panelHnorm * figHcm/figWcm;
end
blockWnorm = 2*panelWnorm + gapX;
blockHnorm = 2*panelHnorm + gapY;
blockLeft = gridLeft + (gridWidth - blockWnorm) / 2;   % re-center 
blockBottom = gridBottom + (gridHeight - blockHnorm) / 2;
pos3D = { ...
    [blockLeft blockBottom+panelHnorm+gapY panelWnorm panelHnorm], ...   % top-left
    [blockLeft+panelWnorm+gapX blockBottom+panelHnorm+gapY panelWnorm panelHnorm], ...   % top-right
    [blockLeft blockBottom panelWnorm panelHnorm], ...   % bottom-left
    [blockLeft+panelWnorm+gapX blockBottom panelWnorm panelHnorm] ...   % bottom-right
    };
ax3D = gobjects(4, 1);
for k = 1:4
    ax3D(k) = axes('Position', pos3D{k}); %#ok<LAXES>
    plot3DTrack(ax3D(k), allTracks, exampleTrackIdx(k), cmaxDown, divCmap);
end

forceFont(fig);
% scale all text 20% larger, applied last so relative sizing is preserved
fsObjs = findall(fig, '-property', 'FontSize');
for i = 1:numel(fsObjs)
    set(fsObjs(i), 'FontSize', get(fsObjs(i), 'FontSize')*1.2);
end
cleanBox(axT);  

%% ===================== helpers =====================

function plotRawEchoPanel(ax, rawEcho, t0, winStartMin, winEndMin, rlim, imgPos, cbPos)
% Processing-chain panel 1: the raw echogram (Sp) alone, before
% detection/tracking. imgPos/cbPos place the image and colorbar explicitly
% so the colorbar stays inside this panel's own box.
    tEchoMin = minutes(rawEcho.time_abs - t0);
    imagesc(ax, tEchoMin, rawEcho.range, rawEcho.Sp);
    colormap(ax, parula);
    clim(ax, [-95 -55]);
    set(ax, 'YDir', 'reverse');
    xlim(ax, [winStartMin winEndMin]);
    ylim(ax, rlim);
    set(ax, 'YTick', 10:10:60);
    xlabel(ax, 'time (min)');
    ylabel(ax, 'range (m)');
    set(ax, 'FontSize', 6, 'FontName', 'Times New Roman');
    box(ax, 'on');
    set(ax, 'Position', imgPos);
    cb = colorbar(ax);
    cb.Label.String = 'S_p (dB re 1 m^2)';
    cb.FontSize = 6;
    cb.Label.FontSize = 6;
    cb.Ticks = -95:10:-55;
    set(ax, 'Position', imgPos);   % colorbar() shrinks ax on creation -- reassert
    set(cb, 'Position', cbPos);
    labelPos = cb.Label.Position;
    cb.Label.Position = [labelPos(1)*0.9, labelPos(2), labelPos(3)];   % pull the label closer to the bar
end


function plotDetectionsPanel2D(ax, rawDets, t0, winStartMin, winEndMin, rlim, imgPos)
% Processing-chain panel 2: raw (pre-tracking) detections, coloured by Sp.
% Shares panel 1's colorbar (same clim/colormap).
    detTimeMin = minutes(datetime(rawDets.time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC') - t0);
    scatter(ax, detTimeMin, rawDets.range, 3.5, rawDets.Sp, 'filled');
    colormap(ax, parula);
    clim(ax, [-95 -55]);
    set(ax, 'YDir', 'reverse');
    xlim(ax, [winStartMin winEndMin]);
    ylim(ax, rlim);
    set(ax, 'YTick', 10:10:60);
    xlabel(ax, 'time (min)');
    ylabel(ax, 'range (m)');
    set(ax, 'FontSize', 6, 'FontName', 'Times New Roman', 'Color', 'w');
    box(ax, 'on');
    set(ax, 'Position', imgPos);
    set(ax, 'Layer', 'top');   % draw the box on top of the data, not underneath
end


function plotDetectionsSp3D(ax, rawDets, t0, winStartMin, winEndMin, imgPos)
% Processing-chain panel 3: the same detections as panel 2, in body-frame
% 3D (x,y,z), coloured by Sp, z axis reversed.
    detTimeMin = minutes(datetime(rawDets.time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC') - t0);
    inWin = detTimeMin >= winStartMin & detTimeMin <= winEndMin;
    scatter3(ax, rawDets.x(inWin), rawDets.y(inWin), rawDets.z(inWin), 5, rawDets.Sp(inWin), 'filled');
    colormap(ax, parula);
    clim(ax, [-95 -55]);
    set(ax, 'Color', [0.9 0.9 0.9], 'FontSize', 6, 'FontName', 'Times New Roman');
    view(ax, [136.5817 14.9388]);
    set(ax, 'ZDir', 'reverse');
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'x (m)');
    ylabel(ax, 'y (m)');
    zlabel(ax, 'z (m)');
    set(ax, 'Position', imgPos);
    zl = zlim(ax);
    set(ax, 'ZTick', ceil(zl(1)/10)*10 : 10 : floor(zl(2)/10)*10);
end


function plotRTSTracks3D(ax, allTracks, t0, winStartMin, winEndMin, imgPos, matchLims)
% Processing-chain panel 4: RTS-smoothed tracks (one line per track,
% coloured by track ID) for tracks overlapping the shown window, in the
% same body-frame 3D view as panel 3. matchLims (optional) is panel 3's
% captured limits/ticks/aspect, forced onto this panel so the two match exactly.
    hold(ax, 'on');
    trackIdxList = [];
    for i = 1:numel(allTracks)
        tt_i_min = minutes(vertcat(allTracks(i).posLog{:, 5}) - t0);
        if all(tt_i_min < winStartMin | tt_i_min > winEndMin), continue; end
        trackIdxList(end+1) = i; %#ok<AGROW>
    end
    nTracks = numel(trackIdxList);
    trackColors = lines(max(nTracks, 1));
    for k = 1:nTracks
        i = trackIdxList(k);
        tt_i_min = minutes(vertcat(allTracks(i).posLog{:, 5}) - t0);
        if isfield(allTracks, 'posSmooth') && ~isempty(allTracks(i).posSmooth)
            pos = allTracks(i).posSmooth;   % RTS-smoothed body frame
        else
            pos = cell2mat(allTracks(i).posLog(:, 4));   % Kalman-updated (filtered) fallback
        end
        inWin = tt_i_min >= winStartMin & tt_i_min <= winEndMin;
        plot3(ax, pos(inWin, 1), pos(inWin, 2), pos(inWin, 3), '-', ...
            'Color', trackColors(k, :), 'LineWidth', 1);
    end
    set(ax, 'Color', [0.9 0.9 0.9], 'FontSize', 6, 'FontName', 'Times New Roman');
    view(ax, [136.5817 14.9388]);
    set(ax, 'ZDir', 'reverse');
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'x (m)');
    ylabel(ax, 'y (m)');
    zlabel(ax, 'z (m)');
    set(ax, 'Position', imgPos);
    if nargin >= 7 && ~isempty(matchLims)
        set(ax, 'DataAspectRatio', matchLims.DataAspectRatio, ...
            'PlotBoxAspectRatio', matchLims.PlotBoxAspectRatio);
        set(ax, 'XLim', matchLims.XLim, 'YLim', matchLims.YLim, 'ZLim', matchLims.ZLim, ...
            'XTick', matchLims.XTick, 'YTick', matchLims.YTick, 'ZTick', matchLims.ZTick);
    end
end


function plotWorldTracks3D(ax, allTracks, t0, winStartMin, winEndMin, imgPos)
% Processing-chain panel 5: the same tracks projected into the world frame
% (posWorld -- East, North, depth, from project.m), coloured by track ID,
% with the AUV's own navigated path overlaid as a thick black line (from
% each track's .auvPos). Everything is plotted relative to the AUV's own
% position at the start of the shown window, since absolute East/North/depth
% are large, not-very-meaningful numbers for a few-minute window.
    hold(ax, 'on');
    trackIdxList = [];
    for i = 1:numel(allTracks)
        tt_i_min = minutes(vertcat(allTracks(i).posLog{:, 5}) - t0);
        if all(tt_i_min < winStartMin | tt_i_min > winEndMin), continue; end
        trackIdxList(end+1) = i; %#ok<AGROW>
    end
    nTracks = numel(trackIdxList);
    trackColors = lines(max(nTracks, 1));

    % pool the AUV's own position across all tracks' hits in the window, to
    % find its start-of-window reference point before plotting anything
    auvT = NaT(0, 1, 'TimeZone', 'UTC');
    auvPos = zeros(0, 3);
    posByTrack = cell(nTracks, 1);
    inWinByTrack = cell(nTracks, 1);
    for k = 1:nTracks
        i = trackIdxList(k);
        tt_i_min = minutes(vertcat(allTracks(i).posLog{:, 5}) - t0);
        inWin = tt_i_min >= winStartMin & tt_i_min <= winEndMin;
        posByTrack{k} = allTracks(i).posWorld;
        inWinByTrack{k} = inWin;

        if isfield(allTracks, 'auvPos') && ~isempty(allTracks(i).auvPos)
            hitTimes = vertcat(allTracks(i).posLog{:, 5});
            auvT = [auvT; hitTimes(inWin)]; %#ok<AGROW>
            auvPos = [auvPos; allTracks(i).auvPos(inWin, :)]; %#ok<AGROW>
        end
    end

    if ~isempty(auvT)
        [auvT, sortIdx] = sort(auvT);
        auvPos = auvPos(sortIdx, :);
        [~, uIdx] = unique(auvT);   % drop exact-duplicate timestamps 
        auvPos = auvPos(uIdx, :);
        refPos = auvPos(1, :);
    else
        refPos = [0 0 0];   % fallback (no AUV path data)
    end

    for k = 1:nTracks
        pos = posByTrack{k} - refPos;
        inWin = inWinByTrack{k};
        plot3(ax, pos(inWin, 1), pos(inWin, 2), pos(inWin, 3), '-', ...
            'Color', trackColors(k, :), 'LineWidth', 1.3);
    end
    if ~isempty(auvT)
        auvPos = auvPos - refPos;
        plot3(ax, auvPos(:, 1), auvPos(:, 2), auvPos(:, 3), '-k', 'LineWidth', 2.5);
    end

    set(ax, 'Color', [0.9 0.9 0.9], 'FontSize', 6, 'FontName', 'Times New Roman');
    view(ax, [243.7936 10.4501]);
    set(ax, 'ZDir', 'reverse');
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '\DeltaE (m)');
    ylabel(ax, '\DeltaN (m)');
    zlabel(ax, '\DeltaD (m)');
    set(ax, 'Position', imgPos);
end


function plotEchoTimeSeries(axE, axT, allTracks, t0, winStartMin, winEndMin, shotsInWinMin, cmax, cmap, posE, posC, exampleTrackIdx, rlim)
% Echogram (range vs. time) in axE, tracks coloured by hit-level v_Down via
% the surface() per-vertex-colour trick, plus a hit-level v_Down time series
% (median +/- IQR, 1 s bins) in axT below it, sharing axE's x-axis. Tracks
% shown in the 3D 2x2 grid (exampleTrackIdx) get a black outline 
    axes(axE); hold(axE, 'on'); box(axE, 'on'); %#ok<LAXES>
    obsTloc = NaT(0, 1, 'TimeZone', 'UTC');
    obsVloc = [];
    chosen = struct('x', {}, 'y', {}, 'vHit', {}, 'gridPos', {});   % chosen tracks drawn in a second pass, on top of everything else

    for i = 1:numel(allTracks)
        tt_i_min = minutes(vertcat(allTracks(i).posLog{:, 5}) - t0);
        if all(tt_i_min < winStartMin | tt_i_min > winEndMin), continue; end

        if isfield(allTracks, 'posSmooth') && ~isempty(allTracks(i).posSmooth)
            pos = allTracks(i).posSmooth;
        else
            pos = cell2mat(allTracks(i).posLog(:, 4));
        end
        rngI = sqrt(sum(pos.^2, 2));

        vpw = allTracks(i).velPtWorld;
        if isempty(vpw) || any(isnan(vpw(:))) || numel(tt_i_min) ~= size(vpw, 1)
            continue
        end
        vHit = vpw(:, 3) * 100;   % cm/s, positive = down

        obsTloc = [obsTloc; vertcat(allTracks(i).posLog{:, 5})]; %#ok<AGROW>
        obsVloc = [obsVloc; vHit(:)]; %#ok<AGROW>

        x = tt_i_min(:);
        y = rngI(:);
        if ismember(i, exampleTrackIdx)
            chosen(end+1) = struct('x', x, 'y', y, 'vHit', vHit(:), 'gridPos', find(exampleTrackIdx == i, 1)); %#ok<AGROW>
            continue
        end
        surface(axE, [x x], [y y], zeros(numel(x), 2), [vHit(:) vHit(:)], ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1.9);
    end

    % chosen tracks' black backing + coloured line, drawn last so they layer above everything else
    tinyEps = (winEndMin - winStartMin) * 0.0005;   % nudges the true endpoint 
    for c = 1:numel(chosen)
        xC = chosen(c).x; yC = chosen(c).y;
        if numel(xC) >= 2
            slopeStart = (yC(2)-yC(1)) / (xC(2)-xC(1));
            slopeEnd = (yC(end)-yC(end-1)) / (xC(end)-xC(end-1));
        else
            slopeStart = 0; slopeEnd = 0;
        end
        xExt = [xC(1)-tinyEps; xC; xC(end)+tinyEps];
        yExt = [yC(1)-slopeStart*tinyEps; yC; yC(end)+slopeEnd*tinyEps];
        hOutline = plot(axE, xExt, yExt, '-k', 'LineWidth', 5.5, 'HandleVisibility', 'off');
        hOutline.LineJoin = 'round';
    end
    for c = 1:numel(chosen)
        surface(axE, [chosen(c).x chosen(c).x], [chosen(c).y chosen(c).y], zeros(numel(chosen(c).x), 2), ...
            [chosen(c).vHit chosen(c).vHit], 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3.2);
    end

    % number each chosen track next to its start point, matching the "N:"
    % label on that track's 3D panel
    timeOffset = (winEndMin - winStartMin) * 0.015;
    rangeOffset = 0.4;   % meters
    for c = 1:numel(chosen)
        xLbl = chosen(c).x(1) - timeOffset;
        if chosen(c).gridPos == 1
            xLbl = xLbl - (winEndMin - winStartMin) * 0.01;
        end
        if chosen(c).gridPos == 4
            yLbl = chosen(c).y(1) + rangeOffset;   % below (YDir reversed)
        else
            yLbl = chosen(c).y(1) - rangeOffset;   % above
        end
        text(axE, xLbl, yLbl, num2str(chosen(c).gridPos), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'Color', 'k');
    end
    set(axE, 'YDir', 'reverse');
    ylim(axE, rlim);
    shotsBehindMin(axE, shotsInWinMin, 1, [0 0 0]);
    hShotLeg = plot(axE, NaN, NaN, '-k', 'LineWidth', 1.5, 'DisplayName', 'airgun shot');
    legend(axE, hShotLeg, 'Location', 'southeast', 'FontSize', 7);
    yl2 = ylabel(axE, 'range (m)');
    colormap(axE, cmap);
    set(axE, 'CLim', [-cmax cmax]);
    cb = colorbar(axE);
    cb.Label.String = 'hit v_D (cm s^{-1})';
    cb.Label.FontSize = 18;
    set(axE, 'Position', posE);   % colorbar() shrinks axE; restore it
    set(cb, 'Position', posC);
    set(axE, 'FontSize', 8, 'FontName', 'Times New Roman', 'XTickLabel', [], 'Color', [0.6 0.6 0.6], ...
        'XColor', 'k', 'YColor', 'k', 'Layer', 'top');   % Layer=top so the box draws over the track data
    set(yl2, 'Color', 'k');
    xlim(axE, [winStartMin winEndMin]);

    % --- time series panel: hit-level median +/- IQR, 1 s bins ---
    axes(axT); hold(axT, 'on'); box(axT, 'on'); %#ok<LAXES>
    tminLoc = minutes(obsTloc - t0);
    vLoc = obsVloc(:);
    tb = winStartMin : (1/60) : winEndMin;
    tc = nan(numel(tb)-1, 1); mb = tc; q1b = tc; q3b = tc;
    for b = 1:numel(tb)-1
        sel = tminLoc >= tb(b) & tminLoc < tb(b+1);
        tc(b) = (tb(b)+tb(b+1))/2;
        if sum(sel) >= 3
            mb(b) = median(vLoc(sel), 'omitnan');
            q1b(b) = prctile(vLoc(sel), 25);
            q3b(b) = prctile(vLoc(sel), 75);
        end
    end
    okb = ~isnan(mb);
    fill(axT, [tc(okb); flipud(tc(okb))], [q3b(okb); flipud(q1b(okb))], [0.856 0.856 0.964], ...
        'FaceAlpha', 1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(axT, tc(okb), mb(okb), '-', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.5, 'DisplayName', 'median \pm IQR');
    yline(axT, 0, 'k:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
    ylabel(axT, 'hit v_D (cm s^{-1})');
    xlabel(axT, 'time (min)');
    legend(axT, 'Location', 'northeast', 'FontSize', 7);
    set(axT, 'FontSize', 8, 'FontName', 'Times New Roman', 'Color', [1 1 1]);
    xlim(axT, [winStartMin winEndMin]);
    yl0 = ylim(axT);
    yPad = range(yl0) * 0.06;
    ylim(axT, [yl0(1)-yPad, yl0(2)+yPad]);
    shotsBehindMin(axT, shotsInWinMin, 1, [0 0 0]);
    linkaxes([axE axT], 'x');
end


function plot3DTrack(ax, allTracks, trkIdx, cmax, cmap)
% One track's 3D world-frame trajectory (East/North/Depth), coloured by
% hit-level v_Down on the same scale as the echogram (no separate colorbar
% -- the echogram's colorbar already covers it).
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    tr = allTracks(trkIdx);
    pw = tr.posWorld;
    vpw = tr.velPtWorld;
    if isempty(pw) || isempty(vpw) || any(isnan(pw(:))) || any(isnan(vpw(:)))
        axis(ax, 'off');
        text(ax, 0.5, 0.5, sprintf('track %d: no data', trkIdx), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'FontName', 'Times New Roman');
        return
    end
    x = pw(:, 1); y = pw(:, 2); z = pw(:, 3);   % z positive = down
    x = x - x(1); y = y - y(1); z = z - z(1);   % offset from the track's own start, not the world origin
    vHit = vpw(:, 3) * 100;   % cm/s, positive = down
    vFit = tr.velFitWorld;
    vDownFit = NaN;
    if ~isempty(vFit) && ~any(isnan(vFit)), vDownFit = vFit(3) * 100; end

    pad = 0.22;
    xr = range(x); if xr == 0, xr = 1; end
    yr = range(y); if yr == 0, yr = 1; end
    zr = range(z); if zr == 0, zr = 1; end
    xl = [min(x)-pad*xr, max(x)+pad*xr];
    yl = [min(y)-pad*yr, max(y)+pad*yr];
    zl = [min(z)-pad*zr, max(z)+pad*zr];

    shadowCol = [0.6 0.6 0.6];
    plot3(ax, x, y, zl(2)*ones(size(z)), '-', 'Color', shadowCol, 'LineWidth', 1.6);   % floor shadow

    surface(ax, [x x], [y y], [z z], [vHit vHit], 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);
    plot3(ax, x(1), y(1), z(1), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 5);   % start marker
    colormap(ax, cmap);
    set(ax, 'CLim', [-cmax cmax]);
    xlabel(ax, '\DeltaE (m)', 'Color', 'k');
    ylabel(ax, '\DeltaN (m)', 'Color', 'k');
    zlabel(ax, '\DeltaD (m)', 'Color', 'k');
    title(ax, sprintf('v_D = %.2f cm s^{-1}', vDownFit), ...
        'FontSize', 8, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'Color', 'k');
    view(ax, [35.7782 12.1090]);
    axis(ax, 'square');
    set(ax, 'ZDir', 'reverse');
    xlim(ax, xl); ylim(ax, yl); zlim(ax, zl);
    set(ax, 'FontSize', 7, 'FontName', 'Times New Roman', 'Color', [0.78 0.78 0.78], ...
        'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
end


function shotsBehindMin(ax, xsMin, a, col)
% Faint vertical shot markers (in minutes), drawn behind the data.
    if nargin < 4, col = [0 0 0]; end
    yl = ax.YLim;
    h = gobjects(numel(xsMin), 1);
    for s = 1:numel(xsMin)
        h(s) = plot(ax, [xsMin(s) xsMin(s)], yl, '-', 'Color', [col a], ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    uistack(h, 'bottom');
    ax.YLim = yl;
end


function cmap = buildDivColormap(n)
% blue (negative) - white (zero) - red (positive) diverging colormap.
    half = ceil(n/2);
    blue = [0.15 0.35 0.75]; mid = [1 1 1]; red = [0.75 0.15 0.15];
    top = [linspace(blue(1), mid(1), half)', linspace(blue(2), mid(2), half)', linspace(blue(3), mid(3), half)'];
    bot = [linspace(mid(1), red(1), n-half+1)', linspace(mid(2), red(2), n-half+1)', linspace(mid(3), red(3), n-half+1)'];
    cmap = [top; bot(2:end, :)];
end
