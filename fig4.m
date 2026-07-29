% Figure 4: shot-locked and population-level vertical velocity response
% across exposure periods + statistical tests

close all; clear; clc;

here = fileparts(mfilename('fullpath'));
data_dir = fullfile(here, 'data');
addpath(fullfile(here, 'utils'));

%% ===================== load =====================
TRACKS = load(fullfile(data_dir, 'auv', 'tracks.mat'), 'allTracks');
SHOTS = load(fullfile(data_dir, 'airgun', 'shots.mat'), 'shots');

%% ===================== compute =====================
allTracks = TRACKS.allTracks;
shotTimes = sort(SHOTS.shots.time);

% firing periods, restricted to the track-data time window
allHitTimes = [];
for i = 1:length(allTracks)
    allHitTimes = [allHitTimes; vertcat(allTracks(i).posLog{:, 5})]; %#ok<AGROW>
end
tData_start = min(allHitTimes);
tData_end = max(allHitTimes);
shotTimes = shotTimes(shotTimes >= tData_start & shotTimes <= tData_end);
gap = seconds(diff(shotTimes));
splitIdx = find(gap > 300, 1);
tF1_start = shotTimes(1);
tF1_end = shotTimes(splitIdx);
tF2_start = shotTimes(splitIdx+1);
tF2_end = shotTimes(end);

% per-track midpoint time
nTracks = length(allTracks);
tMid = NaT(nTracks, 1, 'TimeZone', 'UTC');
for i = 1:nTracks
    times = vertcat(allTracks(i).posLog{:, 5});
    tMid(i) = times(round(end/2));
end

t0 = tData_start;

periodNames = {'S1', 'F1', 'S2', 'F2'};   % S = silent, F = firing
firingRed = [0.8 0.1 0.1];
periodColors = [0.5 0.5 0.5; firingRed; 0.5 0.5 0.5; firingRed];

% ===== shot-level analysis =====
% hit-level v_Down (raw, per-hit, not track-aggregated) pooled across all
% tracks. Shots (~20 s apart) are the unit of replication: well beyond a
% track's ~5-20 s lifespan, so shot windows don't overlap and are independent.
obsT = NaT(0, 1, 'TimeZone', 'UTC');
obsDown = [];
for i = 1:nTracks
    tt = vertcat(allTracks(i).posLog{:, 5});
    vpw = allTracks(i).velPtWorld;
    vd = vpw(:, 3);
    hh = allTracks(i).headingHit;
    if isempty(vd) || isempty(hh) || any(isnan(vd(:))) || numel(tt) ~= numel(vd)
        continue
    end
    obsT = [obsT; tt(:)]; %#ok<AGROW>
    obsDown = [obsDown; vd(:) * 100]; %#ok<AGROW>   % cm/s, positive = down
end

tsec = seconds(obsT - t0);
tsec = tsec(:);
downV = obsDown(:);
ssec = seconds(shotTimes - t0);
ssec = ssec(:);
f1s = seconds(tF1_start - t0);
f1e = seconds(tF1_end - t0);
f2s = seconds(tF2_start - t0);
f2e = seconds(tF2_end - t0);

isi = diff(ssec);
isi = isi(isi < 300);   % drop the F1->F2 gap
isiMed = median(isi);
maxLag = isiMed / 2;   % peri-shot profile half-window (s): half the median inter-shot
                        % interval, so neighboring shots' windows don't overlap
lagBin = 0.5;      % profile bin width (s)
nearRadius = 2;    % "near the shot" half-window (s) for the per-shot peak-above-floor stat below
fprintf('Shots: %d (%d in F1, %d in F2)\n', numel(ssec), ...
    sum(ssec>=f1s & ssec<=f1e), sum(ssec>=f2s & ssec<=f2e));
fprintf('Windows: maxLag=%.1f s, lagBin=%.1f s, nearRadius=%.1f s\n', maxLag, lagBin, nearRadius);

% period windows
perWinShot = {[0 f1s], [f1s f1e], [f1e f2s], [f2s f2e]};
perIsFireShot = [false, true, false, true];

% per-shot peak-above-floor (shape statistic): one value per shot, median
% hits near the shot minus median hits in the surrounding flank.
shotPeakVals = cell(1, 4);
for k = 1:4
    refsK = periodRefs(ssec, perWinShot{k}(1), perWinShot{k}(2), perIsFireShot(k), maxLag, isiMed);
    shotPeakVals{k} = perShotPeakAboveFloor(tsec, downV, refsK, maxLag, nearRadius);
end

% ===== period-level analysis =====
% net-displacement/track (linear OLS slope, velFitWorld): did the
% population's overall level shift across the whole dive?
vel = vertcat(allTracks.velFitWorld);
v = vel(:, 3) * 100;   % cm/s, positive = down

tMid_min = minutes(tMid - t0);
f1a_start = minutes(tF1_start - t0);
f1a_end = minutes(tF1_end - t0);
f2a_start = minutes(tF2_start - t0);
f2a_end = minutes(tF2_end - t0);

periodIdx = zeros(nTracks, 1);
periodIdx(tMid < tF1_start) = 1;   % S1
periodIdx(tMid >= tF1_start & tMid <= tF1_end) = 2;   % F1
periodIdx(tMid > tF1_end & tMid < tF2_start) = 3;   % S2
periodIdx(tMid >= tF2_start & tMid <= tF2_end) = 4;   % F2

% 30 s time bins
binWidth_min = 30/60;
binEdges_min = 0 : binWidth_min : max(tMid_min) + binWidth_min;
nBins = length(binEdges_min) - 1;
binCenters = binEdges_min(1:end-1) + binWidth_min/2;

velBinMed = NaN(nBins, 1);
velBinQ1 = NaN(nBins, 1);
velBinQ3 = NaN(nBins, 1);
for b = 1:nBins
    mask = tMid_min >= binEdges_min(b) & tMid_min < binEdges_min(b+1);
    n_b = sum(mask);
    x_b = v(mask);
    if n_b >= 2
        velBinMed(b) = median(x_b, 'omitnan');
        velBinQ1(b) = prctile(x_b, 25);
        velBinQ3(b) = prctile(x_b, 75);
    elseif n_b == 1
        velBinMed(b) = x_b;
        velBinQ1(b) = x_b;
        velBinQ3(b) = x_b;
    end
end

binPerAll = zeros(nBins, 1);
binPerAll(binCenters' < f1a_start) = 1;
binPerAll(binCenters' >= f1a_start & binCenters' <= f1a_end) = 2;
binPerAll(binCenters' > f1a_end & binCenters' < f2a_start) = 3;
binPerAll(binCenters' >= f2a_start & binCenters' <= f2a_end) = 4;

tsV = ~isnan(velBinMed);
bcT = binCenters(tsV)';
vMed = velBinMed(tsV);
vQ1 = velBinQ1(tsV);
vQ3 = velBinQ3(tsV);
tLim = [min(bcT) max(bcT)];

statV = binPerAll > 0;
sPer = binPerAll(statV);
sMed = velBinMed(statV);
uG = unique(sPer);

contrastPairs = [1 2; 2 3; 3 4; 1 3; 2 4];
nContrasts = size(contrastPairs, 1);

fprintf('\nBins: %d total, %d with data, %d in stat units\n', ...
    nBins, sum(tsV), sum(statV));
for pp = 1:4
    fprintf('  %-8s : %d bins\n', periodNames{pp}, sum(sPer==pp));
end

medCol = [0.2 0.2 0.8];

%% ===================== plot: 2 rows x 5 tiles =====================
fig = figure('Units', 'centimeters', 'Position', [2 2 32 18], 'Color', 'w');
tl = tiledlayout(fig, 2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% shot-level (peri-shot) profiles/stat box on top row, period-level
% (whole-dive evolution) box on bottom row
profTile0 = 1;
statTile = 5;
evoTile0 = 6;
evoBoxTile = 10;

% --- period-level: evolution time series (spanned, 4 tiles) + stat box (1 tile) ---
axEvo = nexttile(tl, evoTile0, [1 4]); hold(axEvo, 'on'); box(axEvo, 'on');
ylimEvo = [-3.5 6];

fill(axEvo, [bcT; flipud(bcT)], [vQ3; flipud(vQ1)], medCol, ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(axEvo, bcT, vMed, '-', 'Color', medCol, 'LineWidth', 1.5, 'DisplayName', 'median \pm IQR');
plot(axEvo, tLim, [0 0], 'k:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
xlim(axEvo, tLim);
ylim(axEvo, ylimEvo);
xlabel(axEvo, 'time (min)');
ylabel(axEvo, 'v_D (cm s^{-1})');
legend(axEvo, 'Location', 'southeast', 'FontSize', 7);
set(axEvo, 'FontSize', 8, 'FontName', 'Times New Roman');
drawFiringBar(axEvo, f1a_start, f1a_end, f2a_start, f2a_end, tLim, periodColors);

axBox = nexttile(tl, evoBoxTile); hold(axBox, 'on'); box(axBox, 'on');
plot(axBox, [0 5], [0 0], 'k:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
drawPeriodBox(axBox, sMed, sPer, uG, periodColors, periodNames);
ylabel(axBox, 'per-bin median v_D (cm s^{-1})');
ylim(axBox, ylimEvo);   % match the evolution panel's ylim so all axes align
addNlabels(axBox, sPer, uG);
[pk1, tk1] = kruskalwallis(sMed, sPer, 'off');
fprintf('\n[period-level] per-bin median vD KW: H(%d)=%.2f, p=%.4g\n', tk1{2,3}, tk1{2,5}, pk1);
by1 = bracketY(axBox, sMed, contrastPairs);
printBoxContrasts(contrastPairs, by1, uG, sMed, sPer, periodNames, nContrasts, axBox);
set(axBox, 'FontSize', 8, 'FontName', 'Times New Roman');

% --- shot-level: peri-shot profile panels (4 tiles) + stat box (1 tile), top row ---
axProf1 = plotProfileRow(tl, profTile0, tsec, downV, ssec, perWinShot, periodNames, perIsFireShot, periodColors, maxLag, lagBin, isiMed, ...
    true, 'v_D (cm s^{-1})');
set(axProf1, 'YLim', ylimEvo);

% per-shot peak-above-floor v_Down, compared across periods with the same
% KW omnibus + pairwise Bonferroni Mann-Whitney machinery as the
% period-level box.
shotVals = [];
shotGrp = [];
for k = 1:4
    shotVals = [shotVals; shotPeakVals{k}]; %#ok<AGROW>
    shotGrp = [shotGrp; k*ones(numel(shotPeakVals{k}), 1)]; %#ok<AGROW>
end
shotUG = unique(shotGrp);

axStat1 = nexttile(tl, statTile); hold(axStat1, 'on'); box(axStat1, 'on');
plot(axStat1, [0 5], [0 0], 'k:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
drawPeriodBox(axStat1, shotVals, shotGrp, shotUG, periodColors, periodNames);
ylabel(axStat1, 'per-shot peak - floor v_D (cm s^{-1})');
ylim(axStat1, ylimEvo);
addNlabels(axStat1, shotGrp, shotUG);
[pk0, tk0] = kruskalwallis(shotVals, shotGrp, 'off');
fprintf('\n[shot-level] per-shot peak-above-floor vD KW across periods: H(%d)=%.2f, p=%.4g\n', tk0{2,3}, tk0{2,5}, pk0);
by0 = bracketY(axStat1, shotVals, contrastPairs, [0.2 0.65 1.15]);   % both n.s.
printBoxContrasts(contrastPairs, by0, shotUG, shotVals, shotGrp, periodNames, nContrasts, axStat1);
set(axStat1, 'FontSize', 8, 'FontName', 'Times New Roman');

forceFont(fig);

allAxes = [axEvo; axBox; axProf1(:); axStat1];
for i = 1:numel(allAxes)
    cleanBox(allAxes(i));
end

%% ===================== helpers =====================

function [ax, P] = plotProfileRow(tl, tileStart, tsec, v, ssec, perWin, perName, perIsFire, periodColors, maxLag, lagBin, isiMed, ...
        signed, ylab)
% N-panel peri-shot profile of v vs. lag to nearest shot, one panel per
% period window. Also returns P (per-period pooled profile), used for
% display only -- the stat box uses per-shot data instead, to avoid
% treating non-independent lag bins as independent replicates.
    nP = numel(perWin);
    legIdx = find(strcmp(perName, 'F2'), 1);   % legend goes on this panel

    P = cell(nP, 1);
    for p = 1:nP
        a = perWin{p}(1);
        b = perWin{p}(2);
        mh = tsec >= a & tsec <= b;
        if perIsFire(p)
            refs = ssec(ssec >= a & ssec <= b);
        else
            if b - a < 2*maxLag + isiMed, continue; end
            refs = (a + maxLag : isiMed : b - maxLag)';   % pseudo 20 s grid
        end
        if numel(refs) < 3 || sum(mh) < 20, continue; end
        [c, m, q1, q3] = periProfile(tsec(mh), v(mh), refs, maxLag, lagBin);
        P{p} = struct('c', c, 'm', m, 'q1', q1, 'q3', q3);
    end

    ax = gobjects(nP, 1);
    for p = 1:nP
        ax(p) = nexttile(tl, tileStart + p - 1); hold(ax(p), 'on'); box(ax(p), 'on');
        if isempty(P{p})
            addCornerLabel(ax(p), perName{p});
            continue
        end
        col = periodColors(p, :);
        ok = ~isnan(P{p}.m);
        fill(ax(p), [P{p}.c(ok) fliplr(P{p}.c(ok))], [P{p}.q3(ok) fliplr(P{p}.q1(ok))], col, ...
            'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax(p), P{p}.c(ok), P{p}.m(ok), '-', 'Color', col, 'LineWidth', 1.8, ...
            'DisplayName', 'median \pm IQR');
        xline(ax(p), 0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        if signed, yline(ax(p), 0, 'k:', 'LineWidth', 0.6, 'HandleVisibility', 'off'); end
        addCornerLabel(ax(p), perName{p});   % top-right corner, instead of a title
        if p == 1, ylabel(ax(p), sprintf('hit %s', ylab)); end
        xlabel(ax(p), 'time (s) relative to shot');
        xlim(ax(p), [-maxLag maxLag]);
        set(ax(p), 'FontSize', 8, 'FontName', 'Times New Roman');
        if ~isempty(legIdx) && p == legIdx, legend(ax(p), 'Location', 'southeast', 'FontSize', 7); end
    end
    yls = cell2mat(get(ax, 'YLim'));
    set(ax, 'YLim', [min(yls(:,1)) max(yls(:,2))]);
end


function addCornerLabel(ax, name)
% period label in the top-right corner of the panel, instead of a title
    text(ax, 0.95, 0.98, name, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'Color', 'k');
end


function [ctrs, med, q1, q3] = periProfile(tsec, v, refsec, maxLag, lagBin)
% peri-event median of v vs lag-to-nearest-reference, with IQR
    lag = nearestLagFast(tsec, refsec);
    edges = -maxLag : lagBin : maxLag;
    ctrs = edges(1:end-1) + lagBin/2;
    med = nan(1, numel(ctrs));
    q1 = med; q3 = med;
    for b = 1:numel(ctrs)
        sel = lag >= edges(b) & lag < edges(b+1);
        n_b = sum(sel);
        if n_b >= 2
            med(b) = median(v(sel), 'omitnan');
            q1(b) = prctile(v(sel), 25);
            q3(b) = prctile(v(sel), 75);
        end
    end
end


function lag = nearestLagFast(tsec, ssec)
% signed lag to nearest shot (>0 = after), O(n log m).
    tsec = tsec(:);
    ss = unique(ssec(:));
    if numel(ss) == 1
        lag = tsec - ss(1);
        return
    end
    iSorted = interp1(ss, (1:numel(ss))', tsec, 'nearest', 'extrap');
    lag = tsec - ss(iSorted);
end


function refs = periodRefs(ssec, a, b, isFire, maxLag, isiMed)
% reference times (real shots or, in quiet periods, a pseudo 20 s grid) for
% one period window [a b].
    if isFire
        refs = ssec(ssec >= a & ssec <= b);
    else
        if b - a < 2*maxLag + isiMed
            refs = [];
        else
            refs = (a + maxLag : isiMed : b - maxLag)';
        end
    end
end


function vals = perShotPeakAboveFloor(tsec, v, refs, maxLag, nearRadius)
% one value per shot: median(hits within +/-nearRadius, the "peak") minus
% median(hits in the (nearRadius, maxLag] flank, the "floor").
    if isempty(refs)
        vals = [];
        return
    end
    vals = nan(numel(refs), 1);
    for i = 1:numel(refs)
        lag = tsec - refs(i);
        inWin = abs(lag) <= maxLag & ~isnan(v);
        if ~any(inWin), continue; end
        lagW = lag(inWin);
        vW = v(inWin);
        nearMask = abs(lagW) <= nearRadius;
        farMask = ~nearMask;   % already restricted to |lag| <= maxLag above
        vals(i) = median(vW(nearMask), 'omitnan') - median(vW(farMask), 'omitnan');
    end
    vals = vals(~isnan(vals));
end


function drawPeriodBox(ax, vals, grp, uG, periodColors, periodNames)
% box plot of bin-level VALS grouped by period, coloured per period
    axes(ax);
    boxplot(vals, grp, 'Labels', periodNames(uG), 'Widths', 0.7, ...
        'Symbol', '+', 'MedianStyle', 'line', 'Notch', 'off');
    hb = findobj(ax, 'Tag', 'Box');
    nB = numel(hb);
    for gi = 1:numel(uG)
        pp = uG(gi);
        patch(get(hb(nB-gi+1), 'XData'), get(hb(nB-gi+1), 'YData'), ...
            periodColors(pp, :)*0.4 + 0.6, 'FaceAlpha', 1);
    end
    medL = findobj(ax, 'Tag', 'Median');
    for mi = 1:numel(medL)
        plot(ax, get(medL(mi), 'XData'), get(medL(mi), 'YData'), '-k', 'LineWidth', 0.5);
    end
    outL = findobj(ax, 'Tag', 'Outliers');
    set(outL, 'Color', 'k', 'MarkerEdgeColor', 'k');
    for gi = 1:numel(uG)
        plot(ax, gi, mean(vals(grp==uG(gi)), 'omitnan'), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
    end
    xlim(ax, [0.4 numel(uG)+0.6]);
end


function addNlabels(ax, grp, uG)
% n = #bins per period, placed near the bottom of the axis
    yl = get(ax, 'YLim');
    for gi = 1:numel(uG)
        text(ax, gi, yl(1) + 0.03*(yl(2)-yl(1)), sprintf('n=%d', sum(grp==uG(gi))), ...
            'HorizontalAlignment', 'center', 'FontSize', 6, 'FontName', 'Times New Roman');
    end
end


function by = bracketY(ax, vals, contrastPairs, fracOverride)
% significance-bracket heights
    if nargin < 4, fracOverride = []; end
    n = size(contrastPairs, 1);
    tier = zeros(n, 1);
    placed = {};
    for i = 1:n
        g1 = min(contrastPairs(i, :));
        g2 = max(contrastPairs(i, :));
        t = 1;
        while true
            if t > numel(placed)
                placed{t} = [g1 g2]; %#ok<AGROW>
                tier(i) = t;
                break
            end
            ivs = placed{t};
            if ~any(g1 < ivs(:,2) & g2 > ivs(:,1))
                placed{t} = [ivs; g1 g2];
                tier(i) = t;
                break
            end
            t = t + 1;
        end
    end
    nTiers = max(tier);

    yl = get(ax, 'YLim');
    dataMax = max(vals(~isnan(vals)));
    room = yl(2) - dataMax;
    if room <= 0, room = 0.1*(yl(2)-yl(1)); end   % fallback if data reaches the axis top

    blockFrac = 0.55;   % how much of the room the tiers span
    blockHeight = blockFrac*room;
    midY = dataMax + 0.4*room;   % tier block centered here, closer to the axis top than mid-room
    lowestY = midY - blockHeight/2;
    if ~isempty(fracOverride)
        frac = fracOverride;
    elseif nTiers == 3
        frac = [0.2 0.55 0.89];   % pulls the two n.s. tiers closer together, adjacent/*** tier lower
    else
        frac = (0:nTiers-1) / max(nTiers-1, 1);
    end
    by = lowestY + reshape(frac(tier), size(tier)) * blockHeight;
end


function drawFiringBar(ax, t1s, t1e, t2s, t2e, tlim, periodColors)
% period bar above the evolution plot 
    yl = get(ax, 'YLim');
    barTop = yl(2);
    barBot = yl(2) - 0.08*(yl(2)-yl(1));
    barMid = (barTop + barBot) / 2;
    segEdges = [tlim(1) t1s t1e t2s t2e tlim(2)];
    segNames = {'S1', 'F1', 'S2', 'F2'};
    for k = 1:4
        a = segEdges(k); b = segEdges(k+1);
        fill(ax, [a b b a], [barBot barBot barTop barTop], periodColors(k, :), ...
            'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        text(ax, mean([a b]), barMid, segNames{k}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'FontName', 'Times New Roman', 'Color', 'k', 'FontWeight', 'bold');
    end
end


function printBoxContrasts(pairs, bracketYv, uGroups, v, grpIdx, pNames, nK, ax)
% pairwise Mann-Whitney on bin-level values (Bonferroni x nK)
    yl = get(ax, 'YLim');
    yRng = yl(2)-yl(1);
    tickH = 0.012*yRng;
    textOffStar = -0.01*yRng;   % asterisk glyphs sit high above baseline even with
                                % VerticalAlignment='bottom', so they need a negative offset
    textOffNS = 0;
    for sp_i = 1:size(pairs, 1)
        g1 = pairs(sp_i, 1);
        g2 = pairs(sp_i, 2);
        x1 = find(uGroups==g1);
        x2 = find(uGroups==g2);
        if isempty(x1) || isempty(x2), continue; end
        v1 = v(grpIdx==g1); v1 = v1(~isnan(v1));
        v2 = v(grpIdx==g2); v2 = v2(~isnan(v2));
        if numel(v1) < 1 || numel(v2) < 1, continue; end
        [p_rs, ~, stats] = ranksum(v1, v2);
        p_adj = min(p_rs*nK, 1);
        n1 = length(v1);
        n2 = length(v2);
        s_pool = sqrt(((n1-1)*var(v1)+(n2-1)*var(v2))/max(n1+n2-2, 1));
        cohens_d = 0;
        if s_pool > 0, cohens_d = (mean(v1)-mean(v2))/s_pool; end
        U = stats.ranksum - n1*(n1+1)/2;
        r_rb = 1 - 2*U/(n1*n2);
        fprintf('  %s vs %s: W=%.0f, p_adj=%.4g, r_rb=%.3f, d=%.3f (med %.2f vs %.2f, n=%d/%d)\n', ...
            pNames{g1}, pNames{g2}, stats.ranksum, p_adj, r_rb, cohens_d, median(v1), median(v2), n1, n2);
        if p_adj<0.001, sigStr='***'; elseif p_adj<0.01, sigStr='**'; elseif p_adj<0.05, sigStr='*'; else, sigStr='n.s.'; end
        textOff = textOffNS;
        if ~strcmp(sigStr, 'n.s.'), textOff = textOffStar; end
        yPos = bracketYv(sp_i);
        xIn = 0.05;
        plot(ax, [x1+xIn x2-xIn], [yPos yPos], '-k', 'LineWidth', 0.8, 'HandleVisibility', 'off', 'Clipping', 'off');
        plot(ax, [x1+xIn x1+xIn], [yPos-tickH yPos], '-k', 'LineWidth', 0.8, 'HandleVisibility', 'off', 'Clipping', 'off');
        plot(ax, [x2-xIn x2-xIn], [yPos-tickH yPos], '-k', 'LineWidth', 0.8, 'HandleVisibility', 'off', 'Clipping', 'off');
        text(ax, mean([x1 x2]), yPos+textOff, sigStr, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 7, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'Clipping', 'off');
    end
end
