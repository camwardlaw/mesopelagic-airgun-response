% Figure 5: pooled Silent (S1+S2) vs Firing (F1+F2) v_D distributions --
% density comparison flanked by ascent/descent tail-exceedance percentage
% panels. Track-level v_D only (velFitWorld)

close all; clear; clc;

here = fileparts(mfilename('fullpath'));
data_dir = fullfile(here, 'data');
addpath(fullfile(here, 'utils'));

%% ===================== load =====================
TRACKS = load(fullfile(data_dir, 'auv', 'tracks.mat'), 'allTracks');
SHOTS = load(fullfile(data_dir, 'airgun', 'shots.mat'), 'shots');

%% ===================== compute =====================
allTracks = TRACKS.allTracks;

% shot times and firing-period boundaries
shotTimes = sort(SHOTS.shots.time);

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
keep = tMid <= tF2_end;
allTracks = allTracks(keep);
tMid = tMid(keep);
nTracks = numel(allTracks);

% assign each track to S1/F1/S2/F2 
periodIdx = zeros(nTracks, 1);
periodIdx(tMid < tF1_start) = 1;   % S1
periodIdx(tMid >= tF1_start & tMid <= tF1_end) = 2;   % F1
periodIdx(tMid > tF1_end & tMid < tF2_start) = 3;   % S2
periodIdx(tMid >= tF2_start & tMid <= tF2_end) = 4;   % F2

% per-track signed downward velocity, linear-fit 
vel = vertcat(allTracks.velFitWorld);
v = vel(:, 3) * 100;   % cm/s, positive = down, negative = up

vS1 = v(periodIdx == 1 & ~isnan(v));
vF1 = v(periodIdx == 2 & ~isnan(v));
vS2 = v(periodIdx == 3 & ~isnan(v));
vF2 = v(periodIdx == 4 & ~isnan(v));

vAll = v(~isnan(v));
fprintf('Pooled across all tracks: IQR = %.2f to %.2f cm/s\n', prctile(vAll, 25), prctile(vAll, 75));
fprintf('Fastest descent: %.2f cm/s, fastest ascent: %.2f cm/s\n', max(vAll), min(vAll));

% pooled Silent (S1+S2) and Firing (F1+F2) distributions
vQuiet = [vS1; vS2];
vFiring = [vF1; vF2];

% ascent/descent tail thresholds
% Tukey boxplot fence of the pooled Silent (S1+S2) distribution 
tukeyK = 1.5;
q1Quiet = prctile(vQuiet, 25);
q3Quiet = prctile(vQuiet, 75);
iqrQuiet = q3Quiet - q1Quiet;
threshDescQuiet = q3Quiet + tukeyK*iqrQuiet;
threshAscQuiet = q1Quiet - tukeyK*iqrQuiet;
fprintf('Silent (S1+S2) Tukey thresholds: descent > %.2f cm/s, ascent < %.2f cm/s\n', ...
    threshDescQuiet, threshAscQuiet);

% shared-bandwidth density estimate 
kdeBW = 0.4;
xg = linspace(-15, 15, 600);
densQuiet = ksdensity(vQuiet, xg, 'Bandwidth', kdeBW);
densFiring = ksdensity(vFiring, xg, 'Bandwidth', kdeBW);
maxYPooled = max([densQuiet, densFiring]) * 1.1;

% --- pooled tail statistics ---
% Tail-exceedance rate, Quiet vs Firing: two-proportion z-test.
pctDescPooled = [100*mean(vQuiet > threshDescQuiet), 100*mean(vFiring > threshDescQuiet)];
pctAscPooled = [100*mean(vQuiet < threshAscQuiet), 100*mean(vFiring < threshAscQuiet)];

x1dP = sum(vQuiet > threshDescQuiet);
n1P = numel(vQuiet);
x2dP = sum(vFiring > threshDescQuiet);
n2P = numel(vFiring);
[~, pDescPooled] = two_prop_test(x1dP, n1P, x2dP, n2P);
x1aP = sum(vQuiet < threshAscQuiet);
x2aP = sum(vFiring < threshAscQuiet);
[~, pAscPooled] = two_prop_test(x1aP, n1P, x2aP, n2P);
fprintf('\nTail exceedance, Silent (S1+S2) vs Firing (F1+F2):\n');
fprintf('  descent: %.1f%% -> %.1f%% (p=%.4g)\n', pctDescPooled(1), pctDescPooled(2), pDescPooled);
fprintf('  ascent: %.1f%% -> %.1f%% (p=%.4g)\n', pctAscPooled(1), pctAscPooled(2), pAscPooled);

% Median v_D of all the flagged (exceeding-threshold) tracks
descAll = [vQuiet(vQuiet > threshDescQuiet); vFiring(vFiring > threshDescQuiet)];
ascAll = [vQuiet(vQuiet < threshAscQuiet); vFiring(vFiring < threshAscQuiet)];
fprintf('\nTail medians (Silent+Firing pooled, tracks exceeding threshold):\n');
fprintf('  descent (v_D > %.2f cm/s): n=%d, median=%.2f cm/s\n', threshDescQuiet, numel(descAll), median(descAll));
fprintf('  ascent (v_D < %.2f cm/s): n=%d, median=%.2f cm/s\n', threshAscQuiet, numel(ascAll), median(ascAll));

groupCols = [0.5 0.5 0.5; 0.80 0.10 0.10];
nameQuiet = 'S1+S2';
nameFiring = 'F1+F2';

%% ===================== plot: pooled Silent (S1+S2) vs Firing (F1+F2) =====================
f5 = figure('Position', [100 100 1450 420]);   % wide/short -- suits the single-row layout
tl5 = tiledlayout(f5, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% shared y-ceiling 
insetYMax = max([pctAscPooled(:); pctDescPooled(:)]) * 1.5;

axAscPct = nexttile(tl5, 1, [1 1]);
draw_tail_pct_single(axAscPct, pctAscPooled, nameQuiet, nameFiring, pAscPooled, ...
    groupCols(1,:), groupCols(2,:), 10, sprintf('up (v_D < %.1f cm s^{-1})', threshAscQuiet), insetYMax);

axPooledDensity = nexttile(tl5, 2, [1 3]);
draw_pdf_pair(axPooledDensity, xg, densQuiet, densFiring, groupCols(1,:), groupCols(2,:), ...
    nameQuiet, nameFiring, median(vQuiet), median(vFiring), maxYPooled, 10, ...
    [threshDescQuiet threshAscQuiet], vQuiet, vFiring, 2*[threshAscQuiet threshDescQuiet]);

axDescPct = nexttile(tl5, 5, [1 1]);
draw_tail_pct_single(axDescPct, pctDescPooled, nameQuiet, nameFiring, pDescPooled, ...
    groupCols(1,:), groupCols(2,:), 10, sprintf('down (v_D > %.1f cm s^{-1})', threshDescQuiet), insetYMax);

forceFont(f5);

%% ===================== helpers =====================

function draw_pdf_pair(ax, xg, dRef, dTest, colRef, colTest, nameRef, nameTest, medRef, medTest, maxY, fontSize, threshLines, rawRef, rawTest, xRange)
% One density-overlay panel (two groups). threshLines is an optional
% vector of x-values (e.g. [threshDesc threshAsc]) drawn as solid gray
% vertical lines marking the tail threshold(s); pass [] to skip.
% rawRef/rawTest are the raw-data vectors for the same two groups being
% plotted -- a thin standard box-and-whisker strip for each is drawn above
% the density curves, colored to match, so the threshold lines can be
% visually checked against the actual box-plot whiskers (the ref group's
% whisker should coincide with the threshold line, by construction). xRange
% (optional, default [-8 8]) sets the x-axis crop.
    if nargin < 16 || isempty(xRange), xRange = [-8 8]; end
    hold(ax, 'on');
    area(ax, xg, dRef, 'FaceColor', colRef, 'FaceAlpha', 0.20, 'EdgeColor', colRef, 'LineWidth', 2, 'DisplayName', nameRef);
    area(ax, xg, dTest, 'FaceColor', colTest, 'FaceAlpha', 0.20, 'EdgeColor', colTest, 'LineWidth', 2, 'DisplayName', nameTest);
    yAtMedRef = interp1(xg, dRef, medRef);
    yAtMedTest = interp1(xg, dTest, medTest);
    plot(ax, [medRef medRef], [0 yAtMedRef], '--', 'Color', colRef, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    plot(ax, [medTest medTest], [0 yAtMedTest], '--', 'Color', colTest, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    for tIdx = 1:numel(threshLines)
        plot(ax, [threshLines(tIdx) threshLines(tIdx)], [0 maxY], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.4, 'HandleVisibility', 'off');
    end
    xlabel(ax, 'v_D (cm s^{-1})');
    ylabel(ax, 'probability density');
    xlim(ax, xRange);
    yHi = maxY;
    % two thin box-and-whisker strips, stacked above the density curves
    boxH = maxY * 0.06;
    gapH = maxY * 0.03;
    yCur = maxY + gapH;
    for grp = 1:2
        if grp == 1, gData = rawRef; gCol = colRef; else, gData = rawTest; gCol = colTest; end
        yBot = yCur;
        yTop = yCur + boxH;
        yMid = (yBot+yTop)/2;
        q1 = prctile(gData, 25);
        q3 = prctile(gData, 75);
        medB = median(gData);
        iqrB = q3 - q1;
        inFence = gData(gData >= q1 - 1.5*iqrB & gData <= q3 + 1.5*iqrB);
        wLo = min(inFence);
        wHi = max(inFence);
        plot(ax, [wLo q1], [yMid yMid], '-', 'Color', gCol, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [q3 wHi], [yMid yMid], '-', 'Color', gCol, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [wLo wLo], [yMid-boxH*0.2 yMid+boxH*0.2], '-', 'Color', gCol, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [wHi wHi], [yMid-boxH*0.2 yMid+boxH*0.2], '-', 'Color', gCol, 'LineWidth', 1, 'HandleVisibility', 'off');
        rectangle(ax, 'Position', [q1, yBot, q3-q1, boxH], 'FaceColor', gCol + (1-gCol)*0.8, 'EdgeColor', gCol, 'LineWidth', 1);
        plot(ax, [medB medB], [yBot yTop], '-', 'Color', gCol, 'LineWidth', 1, 'HandleVisibility', 'off');
        yCur = yTop + gapH;
        yHi = yCur;
    end
    ylim(ax, [0 yHi]);
    legend(ax, 'Location', 'southwest', 'FontSize', fontSize-1, 'Box', 'off');
    set(ax, 'TickDir', 'out', 'FontSize', fontSize, 'XColor', 'k', 'YColor', 'k');
    cleanBox(ax);
    for tIdx = 1:numel(threshLines)
        text(ax, threshLines(tIdx), maxY, sprintf('%.1f', threshLines(tIdx)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', fontSize-1, 'FontWeight', 'bold', 'Color', [0.5 0.5 0.5], 'BackgroundColor', 'w', 'Margin', 1);
    end
end

function draw_tail_pct_single(ax, pct, nameRef, nameTest, pVal, colRef, colTest, fontSize, dirLabel, yMax)
% Grouped-bar tail-exceedance panel: one cluster (Up or Down) in its own
% axes -- % of tracks past a threshold, ref bar + test bar, with a
% significance bracket 
% (Ascent% | Density | Descent%) 
    hold(ax, 'on');
    barWidthP = 0.42;
    xOffsetP = 0.26;
    xC = 1;
    xRef = xC - xOffsetP;
    xTest = xC + xOffsetP;
    fillRef = colRef + (1 - colRef) * 0.8;
    fillTest = colTest + (1 - colTest) * 0.8;
    bar(ax, xRef, pct(1), barWidthP, 'FaceColor', fillRef, 'EdgeColor', colRef, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    bar(ax, xTest, pct(2), barWidthP, 'FaceColor', fillTest, 'EdgeColor', colTest, 'LineWidth', 1.5, 'HandleVisibility', 'off');

    baseLabelPad = max(pct) * 0.03;
    text(ax, xRef, baseLabelPad, nameRef, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize-1, 'FontWeight', 'bold', 'Color', 'k');
    text(ax, xTest, baseLabelPad, nameTest, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize-1, 'FontWeight', 'bold', 'Color', 'k');

    pctLabelPad = max(pct) * 0.03;
    text(ax, xRef, pct(1) + pctLabelPad, sprintf('%.1f%%', pct(1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize-1, 'Color', 'k');
    text(ax, xTest, pct(2) + pctLabelPad, sprintf('%.1f%%', pct(2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize-1, 'Color', 'k');

    tickH = 0.3;
    bracketY = max(pct) + 2;
    draw_sig_bracket(ax, xRef, xTest, bracketY, tickH, pVal, 'k');

    ylim(ax, [0 yMax]);
    xlim(ax, [xC-0.65, xC+0.65]);
    set(ax, 'XTick', xC, 'XTickLabel', {dirLabel}, 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k', 'FontSize', fontSize);
    cleanBox(ax);
    ylabel(ax, '% exceeding v_D threshold');
end

function [z, p] = two_prop_test(x1, n1, x2, n2)
% Two-proportion z-test, pooled SE.
    p1 = x1/n1; p2 = x2/n2;
    pPool = (x1+x2)/(n1+n2);
    sePool = sqrt(pPool*(1-pPool)*(1/n1+1/n2));
    z = (p2-p1)/sePool;
    p = erfc(abs(z)/sqrt(2));   % = 2*(1-normcdf(abs(z))), but stable when z is large
end

function draw_sig_bracket(ax, x1, x2, y, tickH, pval, col)
% Draws a horizontal significance bracket between x1 and x2 with short
% downward ticks at each end, labeled with conventional stars.
    plot(ax, [x1 x1 x2 x2], [y-tickH y y y-tickH], '-', 'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');
    lbl = sig_stars(pval);
    if strcmp(lbl, 'n.s.')
        yOff = tickH*1.5;
        fw = 'bold';
    else
        yOff = tickH*0.5;
        fw = 'normal';
    end
    text(ax, (x1+x2)/2, y + yOff, lbl, 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'FontWeight', fw, 'Color', col, 'HandleVisibility', 'off');
end

function s = sig_stars(p)
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end
