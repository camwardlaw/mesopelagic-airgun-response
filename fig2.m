% Figure 2: AUV 38 kHz Sv shifted to absolute depth (AUV depth + range),
% with paired 70 kHz Sp / airgun-waveform / spectrogram insets.

close all; clear; clc;

here = fileparts(mfilename('fullpath'));
data_dir = fullfile(here, 'data');
addpath(fullfile(here, 'utils'));

%% ===================== load =====================
SV = load(fullfile(data_dir, 'auv', 'auv_sv_38khz.mat'), 'R', 'Sv', 'time');
AUV = load(fullfile(data_dir, 'auv', 'auv_track.mat'));
NAV = load(fullfile(data_dir, 'ship', 'ship_track.mat'), 'nav');
SHOTS = load(fullfile(data_dir, 'airgun', 'shots.mat'));
TRACKS = load(fullfile(data_dir, 'auv', 'tracks.mat'), 'allTracks');

% Sp70 insets: pre-computed Sp arrays for each inset's time window
echoL = load(fullfile(data_dir, 'auv', 'sp_161243.mat'));
echoR = load(fullfile(data_dir, 'auv', 'sp_194843.mat'));

% each inset gets an airgun shot recorded during its own window 
V_pk = 3.0;
sensitivity_dB = -177;
airgun_L = extractAirgunShot(fullfile(data_dir, 'airgun', 'shot_161243.wav'), V_pk, sensitivity_dB);
airgun_R = extractAirgunShot(fullfile(data_dir, 'airgun', 'shot_195027.wav'), V_pk, sensitivity_dB);

%% ===================== compute =====================
Sv = SV.Sv; 
R = SV.R(:);
dt_AUV = AUV.time(:);
depth_AUV = AUV.depth(:);
shots_dt = SHOTS.shots.time(:);

% ship along-track distance (km used to place the AUV/Sv onto a 
% transect-distance axis
[t_nav, cum_km_nav] = shipTrackDistance(NAV.nav);
km_lo = interp1(t_nav, cum_km_nav, min(dt_AUV), 'linear', NaN);

% per-ping AUV distance/depth, and the same on the Sv ping timebase
t_sv = datetime(SV.time(:), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
km_AUV = interp1(t_nav, cum_km_nav, dt_AUV, 'linear', NaN) - km_lo;
km_sv = interp1(t_nav, cum_km_nav, t_sv, 'linear', NaN) - km_lo;
dep_sv = interp1(dt_AUV, depth_AUV, t_sv, 'linear', NaN);

% decimate pings for plotting and drop NaNs
dec = 10;
idx = 1:dec:numel(t_sv);
Sv = Sv(:, idx);
km_sv = km_sv(idx);
dep_sv = dep_sv(idx);
ok = ~isnan(km_sv) & ~isnan(dep_sv);
Sv = Sv(:, ok);
km_sv = km_sv(ok);
dep_sv = dep_sv(ok);

% regrid each ping onto a shared absolute-depth axis so imagesc can be used
z_step = median(diff(R));
z_grid = (floor(min(dep_sv) + min(R)) : z_step : ceil(max(dep_sv) + max(R))).';
Nz = numel(z_grid);
Np = numel(km_sv);
Sv_z = nan(Nz, Np);
for k = 1:Np
    if isnan(dep_sv(k)), continue; end
    Sv_z(:, k) = interp1(R + dep_sv(k), Sv(:, k), z_grid, 'linear', NaN);
end
Sv_z = fillmissing(Sv_z, 'linear', 2, 'EndValues', 'none', 'MaxGap', 5);   

% resample onto a uniform km grid so imagesc's uniform-spacing assumption is correct
km_u = linspace(min(km_sv, [], 'omitnan'), max(km_sv, [], 'omitnan'), numel(km_sv));
Sv_u = nan(Nz, numel(km_u));
for row = 1:Nz
    Sv_u(row, :) = interp1(km_sv, Sv_z(row, :), km_u, 'nearest', NaN);
end

% Sp70 insets -- Sp_L/r_L/t_L and Sp_R/r_R/t_R are each inset's window
Sp_L = echoL.Sp;
r_L = echoL.range;
t_L = echoL.time;
t_L.TimeZone = 'UTC';   

% crop to just the T195027 portion (this file spans a wider window)
t_R_full = echoR.time;
t_R_full.TimeZone = 'UTC';
inR = t_R_full >= datetime(2024, 9, 26, 19, 50, 27, 'TimeZone', 'UTC');
Sp_R = echoR.Sp(:, inR);
r_R = echoR.range;
t_R = t_R_full(inR);

% distance (m) at each inset ping, referenced to each file's first ping
km_L = interp1(dt_AUV, km_AUV, t_L, 'linear');
km_R = interp1(dt_AUV, km_AUV, t_R, 'linear');
m_L = (km_L - km_L(1)) * 1000;
m_R = (km_R - km_R(1)) * 1000;

% AUV depth at each inset ping (for the rectangle drawn on the main echogram)
dep_L = interp1(dt_AUV, depth_AUV, t_L, 'linear');
dep_R = interp1(dt_AUV, depth_AUV, t_R, 'linear');

r_lo = 5;
r_hi = 40;

% SPL summary 
spl_dB = SHOTS.shots.SPLpp_dB_re_1uPa;
fprintf('SPLpp (dB re 1 uPa): mean = %.2f, SD = %.2f, min = %.2f, max = %.2f, n = %d\n', ...
    mean(spl_dB, 'omitnan'), std(spl_dB, 'omitnan'), min(spl_dB, [], 'omitnan'), max(spl_dB, [], 'omitnan'), sum(~isnan(spl_dB)));

xAll = [4.48, 27.8, 20.37, 20.7];   % survey start / end + 2 intermediate markers

% all firing bursts across the full shots.mat record
dShotAll = seconds(diff(shots_dt));
breaksAll = find(dShotAll > median(dShotAll)*5);
period_start = shots_dt([1; breaksAll+1]);
period_end = shots_dt([breaksAll; numel(shots_dt)]);

% S1/F1/S2/F2 boundaries, restricted to 300m analyzed window 
allHitTimes = [];
for i = 1:numel(TRACKS.allTracks)
    allHitTimes = [allHitTimes; vertcat(TRACKS.allTracks(i).posLog{:, 5})]; %#ok<AGROW>
end
tData_start = min(allHitTimes);
tData_end = max(allHitTimes);

shotTimesA = shots_dt(shots_dt >= tData_start & shots_dt <= tData_end);
gapA = seconds(diff(shotTimesA));
splitIdxA = find(gapA > 300, 1);
tF1_start = shotTimesA(1);
tF1_end = shotTimesA(splitIdxA);
tF2_start = shotTimesA(splitIdxA+1);
tF2_end = shotTimesA(end);

% silent gaps between/around firing bursts, for the rest of the top bar
allSilent_start = [dt_AUV(1); period_end];
allSilent_end = [period_start; dt_AUV(end)];

silentGray = [0.5 0.5 0.5];
firingRed = [0.8 0.1 0.1];
barAlpha = 0.20;
fontSize = 8;

%% ===================== plot =====================
f = figure('Units', 'centimeters', 'Position', [2 2 24 48], 'Color', 'w');
tl = tiledlayout(f, 258, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- main echogram ---
ax = nexttile(tl, 1, [106 3]);
imagesc(ax, km_u, z_grid, Sv_u, 'AlphaData', ~isnan(Sv_u));
set(ax, 'YDir', 'reverse');
clim(ax, [-90 -40]);
colormap(ax, parula);
cb = colorbar(ax);
cb.Label.String = 'S_v (dB re 1 m^{-1})';
cb.FontSize = fontSize;
xlabel(ax, 'transect distance (km)');
ylabel(ax, 'depth (m)');
set(ax, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top', 'FontSize', fontSize);
ylim([225 750]);
yticks(ax, 200:100:700);
hold(ax, 'on');

plot(ax, km_AUV, depth_AUV, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);   % AUV dive line

yAll_d = interp1(km_AUV, depth_AUV, xAll, 'linear');
plot(ax, xAll, yAll_d, 's', 'MarkerFaceColor', [0.2 0.2 0.2], ...
    'MarkerEdgeColor', [0.9 0.9 0.9], 'MarkerSize', 6, 'LineWidth', 0.5);

% thick firing/silent bar across the top of the echogram 
barTop = 225;
barBot = 260;

% background: ALL firing bursts + ALL silent gaps
hBar = gobjects(0, 1);
for i = 1:numel(period_start)
    x1 = interp1(t_nav, cum_km_nav, period_start(i), 'linear', NaN) - km_lo;
    x2 = interp1(t_nav, cum_km_nav, period_end(i), 'linear', NaN) - km_lo;
    if isnan(x1) || isnan(x2), continue; end
    hBar(end+1) = patch(ax, [x1 x2 x2 x1], [barTop barTop barBot barBot], ...
        firingRed, 'FaceAlpha', barAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off'); %#ok<AGROW>
end
for i = 1:numel(allSilent_start)
    x1 = interp1(t_nav, cum_km_nav, allSilent_start(i), 'linear', NaN) - km_lo;
    x2 = interp1(t_nav, cum_km_nav, allSilent_end(i), 'linear', NaN) - km_lo;
    if isnan(x1) || isnan(x2), continue; end
    hBar(end+1) = patch(ax, [x1 x2 x2 x1], [barTop barTop barBot barBot], ...
        silentGray, 'FaceAlpha', barAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off'); %#ok<AGROW>
end

% foreground: S1/F1/S2/F2 only 
labeledStart = [tData_start; tF1_start; tF1_end; tF2_start];
labeledEnd = [tF1_start; tF1_end; tF2_start; tF2_end];
labeledColor = [silentGray; firingRed; silentGray; firingRed];
hBarLabeled = gobjects(0, 1);
for i = 1:4
    x1 = interp1(t_nav, cum_km_nav, labeledStart(i), 'linear', NaN) - km_lo;
    x2 = interp1(t_nav, cum_km_nav, labeledEnd(i), 'linear', NaN) - km_lo;
    if isnan(x1) || isnan(x2), continue; end
    hBarLabeled(end+1) = patch(ax, [x1 x2 x2 x1], [barTop barTop barBot barBot], ...
        labeledColor(i, :), 'FaceAlpha', barAlpha, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off'); %#ok<AGROW>
end

uistack(hBarLabeled, 'bottom');
uistack(hBar, 'bottom');

labelPeriods(ax, (barTop+barBot)/2, t_nav, cum_km_nav, km_lo, ...
    tData_start, tF1_start, tF1_end, tF2_start, tF2_end, fontSize, 'middle');

% rectangles on the main echogram showing each inset file's window
for fi = 1:2
    if fi == 1
        kk = km_L; dd = dep_L;
    else
        kk = km_R; dd = dep_R;
    end
    x1 = min(kk, [], 'omitnan');
    x2 = max(kk, [], 'omitnan');
    y1 = min(dd, [], 'omitnan') + r_lo;
    y2 = max(dd, [], 'omitnan') + r_hi;
    rectangle(ax, 'Position', [x1, y1, x2-x1, y2-y1], 'EdgeColor', 'k', 'LineWidth', 1);
end

% --- L/R columns, each [inset / waveform / spectrogram] ---
botBlock = tiledlayout(tl, 1, 2, 'TileSpacing', 'tight', 'Padding', 'none');
botBlock.Layout.Tile = 322;   % row 108, col 1 of 3 -> (108-1)*3+1
botBlock.Layout.TileSpan = [151 3];

insetAudioColumn(botBlock, 1, m_L, r_L, Sp_L, [-95 -55], shots_dt, t_L, airgun_L, fontSize, false);
insetAudioColumn(botBlock, 2, m_R, r_R, Sp_R, [-95 -55], shots_dt, t_R, airgun_R, fontSize, true);

forceFont(f);

%% ===================== helpers =====================

function shot = extractAirgunShot(wavFile, V_pk, sensitivity_dB)
% Extract the first clear shot from an airgun wav recording and its
% spectrogram 
    [p_raw, fs] = audioread(wavFile);
    sensitivity_linear = 10^(sensitivity_dB/20);
    p_Pa = (p_raw * V_pk) / sensitivity_linear * 1e-6;
    thr = 0.3 * max(abs(p_Pa));
    first_peak_idx = find(abs(p_Pa) > thr, 1, 'first');
    pre_s = round(0.1 * fs);    % 100 ms before the peak
    post_s = round(1.0 * fs);   % 1000 ms after 
    s_lo = max(1, first_peak_idx - pre_s);
    s_hi = min(numel(p_Pa), first_peak_idx + post_s);
    shot.p = p_Pa(s_lo:s_hi);
    shot.t_ms = ((0:numel(shot.p)-1) - (first_peak_idx - s_lo)) / fs * 1000;

    nfft_spec = 512;
    nfft_zeropad = 2048;
    noverlap_spec = round(nfft_spec * 0.99);
    [S, F_spec, T_spec] = spectrogram(shot.p, hann(nfft_spec), noverlap_spec, nfft_zeropad, fs);
    shot.S_dB = 10*log10(abs(S).^2 + eps) + 120;
    shot.F_spec = F_spec;
    shot.T_spec_ms = T_spec * 1000 - (first_peak_idx - s_lo)/fs*1000;
end

function insetAudioColumn(parentTl, parentTile, m_file, r_file, Sp_file, spClim, ...
    shots_dt, t_file, shot, fontSize, showCbar)
% Sp70 inset on top, that same window's airgun waveform in the middle, 
% its spectrogram at the bottom. 
    col = tiledlayout(parentTl, 20, 1, 'TileSpacing', 'compact', 'Padding', 'none');
    col.Layout.Tile = parentTile;

    axIn = nexttile(col, 1, [12 1]);
    imagesc(axIn, m_file, r_file, Sp_file);
    set(axIn, 'YDir', 'reverse', 'Box', 'off', 'TickDir', 'out', ...
        'FontSize', fontSize-1, 'Color', [0.85 0.85 0.85], ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 0.8);
    ylim(axIn, [10 40]);
    yticks(axIn, 10:5:40);
    xlim(axIn, [0 100]);
    xticks(axIn, 0:20:100);
    clim(axIn, spClim);
    colormap(axIn, parula);
    ylabel(axIn, 'range (m)');
    xlabel(axIn, 'distance (m)');
    hold(axIn, 'on');
    shotRug(axIn, shots_dt, t_file, m_file);
    drawTopRight(axIn);
    if showCbar
        cbIn = colorbar(axIn);
        cbIn.Label.String = 'S_p (dB re 1 m^2)';
        cbIn.FontSize = fontSize-1;
        hShotLegR = plot(axIn, NaN, NaN, '-', 'Color', [0.8 0.1 0.1], ...
            'LineWidth', 1.5, 'DisplayName', 'airgun shot');
        legend(axIn, hShotLegR, 'Location', 'southeast', 'FontSize', fontSize-1);
    end

    % waveform + spectrogram 
    wsBlock = tiledlayout(col, 2, 1, 'TileSpacing', 'tight', 'Padding', 'none');
    wsBlock.Layout.Tile = 13;
    wsBlock.Layout.TileSpan = [8 1];

    t_lim = [-100 450];   
    axW = nexttile(wsBlock, 1, [1 1]);
    plot(axW, shot.t_ms, shot.p / 1000, 'k-', 'LineWidth', 0.8);
    set(axW, 'TickDir', 'out', 'Box', 'off', 'FontSize', fontSize, 'XTickLabel', []);
    xlim(axW, t_lim);
    ylim(axW, [-1 1]);
    yticks(axW, -1:0.5:1);
    ylabel(axW, 'pressure (kPa)');

    axA = nexttile(wsBlock, 2, [1 1]);
    imagesc(axA, shot.T_spec_ms, shot.F_spec / 1000, shot.S_dB);
    set(axA, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'off', 'FontSize', fontSize);
    ylim(axA, [0 1]);
    yticks(axA, 0:0.25:1);
    xlim(axA, t_lim);
    clim(axA, [160 210]);
    colormap(axA, parula);
    xlabel(axA, 'time (ms)');
    ylabel(axA, 'frequency (kHz)');
    if showCbar
        cbA = colorbar(axA);
        cbA.Label.String = 'dB re 1 \muPa^2/Hz';
        cbA.FontSize = fontSize;
        cbA.Ticks = 160:20:200;
        cbA.TickDirection = 'out';
        cbA.TickLength = 0.04;
    end
end

function labelPeriods(ax, yPos, t_nav, cum_km_nav, km_lo, ...
    tData_start, tF1_start, tF1_end, tF2_start, tF2_end, fontSize, vAlign)
% Centered "S1"/"F1"/"S2"/"F2" text labels over actual analyzed window
    names = {'S1', 'F1', 'S2', 'F2'};
    t1 = [tData_start; tF1_start; tF1_end; tF2_start];
    t2 = [tF1_start; tF1_end; tF2_start; tF2_end];
    for i = 1:4
        x1 = interp1(t_nav, cum_km_nav, t1(i), 'linear', NaN) - km_lo;
        x2 = interp1(t_nav, cum_km_nav, t2(i), 'linear', NaN) - km_lo;
        if isnan(x1) || isnan(x2), continue; end
        text(ax, (x1+x2)/2, yPos, names{i}, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', vAlign, 'FontSize', fontSize, ...
            'FontWeight', 'bold', 'Color', 'k', 'HandleVisibility', 'off');
    end
end

function drawTopRight(ax)
% draw black top + right edges so the box looks closed without mirror ticks
    xl = xlim(ax);
    yl = ylim(ax);
    plot(ax, xl, [yl(1) yl(1)], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(ax, xl, [yl(2) yl(2)], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(ax, [xl(2) xl(2)], yl, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

function shotRug(ax, shots_dt, t_file, m_file)
% draw short red tick marks at the top of an Sp panel where airgun shots fired.
    good = ~isnat(t_file);
    t_file = t_file(good);
    m_file = m_file(good);
    if isempty(t_file), return; end
    in_win = shots_dt >= min(t_file) & shots_dt <= max(t_file);
    if ~any(in_win), return; end
    m_shots = interp1(t_file, m_file, shots_dt(in_win), 'linear');
    yl = ylim(ax);
    y_top = yl(1) + 0.015 * diff(yl);
    y_bot = yl(1) + 0.08 * diff(yl);
    for k = 1:numel(m_shots)
        plot(ax, [m_shots(k) m_shots(k)], [y_top y_bot], '-', ...
            'Color', [0.8 0.1 0.1], 'LineWidth', 2, 'HandleVisibility', 'off');
    end
end
