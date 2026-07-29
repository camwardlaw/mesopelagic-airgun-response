% Figure 1: experimental site map (regional + local bathymetry) and
% shipboard 38 kHz echogram with AUV/airgun depth profiles overlaid.

close all; clear; clc;

here = fileparts(mfilename('fullpath'));
data_dir = fullfile(here, 'data');
addpath(fullfile(here, 'utils'));

%% ===================== load =====================
SV = load(fullfile(data_dir, 'ship', 'ship_sv_38khz.mat'), 'R', 'Sv', 'time');
NAV = load(fullfile(data_dir, 'ship', 'ship_track.mat'), 'nav');
AUV = load(fullfile(data_dir, 'auv', 'auv_track.mat'));
GUN = load(fullfile(data_dir, 'airgun', 'gun_track.mat'), 't_g', 'gun_depth');
shots_dt = load(fullfile(data_dir, 'airgun', 'shots.mat'), 'shots').shots.time(:);

% GEBCO bathymetry, regional (broad East Coast context) and local (the
% survey site itself)
lon_full = ncread(fullfile(data_dir, 'bathymetry_regional.nc'), 'lon');
lat_full = ncread(fullfile(data_dir, 'bathymetry_regional.nc'), 'lat');
Z_full   = ncread(fullfile(data_dir, 'bathymetry_regional.nc'), 'elevation');
bathLocal.lon = ncread(fullfile(data_dir, 'bathymetry_local.nc'), 'lon');
bathLocal.lat = ncread(fullfile(data_dir, 'bathymetry_local.nc'), 'lat');
bathLocal.Z   = ncread(fullfile(data_dir, 'bathymetry_local.nc'), 'elevation');

%% ===================== compute =====================

% timestamps and per-ping arrays
t_sv = datetime(SV.time(:), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
dt_AUV = AUV.time(:);
depth_AUV = AUV.depth(:);
R = SV.R(:);
Sv = SV.Sv;  

% ship cumulative along-track distance (km), used to place the echogram,
% AUV track, and airgun track all on a common transect-distance axis
[t_nav, cum_km_nav] = shipTrackDistance(NAV.nav);
ship_km = interp1(t_nav, cum_km_nav, t_sv, 'linear', NaN);
auv_km = interp1(t_nav, cum_km_nav, dt_AUV, 'linear', NaN);

% crop the ship echogram to the AUV's own survey window
km_lo = interp1(t_nav, cum_km_nav, min(dt_AUV), 'linear', NaN);
km_hi = interp1(t_nav, cum_km_nav, max(dt_AUV), 'linear', NaN);
in_win = ship_km >= km_lo & ship_km <= km_hi;
Sv = Sv(:, in_win);
ship_km = ship_km(in_win) - km_lo;
auv_km = auv_km - km_lo;

% drop bad pings (dropouts read as a sharp deviation from the local
% running median in mean water-column Sv)
y_bottom = 10;
y_top = 750;
depth_idx = R > y_bottom & R < y_top;
mean_Sv = mean(Sv(depth_idx, :), 1, 'omitnan');
bad_pings = abs(mean_Sv - movmedian(mean_Sv, 200)) > 0.5;
keep_sv = ~bad_pings & isfinite(ship_km)';
Sv_clean = Sv(:, keep_sv);
[ship_km_c, order] = sort(ship_km(keep_sv));
Sv_clean = Sv_clean(:, order);

% airgun track (GPS position + per-shot depth), cropped to the AUV's
% survey window and placed on the same transect-distance axis
in_auv = GUN.t_g >= min(dt_AUV) & GUN.t_g <= max(dt_AUV);
t_g = GUN.t_g(in_auv);
gun_depth = GUN.gun_depth(in_auv);
gun_km = interp1(t_nav, cum_km_nav, t_g, 'linear', NaN) - km_lo;
gun_xlim = [min(gun_km, [], 'omitnan') max(gun_km, [], 'omitnan')];

% firing-period membership of each airgun GPS/depth sample (drives the
% red firing-segment overlay on the airgun track below)
dShot = seconds(diff(shots_dt));
breaks = find(dShot > median(dShot)*5);
period_start = shots_dt([1; breaks+1]);
period_end = shots_dt([breaks; numel(shots_dt)]);
gun_firing = false(size(t_g));
for i = 1:numel(period_start)
    gun_firing(t_g >= period_start(i) & t_g <= period_end(i)) = true;
end
gun_depth_firing = gun_depth;
gun_depth_firing(~gun_firing) = NaN;

% survey start/end + two intermediate markers, placed on the echogram at
% the AUV's own depth at those transect distances
xMarkers = [4.48, 27.8, 20.37, 20.7];
yMarkersDepth = interp1(auv_km, depth_AUV, xMarkers, 'linear');

airgunCol = [0.825 0.825 0.825];   

% regional bathymetry: downsample + lightly smooth 
Z_full = fillmissing(Z_full, 'nearest');       % pre-fill NaN
Z_full = fillmissing(Z_full.', 'nearest').';
target_n = 3000;
sf = target_n / max(size(Z_full));
bathRegional.Z = imgaussfilt(imresize(Z_full, sf, 'bicubic'), 1);   % light smoothing 
bathRegional.lon = linspace(lon_full(1), lon_full(end), size(bathRegional.Z, 1)).';
bathRegional.lat = linspace(lat_full(1), lat_full(end), size(bathRegional.Z, 2)).';

% blue-gray ocean / green land colormap
n = 1024;
n_ocean = round(n * 6000/(6000+1500));
n_land = n - n_ocean;
ocean_cm = [linspace(0.35, 0.80, n_ocean)', linspace(0.55, 0.90, n_ocean)', linspace(0.70, 0.97, n_ocean)'];
land_col = [0.40 0.55 0.35];
cmap_bathRegional = [ocean_cm; repmat(land_col, n_land, 1)];
cmap_bathLocal = [linspace(0.35, 0.80, 256)', linspace(0.55, 0.90, 256)', linspace(0.70, 0.97, 256)'];

fontSize = 8;

%% ===================== plot: echogram (top) + regional/local bathymetry insets (bottom) =====================
f = figure('Units', 'centimeters', 'Position', [2 2 24 25.1], 'Color', 'w');
tl = tiledlayout(f, 7, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- echogram ---
axEcho = nexttile(tl, 1, [5 6]);
imagesc(ship_km_c, R, Sv_clean);
hold on;
h_auv = plot(axEcho, auv_km, depth_AUV, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 3, 'DisplayName', 'AUV');
h_gun = plot(axEcho, gun_km, gun_depth, '-', 'Color', airgunCol, 'LineWidth', 3, 'DisplayName', 'airgun');
h_fire = plot(axEcho, gun_km, gun_depth_firing, '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 3, 'DisplayName', 'firing');
plot(axEcho, xMarkers, yMarkersDepth, 's', 'MarkerFaceColor', [0.2 0.2 0.2], ...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerSize', 7, 'LineWidth', 0.5);
set(axEcho, 'YDir', 'reverse');
ylim([y_bottom y_top]);
xlim(gun_xlim);
clim([-75 -40]);
colormap(axEcho, parula);
cb = colorbar(axEcho);
cb.Label.String = 'S_v (dB re 1 m^{-1})';
cb.FontSize = fontSize;
ylabel(axEcho, 'depth (m)');
xlabel(axEcho, 'transect distance (km)');
set(axEcho, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top', 'FontSize', fontSize);
lg = legend(axEcho, [h_auv h_gun h_fire], 'Location', 'southwest', 'Box', 'on', ...
    'FontSize', fontSize, 'TextColor', 'k', 'Color', 'w');
lg.EdgeColor = 'k';

% --- bottom-left: broad East Coast view with Boston pin + rectangle
% marking the local bathymetry panel's extent ---
axBathyRegional = nexttile(tl, 31, [2 3]);
h_im = imagesc(axBathyRegional, bathRegional.lon, bathRegional.lat, bathRegional.Z.');
h_im.Interpolation = 'bilinear';
set(axBathyRegional, 'YDir', 'normal');
hold(axBathyRegional, 'on');
plot(axBathyRegional, [-64 -61.5 -61.5 -64 -64], [37.7 37.7 38.9 38.9 37.7], 'k-', 'LineWidth', 1.5);
text(axBathyRegional, -59.5, 39.4, 'New England Seamounts', 'FontSize', fontSize-1, 'Color', 'k', ...
    'FontAngle', 'italic', 'HorizontalAlignment', 'center', 'Rotation', -30);
plot(axBathyRegional, -71.06, 42.36, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
text(axBathyRegional, -71.06+0.3, 42.36+0.2, 'Boston, MA', 'FontSize', fontSize-1, 'Color', 'k');
xlim(axBathyRegional, [-77 -55]);
ylim(axBathyRegional, [35 45]);
clim(axBathyRegional, [-6000 1500]);
colormap(axBathyRegional, cmap_bathRegional);
xlabel(axBathyRegional, 'longitude (deg)');
ylabel(axBathyRegional, 'latitude (deg)');
daspect(axBathyRegional, [1 cosd(mean(ylim(axBathyRegional))) 1]);
set(axBathyRegional, 'TickDir', 'out', 'Box', 'on', 'FontSize', fontSize);

% --- bottom-right: local bathymetry with AUV track and named seamounts ---
axBathyLocal = nexttile(tl, 34, [2 3]);
imagesc(axBathyLocal, bathLocal.lon, bathLocal.lat, bathLocal.Z.');
set(axBathyLocal, 'YDir', 'normal');
hold(axBathyLocal, 'on');
plot(axBathyLocal, AUV.longitude, AUV.latitude, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
text(axBathyLocal, -62.99, 38.6, 'Atlantis II', 'FontSize', fontSize-1, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'FontAngle', 'italic');
text(axBathyLocal, -62.24, 38.31, 'Gosnold', 'FontSize', fontSize-1, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'FontAngle', 'italic');
xlim(axBathyLocal, [-64 -61.5]);
ylim(axBathyLocal, [37.7 38.9]);
clim(axBathyLocal, [-5000 0]);
colormap(axBathyLocal, cmap_bathLocal);
xlabel(axBathyLocal, 'longitude (deg)');
ylabel(axBathyLocal, 'latitude (deg)');
daspect(axBathyLocal, [1 cosd(mean(ylim(axBathyLocal))) 1]);
yticks(axBathyLocal, [38 38.5]);
set(axBathyLocal, 'TickDir', 'out', 'Box', 'on', 'FontSize', fontSize, 'YAxisLocation', 'right');

forceFont(f);
