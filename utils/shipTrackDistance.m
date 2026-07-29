function [t_nav, cum_km_nav] = shipTrackDistance(nav)
% Cumulative along-track distance (km) from ship nav lat/lon/time, used to
% place echograms and AUV/airgun tracks on a common transect-distance axis.
    t_nav = datetime(nav.time(:), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
    ok = ~isnat(t_nav) & isfinite(nav.lat) & isfinite(nav.lon);
    t_nav = t_nav(ok);
    lat = nav.lat(ok);
    lon = nav.lon(ok);
    dkm = sqrt(diff(lat).^2 + (diff(lon).*cosd(lat(1:end-1))).^2) * 111.32;
    cum_km_nav = [0; cumsum(dkm)];
end
