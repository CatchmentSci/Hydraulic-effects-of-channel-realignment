%% Pooley vs Goldrill: 5-min interp, synced start, % change, post-peak return time, and volumes
% + Rainfall hyetograph added on right axis (0 at top), x-axis fixed to flows

clearvars; close all; clc

% -------- Input files --------
root_dir = 'C:\_git_local\Hydraulic-effects-of-channel-realignment\data\'; 
fnFlow = 'pooley comp.csv';
fnRain = 'Brotherswater.csv';

%% -------- Read flow table --------
T = readtable(fnFlow, 'VariableNamingRule','preserve');

% Extract columns
tP_raw = T.("dateTime pooley");
qP_raw = T.("Q pooley");

tG_raw = T.("dateTime goldrill");
qG_raw = T.("Q goldrill");

% -------- Convert to datetime --------
tP = datetime(tP_raw, 'InputFormat',"yyyy-MM-dd'T'HH:mm:ss");
tG = datetime(tG_raw, 'InputFormat',"dd/MM/yyyy HH:mm");

% -------- Remove missing rows --------
idxP = ~isnat(tP) & ~isnan(qP_raw);
idxG = ~isnat(tG) & ~isnan(qG_raw);

tP = tP(idxP);  qP = qP_raw(idxP);
tG = tG(idxG);  qG = qG_raw(idxG);

% Sort
[tP, ordP] = sort(tP); qP = qP(ordP);
[tG, ordG] = sort(tG); qG = qG(ordG);

% -------- Build timetables --------
ttP = timetable(tP, qP, 'VariableNames', {'Q_Pooley'});
ttG = timetable(tG, qG, 'VariableNames', {'Q_Goldrill'});

% -------- Interpolate Pooley to regular 5-min grid --------
tGrid = (ttP.Properties.RowTimes(1):minutes(5):ttP.Properties.RowTimes(end))';
ttP_5min = retime(ttP, tGrid, 'linear');

% -------- Sync start times --------
tCommonStart = max(ttP_5min.Properties.RowTimes(1), ttG.Properties.RowTimes(1));

ttP_sync = ttP_5min(ttP_5min.Properties.RowTimes >= tCommonStart, :);
ttG_sync = ttG(ttG.Properties.RowTimes >= tCommonStart, :);

if height(ttP_sync) < 3 || height(ttG_sync) < 2
    error('After syncing start times, one or both series has too few points to proceed.');
end

% -------- % difference in Goldrill: first vs last (after common start) --------
G_first = ttG_sync.Q_Goldrill(1);
G_last  = ttG_sync.Q_Goldrill(end);
pctDiff_G = (G_last - G_first) / G_first * 100;   % percent change

% -------- Target Pooley discharge: same % change relative to Pooley start --------
P_start  = ttP_sync.Q_Pooley(1);
P_target = P_start * (1 + pctDiff_G/100);

% -------- Post-peak only return time in Pooley (blue dashed line) --------
qP_series = ttP_sync.Q_Pooley;
tP_series = ttP_sync.Properties.RowTimes;

% Peak in Pooley
[Qpk, idxPk] = max(qP_series);

% Search only after peak (including peak onwards)
qPost = qP_series(idxPk:end);
tPost = tP_series(idxPk:end);

% Define return criterion on recession:
% If target <= peak, find first time q falls to <= target.
if P_target <= Qpk
    idxRel = find(qPost <= P_target, 1, 'first');
else
    idxRel = [];
end

if isempty(idxRel)
    tReturn = NaT;
    qReturn = NaN;
else
    tReturn = tPost(idxRel);
    qReturn = qPost(idxRel);
end

% -------- Volume calculations --------
% Goldrill volume: from common start to end of plotted Goldrill record
tG_series = ttG_sync.Properties.RowTimes;
qG_series = ttG_sync.Q_Goldrill;

tG_sec = seconds(tG_series - tG_series(1));
volGoldrill_m3 = trapz(tG_sec, qG_series);  % m^3

% Pooley volume: from common start to the post-peak return time (blue dashed line)
if isnat(tReturn)
    volPooley_m3 = NaN;
else
    useP = tP_series <= tReturn;
    tP_sec = seconds(tP_series(useP) - tP_series(find(useP,1,'first')));
    volPooley_m3 = trapz(tP_sec, qP_series(useP)); % m^3
end

% -------- Print results --------
fprintf('\nCommon start time: %s\n', datestr(tCommonStart));

fprintf('\nGoldrill (from common start):\n');
fprintf('  First Q = %.4f m^3/s at %s\n', G_first, datestr(tG_series(1)));
fprintf('  Last  Q = %.4f m^3/s at %s\n', G_last,  datestr(tG_series(end)));
fprintf('  Percent difference (last vs first) = %.2f %%\n', pctDiff_G);
fprintf('  Volume (common start -> end of Goldrill plotted record) = %.3e m^3\n', volGoldrill_m3);

fprintf('\nPooley (from common start):\n');
fprintf('  Start Q = %.4f m^3/s at %s\n', P_start, datestr(tP_series(1)));
fprintf('  Peak  Q = %.4f m^3/s at %s\n', Qpk, datestr(tP_series(idxPk)));
fprintf('  Target Q (same %% change as Goldrill) = %.4f m^3/s\n', P_target);

if isnat(tReturn)
    fprintf('  Post-peak return time: NOT FOUND (so Pooley volume to return cannot be computed)\n');
else
    fprintf('  Post-peak return time: %s (Q = %.4f m^3/s)\n', datestr(tReturn), qReturn);
    fprintf('  Volume (common start -> return time) = %.3e m^3\n', volPooley_m3);
end
fprintf('\n');

volGoldrill_hm3 = volGoldrill_m3 / 1e6;
volPooley_hm3   = volPooley_m3   / 1e6;

fprintf('Goldrill total water volume = %.3f hm^3\n', volGoldrill_hm3);
fprintf('Pooley Bridge total water volume = %.3f hm^3\n', volPooley_hm3);
fprintf('Pooley Bridge total water volume is greater than Goldrill by a factor of %.2f\n', ...
        volPooley_hm3 ./ volGoldrill_hm3);

%% -------- Read rainfall data --------
R = readtable(fnRain, 'VariableNamingRule','preserve');

% ASSUMPTION: datetime in col 1, rainfall (mm) in col 2
% If your file differs, change these two lines to the correct columns/headers.
tR_raw = R{:,1};
p_raw  = R{:,2};

tR = parseDatetimeFlexible(tR_raw);
p  = double(p_raw);

idxR = ~isnat(tR) & ~isnan(p);
tR = tR(idxR);
p  = p(idxR);

[tR, ordR] = sort(tR);
p = p(ordR);

% Keep rainfall points that fall within the FLOW x-axis window (does not change xlim)
xlim_flow = [tCommonStart, max(tP_series(end), tG_series(end))];
useR = (tR >= xlim_flow(1)) & (tR <= xlim_flow(2));
tR_plot = tR(useR);
p_plot  = p(useR);

%% -------- Plot (flows + rainfall) --------
figure('Color','w');
hold on; box on; grid on;

% LEFT AXIS: both flows share discharge axis
yyaxis left
hP = plot(tP_series, qP_series, 'k-', 'LineWidth', 1.2); hold on
hG = plot(tG_series, qG_series, 'r-', 'LineWidth', 1.2);
ylabel('Discharge, Q (m^3 s^{-1})');
ax = gca;
ax.YColor = 'k';

% (optional) show the return time as blue dashed line
if ~isnat(tReturn)
    xline(tReturn, 'b--', 'LineWidth', 1.2);
end

% Fix x-axis to the original flow extent (rainfall cannot change this)
xlim(xlim_flow);

% RIGHT AXIS: rainfall hyetograph, zero at top
yyaxis right
hR = bar(tR_plot, p_plot, 1.0, 'FaceAlpha', 0.25, 'EdgeColor','none');
ylabel('Rainfall (mm)');
ax = gca;
ax.YAxis(2).Direction = 'reverse';   % hyetograph style: 0 at top
ax.YAxis(2).Color = [0 0 1];         % rainfall axis in blue
hR.FaceColor = [0 0 1];

xlabel('Time');
title('Pooley Bridge vs Goldrill Beck discharge with rainfall hyetograph');
xtickangle(45);

legend([hP hG hR], ...
       {'Pooley Bridge Q (5-min interp.)','Goldrill Beck Q','Rainfall'}, ...
       'Location','best');

% --- Export ---
exportgraphics(gcf, 'pooley_comp_with_rain.pdf', 'ContentType','vector');
exportgraphics(gcf, 'pooley_comp_with_rain.png', 'Resolution',600);

%% -------- Helper: flexible datetime parser --------
function t = parseDatetimeFlexible(x)
    if isdatetime(x)
        t = x;
        return
    end

    if isnumeric(x)
        % Try Excel serial first
        try
            t = datetime(x, 'ConvertFrom','excel');
            return
        catch
        end
        % Fallback to datenum-like
        t = datetime(x, 'ConvertFrom','datenum');
        return
    end

    if iscell(x)
        x = string(x);
    end
    x = string(x);

    fmts = [ ...
        "dd/MM/yyyy HH:mm:ss", ...
        "dd/MM/yyyy HH:mm", ...
        "yyyy-MM-dd HH:mm:ss", ...
        "yyyy-MM-dd HH:mm", ...
        "yyyy-MM-dd'T'HH:mm:ss", ...
        "MM/dd/yyyy HH:mm:ss", ...
        "MM/dd/yyyy HH:mm" ...
    ];

    t = NaT(size(x));
    for i = 1:numel(fmts)
        try
            tt = datetime(x, 'InputFormat', fmts(i));
            ok = ~isnat(tt);
            t(ok) = tt(ok);
            if all(~isnat(t))
                return
            end
        catch
        end
    end

    % Last resort
    try
        t = datetime(x);
    catch
        t = NaT(size(x));
    end
end
