% Initial code to reproduce Figure B2
% Edit the input variables below to match the locations of the data on your
% PC. Data variables described can be accessed/downloaded from:
% the GitHub repo at: https://github.com/CatchmentSci/Hydraulic-effects-of-channel-realignment
% or Zenodo repository at: https://zenodo.org/records/16748995.

% clear workspace
clear all; close all; clc;

%% -------------------- LOAD MODEL RUNS --------------------
root_dir = 'C:\_git_local\Hydraulic-effects-of-channel-realignment\';

ninputs = {'run02', 'run06.mat','run07.mat','run08.mat','run09.mat','run10.mat','run11.mat','run12.mat'};

% Preallocate struct array to store results
Results = repmat(struct('RunName', [], 'inchannel_values', [], 'floodplain_values', []), length(ninputs), 1);

for r = 1:length(ninputs)

    % Name based on file (e.g. run06)
    runLabel = erase(ninputs{r}, ".mat");     % "run06"
    Results(r).RunName = runLabel;

    % Load run
    load(fullfile([root_dir 'data\'], ninputs{r}));

    % Store in struct
    Results(r).inchannel_values   = inchannel_values;
    Results(r).floodplain_values  = floodplain_values;
end

%% -------------------- LOAD / BUILD MATCHED FLOW DATA --------------------

if exist(fullfile(root_dir,'data','matched_flow_data.mat'), 'file')
    load(fullfile(root_dir,'data','matched_flow_data.mat'));
else
    %% Bring in the PT data for site 1
    T           = readtable(fullfile(root_dir,'data','Site 1.xlsx'), 'Sheet','Sheet1');
    dateIn      = T{:,1};
    levelIn     = T{:,2}; % offset already applied in spreadsheet
    idx         = ~isnan(levelIn);
    dateUse1    = dateIn(idx);
    levelUse1   = levelIn(idx);
    clear T idx

    %% Bring in the PT data for site 2
    T           = readtable(fullfile(root_dir,'data','Site 2.xlsx'), 'Sheet','Sheet 2');
    dateIn2     = T{:,1};
    levelIn2    = T{:,2}; % offset already applied in spreadsheet
    idx         = ~isnan(levelIn2);
    dateUse2    = dateIn2(idx);
    levelUse2   = levelIn2(idx);

    % match the timestamps
    val_a = NaT(size(dateUse1));
    idx_a = zeros(size(dateUse1));
    for a = 1:length(dateUse1)
        dt = cellfun(@(x) x - dateUse1(a), {dateUse2}, 'un', 0);
        [val_a(a), idx_a(a)] = min(abs(dt{:}));
    end
    use_idx = val_a < duration(0,3,0);  % within 3 minutes

    % remove questionable data from s1
    rem_num             = datenum(['19/02/2020 00:00'; '25/06/2020 23:00'],'dd/mm/yyyy HH:MM');
    dateUse1_datenum    = datenum(dateUse1);
    rem_idx             = find(dateUse1_datenum > rem_num(1) & dateUse1_datenum < rem_num(2));
    levelUse1(rem_idx)  = NaN;

    % create the matched datasets
    s1_date_matched     = dateUse1(use_idx,1);
    s1_level_matched    = levelUse1(use_idx,1);
    s2_date_matched     = dateUse2(idx_a(use_idx));
    s2_level_matched    = levelUse2(idx_a(use_idx));

    save(fullfile(root_dir,'data','matched_flow_data.mat'), ...
        's1_date_matched','s1_level_matched','s2_date_matched','s2_level_matched');
end

%% -------------------- LOAD RATING CURVE --------------------
T2 = table2array(readtable(fullfile(root_dir,'data','rating_curve_h_grid_s2_1.csv'), ...
    'FileType','delimitedtext'));

%% -------------------- SELECT COMPARISON WINDOW --------------------
s2_level_comp = s2_level_matched(42190:42190+(44170-42211));

%% -------------------- PLOT: SCATTER + RATING CURVE + COLOR RAMP --------------------
f0 = figure();
ax1 = gca;
hold(ax1,'on');

n = length(Results);

% --- coolwarm colormap ---
% NOTE: requires coolwarm.m on your MATLAB path.
cmap = coolwarm(n);

% Plot each run
h1 = cell(n,1);

for a = 1:n
    if ~isempty(Results(a).floodplain_values)

        mod_out = Results(a).floodplain_values;

        x = sum(mod_out(:,289:end));
        y = s2_level_comp;

        mask = x ~= 0;

        color = cmap(a,:);

        h1{a} = scatter(ax1, x(mask), y(mask), ...
            '+', ...
            'MarkerEdgeColor', color, ...
            'LineWidth', 1.2);
    end
end

% Rating curve
h2 = plot(ax1, T2(:,4), T2(:,2), '-', 'Color', [0 0 0]);

%% -------------------- AXIS LIMITS / FORMATTING --------------------
xlim(ax1, [0.7 22])
ylim(ax1, [0.8 1.75])

set(ax1,'fontname','Arial')
set(ax1,'fontweight','normal')
set(ax1,'fontsize',11)
xticks(ax1, [0 5 10 15 20 25])
yticks(ax1, [0.8 1 1.2 1.4 1.6])

ylabel(ax1, 'h [m]', 'Interpreter', 'none')
xlabel(ax1, 'Total Q [m^{3} s^{-1}]');

%% -------------------- COLOR RAMP (COLORBAR) WITH RUN LABELS --------------------
runLabels = {'0.04', '0.05', '0.06', '0.07', '0.08', '0.09', '0.10', '0.15'};

colormap(ax1, cmap);

% FIX: dummy scatter must have n points to match CData length n
cbDummy = scatter(ax1, nan(n,1), nan(n,1), 1, (1:n)', 'filled', ...
    'Visible','off'); %#ok<NASGU>

caxis(ax1, [0.5, n + 0.5]);

cb = colorbar(ax1);
cb.Ticks = 1:n;
cb.TickLabels = runLabels;
cb.Label.String = 'Floodplain Manning''s {\it n}';
cb.FontName = 'Arial';
cb.FontSize = 10;


% ---- Legend stays on TOP subplot ----
lgd = legend(h2, ...
    {'Rating curve'},...
    'Box', 'off', ...
    'Position',[0.144131493796606 0.75794763842352 0.224999996008617 0.165476185934884]);


exportgraphics(f0, 'floodplain_calibration_plots.png', 'Resolution', 600)