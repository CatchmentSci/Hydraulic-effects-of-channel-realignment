% Initial code to reproduce Figure B1
% Edit the input variables below to match the locations of the data on your
% PC. Data variables described can be accessed/downloaded from:
% the GitHub repo at: https://github.com/CatchmentSci/Hydraulic-effects-of-channel-realignment
% or Zenodo repository at: https://zenodo.org/records/16748995.

% clear workspace
clear all; close all; clc;

root_dir = 'C:\_git_local\Hydraulic-effects-of-channel-realignment\';

ninputs = {'run13.mat', 'run14.mat', 'run15.mat'};

for r = 1:length(ninputs)

    % Preallocate struct array to store results
    Results(r) = struct('RunName', [], ...
        'inchannel_values', [], ...
        'floodplain_values', []);

    runName = sprintf('Run%02d', r+12);

    load([root_dir 'data\' char(ninputs(r))]);

    % Store in struct
    Results(r).RunName = runName;
    Results(r).inchannel_values = inchannel_values;
    Results(r).floodplain_values = floodplain_values;

end

%% check to see if the data has been processed previously, and bring it in if it has, otherwise process
% this takes some time to run

if exist([root_dir 'data\matched_flow_data.mat'], 'file')
    load([root_dir 'data\matched_flow_data.mat']);
else
    %% Bring in the PT data for site 1
    T           = readtable([root_dir 'data\Site 1.xlsx'], 'Sheet','Sheet1');
    dateIn      = T{:,1};
    levelIn     = T{:,2}; % offset of + 0.3418m has already been applied to the data in the spreadsheet
    idx         = ~isnan(levelIn);
    dateUse1    = dateIn(idx);
    levelUse1   = levelIn(idx);
    clear T idx

    %% Bring in the PT data for site 2
    T           = readtable([root_dir 'data\Site 2.xlsx'], 'Sheet','Sheet 2');
    dateIn2     = T{:,1};
    levelIn2    = T{:,2}; % offset of + 0.3827m has already been applied to the data in the spreadsheet
    idx         = ~isnan(levelIn2);
    dateUse2    = dateIn2(idx);
    levelUse2   = levelIn2(idx);

    % match the timestamps
    for a = 1:length(dateUse1)
        dt = cellfun(@(x) x - dateUse1(a), {dateUse2}, 'un', 0); % find where s1 matches s2
        [val_a(a) idx_a(a)] = min(abs(dt{:}));
    end
    use_idx = val_a < '00:03:00';

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

    save([root_dir 'data\matched_flow_data.mat'], 's1_date_matched', 's1_level_matched', 's2_date_matched', 's2_level_matched');
end

%% bring in the rating curve data
T2 = table2array(readtable([root_dir 'data\rating_curve_h_grid_s2_1.csv'],'FileType','delimitedtext')); % The tab where the data is stored is called 'Data'

%% then subplot B)
s2_level_comp = s2_level_matched(42190:42190+(44170-42211));

f0 = figure();
f0.Position = [1000, 596, 731, 642];


% ---- NEW: two stacked subplots ----
ax1 = subplot(2,1,1);  % top: Q vs h
ax2 = subplot(2,1,2);  % bottom: percent vs Q
hold(ax1,'on');
hold(ax2,'on');

n = length(Results);

for a = 1:n

    if ~isempty(Results(a).inchannel_values)

        mod_out = Results(a).inchannel_values;

        % (kept as you wrote it)
        x = sum(mod_out(:,289:end));
        y = s2_level_comp;

        % Logical mask: keep only non-zero x values
        mask = x ~= 0;

        if a == 1
            color = [0.5 0.5 0.5];
        elseif a == 2
            color = [215,25,28]./255;
        else
            color = [44,123,182]./255;
        end

        % ---- TOP subplot: scatter Q vs h ----
        h1{a} = scatter(ax1, x(mask), y(mask), ...
            '+', ...
            'MarkerEdgeColor', color, ...
            'LineWidth', 1.2);

        % ---- BOTTOM subplot: percent within-channel flow ----
        % (kept as you wrote it)
        Ain = sum(Results(a).inchannel_values(:,289:end));
        Afp = sum(Results(a).floodplain_values(:,289:end));
        percent_channel{a} = (Ain ./ (Afp)) .* 100;
        overLim = percent_channel{a} > 100;
        percent_channel{a}(overLim) = 100;

        h3{a} = plot(ax2, Ain, percent_channel{a}, ...
            '-', ...
            'Color', color, ...
            'LineWidth', 1.2);

    end

    max_inchannel(a) = max(sum(Results(a).inchannel_values));

    % perform NSE analysis
    ratingQ_values = ([1:length(y); ...
        interp1(T2(:,2), T2(:,4), y, 'linear')']');
    modQ_values = ([1:length(y); ...
        sum(Results(a).inchannel_values(:,289:end))]');
    [NSout{a}] = nashsutcliffe(ratingQ_values, modQ_values)

end

% ---- TOP subplot: rating curve ----
h2 = plot(ax1, T2(:,4), T2(:,2), '-', 'Color', [0 0 0]);


% ---- AXIS LIMITS / FORMATTING ----
xlim(ax1, [0.7 22])
ylim(ax1, [0.8 1.75])

set(ax1,'fontname','Arial')
set(ax1,'fontweight','normal')
set(ax1,'fontsize',11)
xticks(ax1, [0 5 10 15 20 25])
ax1.XTickLabel = [];
yticks(ax1, [0.8 1 1.2 1.4 1.6])

ylabel(ax1, 'h [m]', 'Interpreter', 'none')

% ---- BOTTOM subplot formatting ----
xlim(ax2, [0.7 22])                 % match x-range
set(ax2,'fontname','Arial')
set(ax2,'fontweight','normal')
set(ax2,'fontsize',11)
xticks(ax2, [0 5 10 15 20 25])

ylim(ax2, [80 100]);
yticks(ax2, (0:5:100));
ylabel(ax2, '% of total Q');
xlabel(ax2, 'In-channel Q [m^{3} s^{-1}]');

% ---- Legend stays on TOP subplot ----
lgd = legend(ax1, [h1{1} h1{2} h1{3} h2], ...
    {'\it n\rm = 0.03', '\it n\rm = 0.035', '\it n\rm = 0.04', 'Rating curve'},...
    'Box', 'off', ...
    'Position',[0.144131493796606 0.75794763842352 0.224999996008617 0.165476185934884]);

lgd.Title.String = 'Channel roughness';


annotation(f0,'textbox',...
    [0.0461436388508892 0.901869158878505 0.0331997264021888 0.0404984423676013],...
    'String',{'A)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);

annotation(f0,'textbox',...
    [0.0461436388508892 0.515576324869538 0.0331997264021888 0.0404984423676011],...
    'String',{'B)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);


% ---- Reduce vertical spacing between subplots ----
pos1 = ax1.Position;
pos2 = ax2.Position;
ax2.Position(2) = ax2.Position(2) + 0.1 ;


exportgraphics(f0, 'calibration_plots.png', 'Resolution', 600)