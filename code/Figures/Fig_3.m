% Initial code to reproduce Figure 3
% Edit the input variables below to match the locations of the data on your
% PC. Data variables described can be accessed/downloaded from:
% [insert link].

clear all; close all; clc

% directory containing data and scripts downloaded from github/zenodo. 
% should contain subfolders called 'Data' and 'Scripts'.
root_dir = 'D:\OneDrive - Newcastle University\Documents - Goldrill Beck Research\General\_shared\Submission Docs\'; 

% specify the folder where generated data/outputs will be stored to.
data_out = [root_dir 'Data\']; 

% ensure matlab can find required m files
addpath(genpath(root_dir)); 


%% Bring in the cross-section data for site 2
T               = readtable([root_dir 'Data\MS50 Export Site2.xlsx'],'Sheet','Corrected XS'); % The tab where the data is stored is called 'Data'
lengthData      = sum(~isnan(table2array(T(:,11))));
distanceIn      = table2array(T(1:lengthData,9)); % This extracts the data from the first collumn i.e. disctance for site 2
realHgtIn       = table2array(T(1:lengthData,11)); % This extracts the data from the third collumn i.e. corrected hgt for site 2
hgtRange        = [min(realHgtIn), max(realHgtIn)];

if hgtRange(1)<0
    hgtRange(1) = 0;
end
interPoints     = hgtRange(1):0.001:hgtRange(2); %1mm incrememnts


% This loop calculates the area for each 1cm increment in river stage from
% 0m up to the maximum value in the survey.
for a = 1:length(interPoints)
    idx         = realHgtIn <= interPoints(a);
    
    if length(find(idx == 1)) < 2
        areaOut(a)      = NaN;
    else
        xi              = distanceIn(idx);
        yi              = realHgtIn(idx);
        invertedY       = -yi+abs(min(-yi));
        areaOut(a)      = trapz(xi,invertedY);
        hyd_radius(a)   = areaOut(a) ./ perimeter(polyshape([xi,invertedY;xi(1),invertedY(1)]));
    end
end


% The output called areaOut is the xs area for each corresponding stage
% measurement which is found in the array 'realHgtIn'
[interPoints,areaOut] = prepareCurveData(interPoints,areaOut);


%% Bring in the gauging data for site 2
% Site 2: Date Q Area
s2_gauging = {'15/11/2018 14:10', 4.792, 5.45, 1.0910; ... % ADCP: q xsa stage
    '16/11/2018 13:45', 2.667 4.27, 0.9392; ... % ADCP
    '13/03/2019 13:10', 7.88, 7.11, 1.2364;... % ADCP
    '18/02/2021 12:30', 4.852, 5.63, 1.0658;... % ADCP
    '19/02/2021 10:45', 6.945, 7.15, 1.1812;... % ADCP
    '25/02/2021 12:05', 7.154, 6.74, 1.2165;... % ADCP
    '26/02/2021 10:20', 3.503, 4.83, 0.9826;... % ADCP
    '12/03/2021 10:45', 4.599, 5.50, 1.0531;... % ADCP
    '13/03/2021 10:45', 6.613, 6.26, 1.1660;... % ADCP

    '10/10/2018 09:45', 0.79, NaN, 0.7513;... % ECM
    '10/06/2021 12:15', 0.36, NaN, 0.6406;... % ECM
    '11/06/2021 11:20', 0.38, NaN, 0.6549;... % ECM
    '08/07/2021 12:20', 0.36, NaN, 0.6363}; % ECM
s2_gauging_times = datetime(s2_gauging(:,1));


%% check to see if the data has been processed previously, and bring it in if it has, otherwise process
% this takes some time to run

if exist ([root_dir 'Data\matched_flow_data.mat'])
    load([root_dir 'Data\matched_flow_data.mat']);
else
    %% Bring in the PT data for site 1
    T           = readtable([root_dir 'Data\Site 1.xlsx'], 'Sheet','Sheet1');
    dateIn      = T{:,1};
    levelIn     = T{:,2}; % offset of + 0.3418m has already been applied to the data in the spreadsheet
    idx         = ~isnan(levelIn);
    dateUse1    = dateIn(idx);
    levelUse1   = levelIn(idx);
    clear T idx


    %% Bring in the PT data for site 2
    T           = readtable([root_dir 'Data\Site 2.xlsx'], 'Sheet','Sheet 2');
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

    save([root_dir 'Data\matched_flow_data.mat'], 's1_date_matched', 's1_level_matched', 's2_date_matched', 's2_level_matched');

end

%% bring in the rating curve data
T2 = table2array(readtable([root_dir 'Data\rating_curve_h_grid_s2_1.csv'],'FileType','delimitedtext')); % The tab where the data is stored is called 'Data'
   


%% generate the plot (subplot A first)

f0          = figure(); hold on
ax0         = subplot(2,2,1);hold on
ax0.Color   = 'none';
set(ax0,'DefaultTextFontName','Arial')
axes(ax0); hold on

firstHalf   = realHgtIn(1:round(length(realHgtIn)./2));
[r(1)]      = findnearest(1.72,firstHalf);
secondHalf  = realHgtIn(round(length(realHgtIn)./2:end));
[r(2)]      = length(firstHalf) + findnearest(1.72,secondHalf) - 1;

h0 = plot(distanceIn,realHgtIn,...
    'k');

h1 = plot(polyshape([distanceIn(r(1):r(2));distanceIn(r(1))] , [realHgtIn(r(1):r(2));realHgtIn(r(1))]) , ...
    'FaceColor', [254./255,232./255,200./255],... % red
    'FaceAlpha',1.0 ...
    ); hold on;

idx_in  = [744 1801];
h2 = plot(polyshape([distanceIn(idx_in(1):idx_in(2));distanceIn(idx_in(1))] , ...
    [realHgtIn(idx_in(1):idx_in(2));realHgtIn(idx_in(1))]) , ...
    'FaceColor', [224./255,236./255,244./255],... % blue
    'FaceAlpha',1.0 ...
    ); hold on;

axis normal
axis tight
set(ax0,'xLim',[0.29 22.17])
set(ax0,'fontname','Arial')
set(ax0,'fontweight','normal')
set(ax0,'fontsize',11)
xlabel( 'Chainage [m]', 'Interpreter', 'none' );
ylabel( 'h [m]', 'Interpreter', 'none' );


%% then subplot B)

ax1 = subplot(2,2,2);hold on
axes(ax1); hold on
set(ax1,'DefaultTextFontName','Arial')
ax1.Color = 'none';

% load the listflood results
load([root_dir 'Data\run13.mat'])

s2_level_comp = s2_level_matched(42190:42190+(44170-42211));

h3_0 = scatter(sum(inchannel_values(:,289:end)),s2_level_comp,...
    'x',...
    'MarkerEdgeColor', [0.5 0.5 0.5],...;
    'MarkerFaceColor', [0.5 0.5 0.5]);

% median est
h3 = plot(T2(:,4),T2(:,2),...
    'Color', [158./255,202./255,225./255],...
    'LineWidth', 2);

% lower bounds
h4 = plot(T2(:,3),T2(:,2),...
    'Color', [158./255,202./255,225./255],...
    'LineWidth', 1,...
    'LineStyle', '--');

% upper bounds
h5 = plot(T2(:,5),T2(:,2),...
    'Color', [158./255,202./255,225./255],...
    'LineWidth', 1,...
    'LineStyle', '--');

% data scatter
h6 = scatter(cell2mat(s2_gauging(:,2)), cell2mat(s2_gauging(:,4)) , ...
    'MarkerEdgeColor', [0.5 0.5 0.5],...
    'MarkerFaceColor', [224./255,236./255,244./255]...
    ); hold on;


axis normal
axis tight
xlim([0 21])
ylim([0.5 1.7])
set(ax1,'fontname','Arial')
set(ax1,'fontweight','normal')
set(ax1,'fontsize',11)
xlabel('Q [m^{3} s^{-1}]');
ylabel( 'h [m]', 'Interpreter', 'none' )


%% finally subplot C)
ax2 = subplot(2,2,3:4);hold on
axes(ax2); hold on
set(ax2,'DefaultTextFontName','Arial')
ax2.Color = 'none';

yyaxis left
elapsed_hr = 5.*(1:length(s2_level_comp))./60;

h7 = plot(elapsed_hr, sum(floodplain_values(:,289:end)),...
    'Color', [158./255,202./255,225./255],...
    'LineWidth', 2,...
    'LineStyle', '-');
axis normal
axis tight
set(ax2,'fontname','Arial')
set(ax2,'fontweight','normal')
set(ax2,'fontsize',11)
ylabel('Q [m^{3} s^{-1}]');


yyaxis right
h8 = plot(elapsed_hr, s2_level_comp,...
    'Color', [252./255,146./255,114./255],...
    'LineWidth', 2,...
    'LineStyle', '-');

ylim([0.5 2])
ylabel( 'h [m]', 'Interpreter', 'none' )
xlabel( 'Elapsed time [hr]', 'Interpreter', 'none' )

%text(74,1.7,'\delta t = 5-min ','Interpreter','tex') 

annotation(f0,'textbox',...
    [0.0617142857142855 0.907142857142861 0.0615 0.0476190476190478],...
    'String',{'A)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);

annotation(f0,'textbox',...
    [0.517071428571427 0.88571428706249 0.0678571413670268 0.0690476176994188],...
    'String',{'B)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);

annotation(f0,'textbox',...
    [0.0420714285714283 0.388095239443448 0.0678571413670268 0.0690476176994188],...
    'String',{'C)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);


exportgraphics(f0,['calibration_workflow.png'],'Resolution',600)





