% clear the workspace
clear all; close all; clc;

% directory containing data and scripts downloaded from github/zenodo. 
% should contain subfolders called 'data' and 'scripts'.
root_dir = 'C:\_git_local\Hydraulic-effects-of-channel-realignment\'; 

% due to the large file size of model outputs required to produce the plots 
% for Figure 9 (totaling >66 GB), for replication we use files that have already
% been processed to extract the required information for plotting. 
% For original model simulation outputs, please contact the author.
summary = 1; % 1 = summarised inputs, 0 = raw data

% specify the folder where generated data/outputs will be stored to.
data_out = [root_dir 'data\']; 

% create pie chart of biotopes
f0=figure(1); hold on;
f0.Units='pixels';
f0.Units='normalized';
f0_ax1 = subplot(2,1,1);hold on
set(f0_ax1, 'DefaultAxesFontName', 'Arial');
set(f0_ax1,'fontsize', 12);
f0_ax1.Color = 'none';

% pre data
pre_data = [1579 832 994 8192 320];
post_data = [1790 2538 4734 8405 39];
mypalette = validatecolor(["#CA0020" "#F4A582" "#F7F7F7" "#92C5DE" "#0571B0"], "multiple");

% Create pie charts
p1 = pie(f0_ax1,pre_data);
axis equal
axis off

b = 1;
for a = 1:length(p1)
    try
        p1(a).FaceColor = mypalette(b,:);
        p1(a).EdgeColor = mypalette(b,:);
        b = b + 1;
    catch
    end
end

f0_ax2 = subplot(2,1,2);hold on
p2 = pie(f0_ax2,post_data);
axis equal
axis off

b = 1;
for a = 1:length(p2)
    try
        p2(a).FaceColor = mypalette(b,:);
        p2(a).EdgeColor = mypalette(b,:);
        b = b + 1;
    catch
    end
end
exportgraphics(f0,['pie-plots.png'],'Resolution',1000)


% create a plot of inundation areas
if summary == 0 % if using raw model outputs
    root_dir2  = 'X:\Staff\MTP\Goldrill Lisflood Data\local_out\';

    % ensure matlab can find required m files
    addpath(genpath(root_dir2));

    % Get a list of all .wd files in the directory
    fileList = dir(fullfile(root_dir2, '*.wd'));
    idx = contains({fileList.name}, 'Run_pre_b-');
    matchingFiles = fileList(idx);

    % Loop through each file
    for k = 1:length(matchingFiles)
        filePath = fullfile(root_dir2, matchingFiles(k).name);
        fileContent = readmatrix(filePath, "FileType","text" );
        isWetPre(k,1) = sum(fileContent(:) > 0); % in m2
    end

    % Get a list of all .wd files in the directory
    fileList = dir(fullfile(root_dir2, '*.wd'));
    idx = contains({fileList.name}, 'Run_post_b-');
    matchingFiles = fileList(idx);

    % Loop through each file
    for k = 1:length(matchingFiles)
        filePath = fullfile(root_dir2, matchingFiles(k).name);
        fileContent = readmatrix(filePath, "FileType","text" );
        isWetPost(k,1) = sum(fileContent(:) > 0); % in m2
    end


else % bring in summary data
    load ([root_dir 'data\simB_inundation.mat'])

end

% bring in the flow input data file
fileIn      = [root_dir 'data\experimental run inputs1.xlsx'];
opts        = detectImportOptions(fullfile([fileIn]));
ii          = readtable(fileIn,opts);

% write the input hydrograph
idx = [2 6]; % sim B
q_in = table2array(ii(:,idx(1))) + table2array(ii(:,idx(2))); % sim B
q_in = q_in(~isnan(q_in));
q_in = q_in(289:end);


% remove the spin-up period
inun_pre = isWetPre(289:end);
inun_post = isWetPost(289:end);

% only plot increases in q
maxer = 0;
for a = 1:length(inun_pre)
    if q_in(a) > maxer
        idx2(a) = true;
        maxer = q_in(a);
    else
        idx2(a) = false;
    end

end

f1 = figure(); hold on;
pbaspect([1 1 1])
ax1 = gca;
set(ax1,'DefaultTextFontName','Arial')
set(ax1,'fontsize', 12);

axes(ax1); hold on
ax1.Color = 'none';

offset1 = min(inun_post(idx2)) - min(inun_pre(idx2));

h1 = plot(q_in(idx2), inun_pre(idx2),...
    'Color', [215,25,28]./255,...
    'LineWidth', 2,...
    'LineStyle', '-');

h2 = plot(q_in(idx2), inun_post(idx2)-offset1,...
    'Color', [44,123,182]./255,...
    'LineWidth', 2,...
    'LineStyle', '-');


xlabel('Q [m^{3} s^{-1}]');
ylabel( 'Inundation area [m^{2}]', 'Interpreter', 'tex' )
xlim([0 81]);
ylim([27000 409000])

leg1 = legend('Pre-realignment',...
    'Post-realignment',...
    'Interpreter','tex',...
    'fontsize', 12,...
    'Location', 'northwest');

exportgraphics(f1,['inundation_extents.png'],'Resolution',600)

