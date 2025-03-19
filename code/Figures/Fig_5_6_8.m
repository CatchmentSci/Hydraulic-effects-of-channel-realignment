% Initial code to reproduce Figures 5, 6, 8
% Edit the input variables below to match the locations of the data on your
% PC. Data variables described can be accessed/downloaded from:
% [insert link].

clear all; close all; clc

% directory containing data and scripts downloaded from github/zenodo. 
% should contain subfolders called 'data' and 'Scripts'.
root_dir = 'C:\_git_local\Hydraulic-effects-of-channel-realignment\'; 

% specify the folder where generated data/outputs will be stored to.
data_out = [root_dir 'data\']; 

% ensure matlab can find required m files
addpath(genpath(root_dir)); 


%% check to see if the data has been processed previously, and bring it in if it has, otherwise process
% this takes some time to run

if exist ([root_dir 'data\matched_flow_data.mat'])
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

    save([root_dir 'data\matched_flow_data.mat'], 's1_date_matched', 's1_level_matched', 's2_date_matched', 's2_level_matched', 'dateUse1', 'levelUse1');

end


%% perform basflow seperation
interval    = datenum(s1_date_matched(2)) - datenum(s1_date_matched(1));
input1      = transpose(datenum(s1_date_matched(1)):interval:datenum(s1_date_matched(end)));
idx_use     = find(s1_level_matched>0);
Interp      = interp1(datenum(s1_date_matched(idx_use)),s1_level_matched(idx_use),input1);


filter = 0.995; % filter coefficient for baseflow separation
PKThreshold = 0.03; % peak threshold for runoff events
ReRa = 0.1; % return ratio

% baseflow separation using digital filter method
[stormflow, baseflow] = separatebaseflow([input1, Interp], 0.995, 2);

% extract runoff events from stormflow (baseflow-free hydrograph)
[runoffEvents, nRunoffEvent] = extractrunoff(stormflow, PKThreshold, ReRa, 0.001, 0.0001, 4);

% remove poorly constrained runoff events
nEvents = 1:length(runoffEvents);
rem = [];
rem = [7:8,14:24,55:58,71,78:80,104:105,111, 114,133, 135:138, 141:143, 144, 156:157,160, 188:193, 199 ,202:207, 209:212, 214, 232];
keeper = nEvents(~ismember(nEvents,rem));

for a =1:length(keeper)
    runoffEvents_keep{a,1} = runoffEvents{keeper(a),1};
end

% plot the extracted runoff events on the hydrograph
%figure();
%plot(input1(:,end),Interp(:,end),'k'); hold on;
for a = 1:length(runoffEvents_keep)
    %plot(runoffEvents_keep{a,1}(:,1),runoffEvents_keep{a,1}(:,2)); hold on;
    %text(runoffEvents_keep{a,1}(1,1),runoffEvents_keep{a,1}(1,2),num2str(a))
    event_period(a,1:2) = [runoffEvents_keep{a,1}(1,1),runoffEvents_keep{a,1}(end,1)];

end
event_period = datetime(event_period,'ConvertFrom','datenum');



%% bring in the extended rating curve data for upper site
load(['D:\OneDrive - Newcastle University\Documents - Goldrill Beck Research\General\_shared\Lisflood-FP\site1_flow_input.mat'], 'rating' );



%% undertake some xcorr analysis

clear R maxVar peak_stage idx
b = 1;
aa = 1;
array1 = [];
array2 = [];


for aa = 1:length(event_period)
    t1              = find(min(abs(event_period(aa,1) - s1_date_matched)) == abs(event_period(aa,1) - s1_date_matched));
    idx(aa,1)       = t1(1);
    t2              = find(min(abs(event_period(aa,2) - s1_date_matched)) == abs(event_period(aa,2) - s1_date_matched));
    idx(aa,2)       = t2(1);
    array1{aa,1}    = s2_level_matched(idx(aa,1):idx(aa,2));
    array2{aa,1}    = s1_level_matched(idx(aa,1):idx(aa,2));

    % modification to pull out peak values for xcorr
    peak_idx        = find(max(array2{aa,1}) == array2{aa,1});
    ranger = 60;
    if peak_idx < ranger + 1
        array2{aa,1}    = array2{aa,1}(1:ranger);
        array1{aa,1}    = array1{aa,1}(1:ranger);
    elseif peak_idx + ranger > length(array2{aa,1})
        array2{aa,1}    = array2{aa,1}(peak_idx-ranger:end);
        array1{aa,1}    = array1{aa,1}(peak_idx-ranger:end);
    else
        array2{aa,1}    = array2{aa,1}(peak_idx-ranger:peak_idx+ranger);
        array1{aa,1}    = array1{aa,1}(peak_idx-ranger:peak_idx+ranger);
    end

    % remove start data for s1
    padder1 = repmat(NaN,12,1);
    array2{aa,1} = [array2{aa,1}; padder1];
    array1{aa,1} = [padder1; array1{aa,1}];

    temp            = s1_date_matched(idx(aa,1):idx(aa,2));
    array3(aa,1)    = temp(1); % extract the start time of the event
    temp2(aa,1:2)   = [nanmax(s2_level_matched(idx(aa,1):idx(aa,2))), nanmax(s1_level_matched(idx(aa,1):idx(aa,2)))];
    array4(aa,1)    = temp2(aa,2) - temp2(aa,1); % s1 peak minus s2 peak
    array5{aa,1}    = s1_date_matched(idx(aa,1):idx(aa,2));
    

    for a = 1:60
        shifter = a;
        prepData = zeros(1,shifter)';
        prepData = replace_num(prepData,0,NaN);
        xIn = [array1{aa,1};prepData];
        yIn = [prepData;array2{aa,1}];
        [xIn,yIn] = prepareCurveData(xIn,yIn);
        if ~isempty(xIn)
            R(a,aa) = corr(xIn,yIn);
        else
            R(1:100,aa) = NaN;
        end
        %max_corr(aa,1) = max(R(:,aa))
    end
    if ~isnan(R(1,aa)) && max(R(:,aa)) >0.90 && find ( R(:,aa) == max(R(:,aa))) > 11 % need an r2 of >0.90 to be accepted
        peak_stage(aa,1) = nanmax(array1{aa,1}); % peak flow at site 2 (lower)
        maxVar(aa,1) = find ( R(:,aa) == max(R(:,aa)));
        maxVar(aa,1) = maxVar(aa,1) - length(padder1);
    else
        maxVar(aa,1) = NaN;
        peak_stage(aa,1) = NaN;
    end

end


% convert peak stage to peak flow using the imported rating curve
for a = 1:length(peak_stage)
    if ~isnan(peak_stage(a,1))
        idx(a) = find(min(abs(peak_stage(a,1) - rating(:,1))) == abs(peak_stage(a,1) - rating(:,1))) ;
        peak_q(a,1) = rating(idx(a),2);
    else
        peak_q(a,1) = NaN;
    end
end


%% Testing for changes in transmission time with stage in pre- and post-config
first_after = min(find(datenum(array3)-738338 > 0));
%s1 = scatter(maxVar(1:first_after-1)*5,peak_q(1:first_after-1),[], 'r', 'filled'); hold on;
idx = find(peak_q(first_after) == peak_q(find(~isnan(peak_q))));
idx_use = find(~isnan(peak_q));
idx_use = idx_use(idx:end);
%t1 = scatter(maxVar(idx_use)*5,peak_q(idx_use),[], 'b', 'filled'); 
y = maxVar(idx_use)*5;
X = [peak_q(idx_use), ones(length(peak_q(idx_use)),1)]   ;
[b,bint,r,rint,stats] = regress(y,X); % p = 3rd

%% pull out stats of before vs after median values

nanmedian(maxVar(1:first_after-1))
nanmedian(maxVar(first_after:end))
in1 = maxVar(1:first_after-1);
in1 = in1(~isnan(in1));
in2 = maxVar(first_after:end);
in2 = in2(~isnan(in2));
[p,h,stats] = ranksum(in1, in2); %% MW U test



%% Transmission time plot
f0=figure(1); hold on;
f0.Units='pixels';
%set(f0,'Position',[100, 0, 2480./4, 3508./2]); % A4 aspect ratio
f0.Units='normalized';
ax0 = gca;
set(ax0, 'DefaultAxesFontName', 'Arial'); 
set(ax0,'fontsize', 12);
ylabel('Transmission Time [mins]')


s1_1 = scatter(array3(1:first_after-1),maxVar(1:first_after-1)*5,(2*pi*peak_q(1:first_after-1))*0.5,'filled',...
    'MarkerFaceColor', [215,25,28]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

hold on;

s1_2 = scatter(array3(first_after:end),maxVar(first_after:end)*5,(2*pi*peak_q(first_after:end))*0.5,'filled',...
    'MarkerFaceColor',  [44,123,182]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes


leg_min1 = scatter([array3(20)-20], [75], 2*pi*(min(peak_q))*0.5, 'filled',...
    'MarkerFaceColor',  [215,25,28]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_min2 = scatter([array3(20)+20], [75], 2*pi*(min(peak_q))*0.5, 'filled',...
    'MarkerFaceColor',  [44,123,182]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_25_1 = scatter([array3(20)-22], 80, 2*pi*(prctile(peak_q,25))*0.5, 'filled',...
    'MarkerFaceColor',   [215,25,28]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_25_2 = scatter([array3(20)+22], 80, 2*pi*(prctile(peak_q,25))*0.5, 'filled',...
    'MarkerFaceColor',   [44,123,182]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_50_1 = scatter(array3(20)-25, 85, 2*pi*(prctile(peak_q,50))*0.5, 'filled',...
    'MarkerFaceColor',  [215,25,28]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_50_2 = scatter(array3(20)+25, 85, 2*pi*(prctile(peak_q,50))*0.5, 'filled',...
    'MarkerFaceColor',  [44,123,182]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_75_1 = scatter(array3(20)-27, 90, 2*pi*(prctile(peak_q,75))*0.5, 'filled',...
    'MarkerFaceColor',  [215,25,28]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_75_2 = scatter(array3(20)+27, 90, 2*pi*(prctile(peak_q,75))*0.5, 'filled',...
    'MarkerFaceColor',  [44,123,182]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_max1 = scatter(array3(20)-42, 96, 2*pi*nanmax(peak_q)*0.5, 'filled',... % peak
    'MarkerFaceColor',  [215,25,28]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

leg_max2 = scatter(array3(20)+42, 96, 2*pi*nanmax(peak_q)*0.5, 'filled',...  % peak
    'MarkerFaceColor',  [44,123,182]./255,...
    'MarkerEdgeColor','none',...
    'MarkerFaceAlpha', 0.6); % plot the delay time in minutes

annotation(f0,'textbox',...
    [0.2375 0.915 0.3 0.06],...
    'String',['Q [m^3 s^{-1}]'],...
    'LineStyle','none', ...
    'Interpreter', 'tex');

annotation(f0,'textbox',...
    [0.2375 0.87 0.3 0.06],...
    'String',[sprintf('%.0f', nanmax(peak_q)) ' (max)'],...
    'LineStyle','none', ...
    'Interpreter', 'tex');

annotation(f0,'textbox',...
    [0.2375 0.83 0.3 0.06],...
    'String',[sprintf('%.0f', prctile(peak_q,75))  ' (Q75)'],...
    'LineStyle','none', ...
    'Interpreter', 'tex');

annotation(f0,'textbox',...
    [0.2375 0.79 0.3 0.06],...
    'String',[sprintf('%.0f', prctile(peak_q,50)) ' (Q50)'],...
    'LineStyle','none', ...
    'Interpreter', 'tex');

annotation(f0,'textbox',...
    [0.2375 0.75 0.3 0.06],...
    'String',[sprintf('%.0f', prctile(peak_q,25))  ' (Q25)'],...
    'LineStyle','none', ...
    'Interpreter', 'tex');

annotation(f0,'textbox',...
    [0.2375 0.71 0.3 0.06],...
    'String',[sprintf('%.1f', nanmin(peak_q)) ' (min)'],...
    'LineStyle','none', ...
    'Interpreter', 'tex');

xtickangle(ax0,45);

exportgraphics(f0,['transmission_time.png'],'Resolution',600)

close all;

%% bypass flows plot

f00=figure(1); hold on;
f00.Units='pixels';
%set(f0,'Position',[100, 0, 2480./4, 3508./2]); % A4 aspect ratio
f00.Units='normalized';
ax00 = gca;
set(ax00, 'DefaultAxesFontName', 'Arial'); 
set(ax00,'fontsize', 10);

cd = colormap(ax00, parula); % take your pick (doc colormap)
cd = coolwarm(length(cd));
cd = flip(cd,1);
interpIn = linspace(min(datenum(array3),[],'omitnan'),max(datenum(array3)),length(cd));
cd = interp1(interpIn,cd,datenum(array3)); % map color to velocity values
scatter(-array4,peak_stage,[],cd,'filled'); % plot the peak level change time in minutes

c = colorbar(ax00);
clim(ax00,[min(datenum(array3)),max(datenum(array3))]);
colormap(ax00,cd);
c.Label.String = 'Date [yyyy]'; % Set the label for the colorbar
set(c,'fontname','Arial')
set(c,'fontweight','normal')
set(c,'fontsize',10)
set(c,'Ticks',linspace(min(c.Ticks),max(c.Ticks),6))
datetick(c,'y',    'keepticks','keeplimits')

axis equal
axis tight
xlim([0 0.7]);
xticks(min(xlim):0.2:0.6)
ylim([0.63 2.2]);

plot([0 max(xlim)], [1.86 1.86], '--', ...
    'Color', [0.5 0.5 0.5]);

xlabel({'Difference in peak stage'; '(lower minus upper) [m]' })
ylabel('Event peak stage (lower station) [m]')

exportgraphics(f00,['bypass_flows.png'],'Resolution',600);

close all;


%% undertake hysteresis analysis (Fig 6)
clear h h2
b = 1;
aa = 1;
array1 = {};
array2 = {};


for aa = 1:length(event_period)

    t1              = find(min(abs(event_period(aa,1) - s1_date_matched)) == abs(event_period(aa,1) - s1_date_matched));
    idx(aa,1)       = t1(1);
    t2              = find(min(abs(event_period(aa,2) - s1_date_matched)) == abs(event_period(aa,2) - s1_date_matched));
    idx(aa,2)       = t2(1);

    array1{aa,1}    = s2_level_matched(idx(aa,1):idx(aa,2));
    half_way        = floor(length(array1{aa,1})./2);
    array1_1half    = array1{aa,1}(1:half_way);
    array1_2half    = array1{aa,1}(half_way+1:half_way.*2);
    low_value       = max(min([array1_1half,array1_2half])); % only need the low value for one array
    idx2             = [min(find(array1{aa,1} > low_value == 1)); max(find(array1{aa,1} > low_value == 1))]; % min and max idx for data
    
    array2{aa,1}    = s1_level_matched(idx(aa,1):idx(aa,2));
    array2_1half    = array2{aa,1}(1:half_way);
    array2_2half    = array2{aa,1}(half_way+1:half_way.*2);

    % clip the data so that the miniumum values are present on the rising
    % and falling limb of the hydrograph
    h{aa,1}         = array2{aa,1}(idx2(1):idx2(2));
    h2{aa,1}        = array1{aa,1}(idx2(1):idx2(2));
    array6{aa,1}    = array5{aa,1}(idx2(1):idx2(2));

   
    % normalise
    Hmin{aa,1}      = min(h{aa,1});
    Hmax{aa,1}      = max(h{aa,1});
    h{aa,1}         = (h{aa,1}-min(h{aa,1})) ./ (Hmax{aa,1} - Hmin{aa,1});
    
    % normalise
    Hmin2{aa,1}      = min(h2{aa,1});
    Hmax2{aa,1}      = max(h2{aa,1});
    h2{aa,1}         = (h2{aa,1}-min(h2{aa,1})) ./ (Hmax2{aa,1} - Hmin2{aa,1});


    % calculate the hysteresis
    for b = 1:length(h2{aa,1})
        if b == 1
            p0_1(b,1) = h{aa,1}(b,1) - 0;
            p0_2(b,1) = h2{aa,1}(b,1) - 0;
            p1(b,1) = (h2{aa,1}(b,1).*(p0_1(b,1)./300)) - (h{aa,1}(b,1) .* (p0_2(b,1)./300));

        else
            p0_1(b,1) = h{aa,1}(b,1) - h{aa,1}(b-1,1);
            p0_2(b,1) = h2{aa,1}(b,1) - h2{aa,1}(b-1,1);
            p1(b,1) = h2{aa,1}(b,1).*(p0_1(b,1)./300) - h{aa,1}(b,1) .* (p0_2(b,1)./300);
        end
        p2(b,1) = p1(b,1) .* 300;
        if b == length(h2{aa,1})
            p3(aa) = 0.5.*sum(p2);
        end
    end
    clear p1 p0_1 p0_2 p1 p2

    % alternative (graphical) method
    vertices = [h{aa,1},h2{aa,1}];  % Coordinates of vertices of a rectangle
    polygon_area(aa,1) = calculate_polygon_area(vertices);

end


% only use those events where there is a xcorr >0.90
use_events = find(nanmax(R)'>0.9);
peak_stage_plot = peak_stage(use_events);

% compare the outputs from the two different methods
clear ans
hyster = p3';
%hyster(:,2) = polygon_area(use_events); % graphical method

idx = find(use_events == first_after-1);
idx_init = use_events(1:idx);
idx_after = use_events(idx:end);
temp_hmax = cell2mat(Hmax);

% convert peak stage to peak flow using the imported rating curve
% rerun to include events that are not used in the corr analysis
for a = 1:length(temp_hmax)
    if ~isnan(temp_hmax(a,1))
        idx(a) = find(min(abs(temp_hmax(a,1) - rating(:,1))) == abs(temp_hmax(a,1) - rating(:,1))) ;
        peak_q(a,1) = rating(idx(a),2);
    else
        peak_q(a,1) = NaN;
    end
end


f2=figure(1); hold on;
f2.Units='pixels';
f2.Units='normalized';
f2.Position = [1.2563    0.7039    0.5174    0.1641];
f2_ax1 = subplot(1,2,1);hold on
set(f2_ax1, 'DefaultAxesFontName', 'Arial'); 
set(f2_ax1,'fontsize', 12);
f2_ax1.Color = 'none';

% hysteresis distribution
x = hyster(idx_init,1); % differs due to removed events
numBins = 20;
[counts, edges] = histcounts(x, numBins);
counts = counts./length(x);
binWidth = edges(2)-edges(1);
binCenters = edges(1:end-1) + binWidth/2;
b1 = bar(binCenters, counts, 'BarWidth', 1, 'FaceColor', [215,25,28]./255,'edgecolor','none'); %r
set(b1,'FaceAlpha',0.6)
grid on;
hold on;
nanmedian(x)

y = hyster(idx_after,1);
numBins = 20;
[counts, edges] = histcounts(y, numBins);
counts = counts./length(y);
binWidth = edges(2)-edges(1);
binCenters = edges(1:end-1) + binWidth/2;
b2 = bar(binCenters, counts, 'BarWidth', 1, 'FaceColor', [44,123,182]./255,'edgecolor','none'); %b
set(b2,'FaceAlpha',0.6)
nanmedian(y)

xlabel('\eta_i');
ylabel('Normalised frequency')

leg1 = legend('Pre-diversion \ \ ($\mu$ = -0.06)',...
    'Post-diversion \ ($\mu$ = -0.15)',...
    'Interpreter','latex',...
    'fontsize', 12,...
    'Position',[0.148128311352239 0.86768263535444 0.321218434423899 0.102142855780465]);

annotation(f2,'textbox',...
    [0.0469502948196414 0.911914556932469 0.0510028765160404 0.0690318389933632],...
    'String',{'A)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);

% peak q vs hysteresis plot
f2_ax2 = subplot(1,2,2);hold on
set(f2_ax2, 'DefaultAxesFontName', 'Arial'); 
set(f2_ax2,'fontsize', 12);
f2_ax2.Color = 'none';
xlabel('Q [m^{3} s^{-1}]');
ylabel('\eta_i')

s1 = scatter ([peak_q(idx_init)]', hyster(idx_init,1),'filled',...
    'MarkerFaceColor', [215,25,28]./255); hold on;
s1.MarkerFaceAlpha = 0.6;
s2 = scatter ([peak_q(idx_after)]', hyster(idx_after,1),'filled',...
        'MarkerFaceColor', [44,123,182]./255); hold on;
s2.MarkerFaceAlpha = 0.6;

annotation(f2,'textbox',...
    [0.496580926611078 0.911914556932469 0.0510028765160404 0.0690318389933632],...
    'String',{'B)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);


exportgraphics(f2,['hysteresis_analysis.pdf'],'Resolution',600)

