% Initial code to reproduce Figure 2
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


%% Bring in the cross-section data for site 1
T           = readtable([root_dir 'data\MS50 Export Site1.xlsx'],'Sheet','Sheet 1'); % The tab where the data is stored is called 'data'
lengthData  = sum(~isnan(table2array(T(:,11))));
distanceIn  = table2array(T(1:lengthData,11)); % This extracts the distance for site 1
realHgtIn   = table2array(T(1:lengthData,13)); % This extracts the corrected hgt for site 1
hgtRange    = [min(realHgtIn), max(realHgtIn)];

if hgtRange(1)<0
    hgtRange(1) = 0;
end
interPoints  = hgtRange(1):0.001:hgtRange(2); %1mm incrememnts


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


%% Bring in the pressure transducer data for site 1
% Site 1: Date Q Area
s1_gauging = {'15/11/2018 12:40', 4.485 4.93; ... % ADCP
    '16/11/2018 15:10', 2.218 3.97; ... % ADCP
    '13/03/2019 11:50', 6.661, 5.63;... % ADCP
    '18/02/2021 10:25', 4.322, 5.70;... % ADCP - no PT
    '19/02/2021 09:40', 5.303, 5.93;... % ADCP - no PT
    '25/02/2021 10:30', 7.812, 7.41;... % ADCP - no PT
    '26/02/2021 09:30', 3.319, 5.51;... % ADCP - no PT
    '12/03/2021 09:50', 4.003, 5.58;... % ADCP - no PT
    '13/03/2021 09:25', 6.137, 6.58;... % ADCP - no PT

    '10/10/2018 08:30', 0.74, NaN;... % ECM
    '01/11/2018 10:45', 1.03, NaN;... % ECM
    '10/06/2021 10:30', 0.31, NaN;... % ECM
    '11/06/2021 10:15', 0.33, NaN;... % ECM
    '08/07/2021 10:20', 0.31, NaN;... % ECM
    '09/07/2021 13:45', 0.24, NaN}; % ECM
s1_gauging_times = datetime(s1_gauging(:,1));


% Bring in the PT data for site 1
T           = readtable([root_dir 'data\Site 1.xlsx'], 'Sheet','Sheet1');
dateIn      = T{:,1};
levelIn     = T{:,2}; % offset of + 0.3418m has already been applied to the data 

idx         = ~isnan(levelIn);
dateUse     = dateIn(idx);
levelUse    = levelIn(idx);

for a = 1:length(s1_gauging_times)
    dt = cellfun(@(x) x - s1_gauging_times(a), {dateUse}, 'un', 0);
    [val_b(a) idx_b(a)] = min(abs(dt{:}));
end
use_rating_idx  = val_b < '00:15:00';
gauging_level   = levelUse(idx_b(use_rating_idx));
gauging_q       = cell2mat(s1_gauging(use_rating_idx,2));
missing_level   = find(use_rating_idx == 0);


%% use site 3 to compensate for no site 1 data for some gaugings

T               = readtable([root_dir 'data\Site 3.xlsx'], 'Sheet','Sheet1');
dateIn_s3       = T{:,1};
levelIn_s3      = T{:,2};
idx             = ~isnan(levelIn_s3);
dateUse_s3      = dateIn_s3(idx);
levelUse_s3     = levelIn_s3(idx);

dateUse_s1      = dateUse(1:70000); % start thru to 13th May 2019 (before new PT installed at S3)
levelUse_s1     = levelUse(1:70000); 

for a = 1:length(dateUse_s1) 
    dt = cellfun(@(x) x - dateUse_s1(a), {dateUse_s3}, 'un', 0);
    [val_c(a) idx_c(a)] = min(abs(dt{:}));
end

use_s3_idx = val_c < '00:05:00';
%scatter(levelUse_s3(idx_c(use_s3_idx)), levelUse_s1(use_s3_idx,1));
[f, gof] = fit(levelUse_s3(idx_c(use_s3_idx)), levelUse_s1(use_s3_idx,1),'poly5');

for a = missing_level
    dt = cellfun(@(x) x - s1_gauging_times(a), {dateUse_s3}, 'un', 0);
    [val_d(a) idx_d(a)] = min(abs(dt{:}));
end
use_rating_idx = val_d < '00:15:00';
use_rating_idx(1:3) = 0;
levelUse_s3(idx_d(use_rating_idx))

eval_level = feval(f,levelUse_s3(idx_d(use_rating_idx))); % evaluate the fit for the missing river levels
infilled_gauging_level = [gauging_level;eval_level];
infilled_gauging_q = [gauging_q; cell2mat(s1_gauging(missing_level,2))];

%scatter(infilled_gauging_level,infilled_gauging_q)

%% calculate the xsa for each gauging point
for a = 1:length(infilled_gauging_level)

    [r(a)] = findnearest(infilled_gauging_level(a),interPoints);
    infilled_gauging_xsa(a,1) = areaOut(r(a));
    infilled_gauging_1d_v(a,1) = infilled_gauging_q(a)./infilled_gauging_xsa(a,1);
    gauging_hyd_radius(a,1) = hyd_radius(r(a));

    %mannings_n(a,1) = ((gauging_hyd_radius(a).^(2./3)) .* (0.012.^0.5)) ./ infilled_gauging_1d_v(a,1);  % using A/P
    mannings_n(a,1) = ((infilled_gauging_level(a).^(2./3)) .* (0.012.^0.5)) ./ infilled_gauging_1d_v(a,1); 

end

%% calulcate the Manning's n for each flow gauge measurment - need to get the water surface slope

[xData, yData_1]    = prepareCurveData( infilled_gauging_level, mannings_n );

% Set up fittype and options.
ft                  = fittype( 'power1' );
opts                = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display        = 'Off';
opts.StartPoint     = [0.0397081202885953 -2.17306416814885];

% Fit model to data.
[fitresult, gof]    = fit( xData, yData_1, ft, opts );
constant_n          = fitresult(nanmax(xData)); % get the mannins n corresponding to the max gauged flow

% perform a sanity check on the calcs
vel_check           = ((infilled_gauging_level.^(2/3)) .* (0.012.^(0.5))) ./ mannings_n; % this should match the 1d calc obs;
eval_idx            = find(levelUse > max(infilled_gauging_level));

% import the rating curve developed in R Studio
T                   = readtable([root_dir 'data\site1_rating_curve.csv']);

for a = 1:length(levelUse)
    [r(a)]      = findnearest(levelUse(a),table2array(T(:,2)));
    flow(a,1)   = table2array(T(r(a),4));
end

% extrapolate velocities beyond observed
n_extrap    = fitresult(levelUse(eval_idx));
vel_est     = ((levelUse(eval_idx).^(2/3)) .* (0.012.^(0.5))) ./ constant_n; % 

for a = 1:length(eval_idx)
    [r(a)]          = findnearest(levelUse(eval_idx(a)),interPoints);
    extrap_xsa(a,1) = areaOut(r(a));
    extrap_q(a,1)   = vel_est(a).*extrap_xsa(a,1);
end

% infill the flow record using the extrapolated values
flow(eval_idx) =  extrap_q;

% provide extended rating up to 2.2m 
h_ext       = 0.809:0.001:2.20;
vel_est     = ((h_ext.^(2/3)) .* (0.012.^(0.5))) ./ constant_n; % 

for a = 1:length(vel_est)
    [r(a)]              = findnearest(h_ext(a),interPoints);
    extrap_xsa_ext(a,1) = areaOut(r(a));
    extrap_q_ext(a,1)   = vel_est(a).*extrap_xsa_ext(a,1);
end

rating(:,1)     = [T{:,2}; h_ext']; % stage
rating(:,2)     = [T{:,4}; extrap_q_ext]; % q


%% create a figure summarising the workflow:
% Plot the cross-section data (subplot 1, labelled A)
f0          = figure(); hold on
ax0         = subplot(2,2,1);hold on
ax0.Color   = 'none';
set(ax0,'DefaultTextFontName','Arial')
axes(ax0); hold on

firstHalf   = realHgtIn(1:round(length(realHgtIn)./2));
[r(1)]      = findnearest(nanmax(levelUse),firstHalf);
secondHalf  = realHgtIn(round(length(realHgtIn)./2:end));
[r(2)]      = length(firstHalf) + findnearest(nanmax(levelUse),secondHalf) - 1;

h0 = plot(polyshape([distanceIn(r(1):r(2));distanceIn(r(1))] , [realHgtIn(r(1):r(2));realHgtIn(r(1))]) , ...
    'FaceColor', [254./255,232./255,200./255],... % red
    'FaceAlpha',1.0 ...
    ); hold on;


h1 = plot(polyshape([distanceIn(207:1523);distanceIn(207)] , [realHgtIn(207:1523);realHgtIn(207)]) , ...
    'FaceColor', [224./255,236./255,244./255],... % blue
    'FaceAlpha',1.0 ...
    ); hold on;

axis normal
axis tight

set(ax0,'xLim',[0 16.5])
set(ax0,'fontname','Arial')
set(ax0,'fontweight','normal')
set(ax0,'fontsize',11)
xlabel( 'Chainage [m]', 'Interpreter', 'none' );
ylabel( 'h [m]', 'Interpreter', 'none' );


%% Plot the flow gauging data and associated rating curve
% labelled B)
ax1 = subplot(2,2,2);hold on
axes(ax1); hold on
set(ax1,'DefaultTextFontName','Arial')
ax1.Color = 'none';

% median est
h3 = plot(table2array(T(:,4)),table2array(T(:,2)),...
    'Color', [0.5 0.5 0.5],...
    'LineWidth', 1);

% lower bounds
h4 = plot(table2array(T(:,3)),table2array(T(:,2)),...
    'Color', [0.5 0.5 0.5],...
    'LineWidth', 0.5,...
    'LineStyle', '--');

% upper bounds
h5 = plot(table2array(T(:,5)),table2array(T(:,2)),...
    'Color', [0.5 0.5 0.5],...
    'LineWidth', 0.5,...
    'LineStyle', '--');

% data scatter
h6 = scatter(infilled_gauging_q , infilled_gauging_level , ...
    'MarkerEdgeColor', [0.5 0.5 0.5],...
    'MarkerFaceColor', [224./255,236./255,244./255]...
    ); hold on;


axis normal
axis tight
set(ax1,'fontname','Arial')
set(ax1,'fontweight','normal')
set(ax1,'fontsize',11)
xlabel('Q [m^{3} s^{-1}]');
ylabel( 'h [m]', 'Interpreter', 'none' )


%% Plot the manning's calculations
ax2 = subplot(2,2,3);hold on
axes(ax2); hold on
set(ax2,'DefaultTextFontName','Arial')
ax2.Color = 'none';

% Plot fit with data.
h7 = plot(min(xData):0.01:max(xData),fitresult(min(xData):0.01:max(xData)),...
    'LineWidth', 1);
h7.Color = [0.5 0.5 0.5];

h8 = scatter(xData, yData_1,...
    'MarkerEdgeColor', [0.5 0.5 0.5],...
    'MarkerFaceColor', [224./255,236./255,244./255]...
    ); hold on;

xlabel( 'h [m]', 'Interpreter', 'none' );
ylabel( 'Mannings n [s m^{1/3}]');
set(ax2,'fontname','Arial')
set(ax2,'fontweight','normal')
set(ax2,'fontsize',11)


%% plot the extrapolated rating curve

rating_extrap_level = [max(infilled_gauging_level):0.01:max(levelUse)];
rating_extrap_vel = (rating_extrap_level.^(2/3)) .* (0.014.^(0.5)) ./ constant_n; % 

for a = 1:length(rating_extrap_vel)

    [r(a)] = findnearest(rating_extrap_level(a),interPoints);
    rating_extrap_xsa(a,1) = areaOut(r(a));
    rating_extrap_q(a,1)   = rating_extrap_vel(a).*rating_extrap_xsa(a,1);

end

ax3 = subplot(2,2,4);hold on
axes(ax3); hold on
set(ax3,'DefaultTextFontName','Arial')
ax3.Color = 'none';

% Plot fit with data.
h9 = plot(rating_extrap_q,rating_extrap_level,...
    'LineWidth', 1);
h9.Color = [0.5 0.5 0.5];

% median est
h10 = plot(table2array(T(:,4)),table2array(T(:,2)),...
    'Color', [0.5 0.5 0.5],...
    'LineWidth', 1);

% lower bounds
h11 = plot(table2array(T(:,3)),table2array(T(:,2)),...
    'Color', [0.5 0.5 0.5],...
    'LineWidth', 0.5,...
    'LineStyle', '--');

% upper bounds
h12 = plot(table2array(T(:,5)),table2array(T(:,2)),...
    'Color', [0.5 0.5 0.5],...
    'LineWidth', 0.5,...
    'LineStyle', '--');

% data scatter
h13 = scatter(infilled_gauging_q , infilled_gauging_level , ...
    'MarkerEdgeColor', [0.5 0.5 0.5],...
    'MarkerFaceColor', [224./255,236./255,244./255]...
    ); hold on;

axis normal
axis tight
set(ax3,'fontname','Arial')
set(ax3,'fontweight','normal')
set(ax3,'fontsize',11)
xlabel('Q [m^{3} s^{-1}]');
ylabel( 'h [m]', 'Interpreter', 'none' )

annotation(f0,'textbox',...
    [0.0376523551825029 0.926184492512128 0.0678571413670268 0.0690476176994188],...
    'String',{'A)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);

annotation(f0,'textbox',...
    [0.475152355182502 0.926184492512128 0.0678571413670268 0.0690476176994188],...
    'String',{'B)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);

annotation(f0,'textbox',...
    [0.0376523551825029 0.454755921083559 0.0678571413670268 0.0690476176994188],...
    'String',{'C)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);


annotation(f0,'textbox',...
    [0.475152355182502 0.454755921083559 0.0678571413670268 0.0690476176994188],...
    'String',{'D)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);


exportgraphics(f0,[data_out 'rating_workflow.png'],'Resolution',600)
