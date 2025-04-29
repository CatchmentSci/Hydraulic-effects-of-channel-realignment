% Step 1: Load the ASC file
clear all; close all; clc;

% directory containing data and scripts downloaded from github/zenodo. 
% should contain subfolders called 'Data' and 'Scripts'.
root_dir = 'C:\_git_local\Hydraulic-effects-of-channel-realignment\'; 

% specify the folder where generated data/outputs will be stored to.
data_out = [root_dir 'data\']; 

% ensure matlab can find required m files
addpath(genpath(root_dir)); 

% bring in the gauging data input file
fileIn      = [root_dir 'data\experimental run inputs1.xlsx'];
opts        = detectImportOptions(fullfile([fileIn]));
ii          = readtable(fileIn,opts);

% write the input hydrograph
idx = [2 6]; % sim B
q_in = table2array(ii(:,idx(1))) + table2array(ii(:,idx(2))); % sim B
q_in = q_in(~isnan(q_in));
q_in = q_in(289:end);

% specify the matlab data files
processed_dirName = [root_dir 'data\'];
processed_list = {'simB_pre_post.mat', 'simC-D.mat' };

for looper = 1:2

    load([processed_dirName processed_list{looper} ])

    if looper == 1
        f1 = figure(); hold on;
        pbaspect([2 1 1])
        ax1 = gca;
        set(ax1,'DefaultTextFontName','Arial')
        axes(ax1); hold on
        ax1.Color = 'none';
    end

    elapsed_hr = 5.*(1:length(sum(pre_floodplain_values(:,289:end))))./60;

if looper == 1
    h1(looper) = plot(elapsed_hr, sum(post_floodplain_values(:,289:end)),...
        'Color', [0.5 0.5 0.5]./255,...
        'LineWidth', 2,...
        'LineStyle', '-');

    simB_post = sum(post_floodplain_values(:,289:end));

else


    h1(looper) = plot(elapsed_hr, sum(floodplain_values_c(:,289:end)),...
        'Color', [215,25,28]./255,...
        'LineWidth', 2,...
        'LineStyle', '-');
    set(h1(looper), 'Color', [h1(looper).Color, 0.6]); % transparent


    h1(looper+1) = plot(elapsed_hr, sum(floodplain_values_d(:,289:end)),...
        'Color', [44,123,182]./255,...
        'LineWidth', 2,...
        'LineStyle', '-');
    set(h1(looper), 'Color', [h1(looper).Color, 0.6]); % transparent

    simC_post = sum(floodplain_values_c(:,289:end));
    simD_post = sum(floodplain_values_d(:,289:end));

end
end



axis tight
ylabel('Q [m^{3} s^{-1}]');
xlabel( 'Elapsed time [hr]', 'Interpreter', 'none' )
xticks(0:10:100);
yticks(0:20:80)


plot(polyshape([38 38 44 44] , [65 82 82 65 ]) , ...
    'FaceAlpha',0.0 ...
    ); hold on;

ax01 = axes('Position',[0.620 0.531 0.278 0.20]); hold on;

for looper = 1:2


    if looper == 1

        % create inset plot
        set(ax01,'DefaultTextFontName','Arial')
        ax01(looper).Color = 'none';
        box on

        h1(looper) = plot(elapsed_hr, sum(post_floodplain_values(:,289:end)),...
            'Color', [0.5 0.5 0.5]./255,...
            'LineWidth', 2,...
            'LineStyle', '-');

        xlim([39.5 42.5]);
        ylim ([65 82]);
        xticks(40:1:42);
        yticks(65:5:80)
    
    else

        h1(looper) = plot(elapsed_hr, sum(floodplain_values_c(:,289:end)),...
            'Color', [215,25,28]./255,...
            'LineWidth', 2,...
            'LineStyle', '-');
        set(h1(looper), 'Color', [h1(looper).Color, 0.6]); % transparent


        h1(looper+1) = plot(elapsed_hr, sum(floodplain_values_d(:,289:end)),...
            'Color', [44,123,182]./255,...
            'LineWidth', 2,...
            'LineStyle', '-');
        set(h1(looper), 'Color', [h1(looper).Color, 0.6]); % transparent

    end


end


% do xcorr analysis
up_in    = q_in(:,:);
down_in1 = simB_post';
down_in2 = simC_post';
down_in3 = simD_post';

transmission_t(1) = cross_correlation_sims (up_in, down_in1);
transmission_t(2) = cross_correlation_sims (up_in, down_in2);
transmission_t(3) = cross_correlation_sims (up_in, down_in3);

% perform hysteresis analysis
[p3(1), polygon_area(1)] = hysteresis_sims (up_in, down_in2) ; % sim c
[p3(2), polygon_area(2)] = hysteresis_sims (up_in, down_in3) ; % sim d

% calculate phase of signals
% specify points
wave_duration = [50];
starting_point = [15];

% first for simC
stage_file      = [root_dir 'data\Run_post_c.stage'];
opts            = detectImportOptions(stage_file, 'FileType', 'delimitedtext' );
jj              = readtable(stage_file,opts);
stage_c         = table2array(jj(4:end,2));

[val, idx]      = nanmax(stage_c); 
peak_stage_idx  = (idx./12) - (289./12) - starting_point;
[Qp, q_idx]     = nanmax(sum(floodplain_values_c(:,:)));
peak_q_idx      = (q_idx ./ 12) - (289./12) - starting_point;
phi(1)          = calculate_phase(wave_duration, peak_stage_idx, peak_q_idx);


% then for SimD
stage_file      = [root_dir 'data\Run_post_d.stage'];
opts            = detectImportOptions(stage_file, 'FileType', 'delimitedtext' );
jj              = readtable(stage_file,opts);
stage_d         = table2array(jj(4:end,2));

[val, idx]      = nanmax(stage_d); 
peak_stage_idx  = (idx./12) - (289./12) - starting_point;
[Qp, q_idx]     = nanmax(sum(floodplain_values_d(:,:)));
peak_q_idx      = (q_idx ./ 12) - (289./12) - starting_point;
phi(2)          = calculate_phase(wave_duration, peak_stage_idx, peak_q_idx);

% export the figure
exportgraphics(f1,['simC-D_out.png'],'Resolution',600)






