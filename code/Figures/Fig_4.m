% Step 1: Load the ASC file
clear all; close all; clc;

% directory containing data and scripts downloaded from github/zenodo. 
% should contain subfolders called 'Data' and 'Scripts'.
root_dir = 'D:\OneDrive - Newcastle University\Documents - Goldrill Beck Research\General\_shared\Submission Docs\'; 

% specify the folder where generated data/outputs will be stored to.
data_out = [root_dir 'Data\']; 

% ensure matlab can find required m files
addpath(genpath(root_dir)); 

% bring in the gauging data input file
fileIn      = [root_dir 'Data\experimental run inputs1.xlsx'];
opts        = detectImportOptions(fullfile([fileIn]));
ii          = readtable(fileIn,opts);

% specify the matlab data files
processed_list = {'simA_pre_post.mat', 'simB_pre_post.mat' };

for looper = 1:2

    load([root_dir 'Data\' processed_list{looper} ])

    if looper == 1 % sim A
        idx = [13 17];

        % create the figure
        f0 = figure(); hold on
        ax0 = subplot(2,1,1);hold on
        set(ax0,'DefaultTextFontName','Arial')
        axes(ax0); hold on
        ax0.Color = 'none';


    else
        ax1 = subplot(2,1,2);hold on
        set(ax1,'DefaultTextFontName','Arial')
        axes(ax1); hold on
        ax1.Color = 'none';

        idx = [2 6]; % sim B
       
    end

    % write the input hydrograph
    q_in = table2array(ii(:,idx(1))) + table2array(ii(:,idx(2))); % sim B
    q_in = q_in(~isnan(q_in));

    elapsed_hr = 5.*(1:length(q_in(289:end,:)))./60;


    h0(looper) = plot(elapsed_hr, q_in(289:end,:),...
        'Color', [0.5 0.5 0.5],...
        'LineWidth', 2,...
        'LineStyle', '-');

    h1(looper) = plot(elapsed_hr, sum(pre_floodplain_values(:,289:end)),...
        'Color', [215,25,28]./255,...
        'LineWidth', 2,...
        'LineStyle', '-');
    set(h1(looper), 'Color', [h1(looper).Color, 0.6]); % transparent

    h2(looper) = plot(elapsed_hr, sum(post_floodplain_values(:,289:end)),...
        'Color', [44,123,182]./255,...
        'LineWidth', 2,...
        'LineStyle', '-');
    set(h2(looper), 'Color', [h2(looper).Color, 0.6]); % transparent

    axis tight
    ylabel('Q [m^{3} s^{-1}]');
    if looper == 2
    xlabel( 'Elapsed time [hr]', 'Interpreter', 'none' )
    end

    % create inset box
    if looper == 1
        plot(polyshape([147 147 154 154] , [38 43 43 38 ]) , ...
            'FaceAlpha',0.0 ...
            ); hold on;

    ax01(looper) = axes('Position',[.75 .75 .2 .2]); hold on;

    else
        plot(polyshape([38 38 46 46] , [55 82 82 55 ]) , ...
            'FaceAlpha',0.0 ...
            ); hold on;

     ax01(looper) = axes('Position',[0.626 0.283 0.278 0.20]); hold on;

    end

    % create inset plot
    set(ax01(looper),'DefaultTextFontName','Arial')
    ax01(looper).Color = 'none';
    box on
    h00(looper) = plot(elapsed_hr, q_in(289:end,:),...
        'Color', [0.5 0.5 0.5],...
        'LineWidth', 2,...
        'LineStyle', '-');

    h01(looper) = plot(elapsed_hr, sum(pre_floodplain_values(:,289:end)),...
        'Color', [215,25,28]./255,...
        'LineWidth', 2,...
        'LineStyle', '-');

    set(h01(looper), 'Color', [h01(looper).Color, 0.6]); % transparent

    h02(looper) = plot(elapsed_hr, sum(post_floodplain_values(:,289:end)),...
        'Color', [44,123,182]./255,...
        'LineWidth', 2,...
        'LineStyle', '-');
        
    set(h02(looper), 'Color', [h02(looper).Color, 0.6]); % transparent


    if looper == 1
        xlim([149 152]);
        ylim ([38 43]);
    else
        xlim([38 46]);
        ylim ([55 82]);
    end

% pull out the required stats

ranger = ([25, 75, 80, 120, 140, 190 ] *12) + 289; % ignore spin up period
    if looper == 1
        for b = 1:3
            % calculate Qp
            [Qp{looper}(b,1), q_idx{looper}(b,1)] = nanmax(sum(pre_floodplain_values(:,ranger((b*2)-1):ranger(b*2))));
            [Qp{looper}(b,2), q_idx{looper}(b,2)] = nanmax(sum(post_floodplain_values(:,ranger((b*2)-1):ranger(b*2))));
            q_idx{looper} = q_idx{looper} + ranger((b*2)-1)-1;

            % do xcorr analysis
            up_in   = q_in(ranger((b*2)-1):ranger(b*2),:);
            down_in1 = sum(pre_floodplain_values(:,ranger((b*2)-1):ranger(b*2)))';
            down_in2 = sum(post_floodplain_values(:,ranger((b*2)-1):ranger(b*2)))';
            [maxVar{looper}(b,1), peak_stage{looper}(b,1)] = cross_correlation_sims (up_in, down_in1);
            [maxVar{looper}(b,2), peak_stage{looper}(b,2)] = cross_correlation_sims (up_in, down_in2);
           

            if b == 3
                maxVar = cellfun(@(x) x*5,maxVar,'un',0);  % convert to mins

                % calculate the celerity
                maxVar_s = cellfun(@(x) x*60,maxVar,'un',0);  % convert to seconds
                temp = cellfun(@(x) 1350./x, maxVar_s, 'un',0);  % pre calc
                temp2 = cellfun(@(x) 1460./x, maxVar_s, 'un',0);  % post calc
                celerity{1,1}(:,1) = temp{1,1}(:,1) ; % pre sim A
                celerity{1,1}(:,2) = temp2{1,1}(:,2) ; % post sim A

            end

            % perform hysteresis analysis
            [p3{looper}(b,1), polygon_area{looper}(b,1)] = hysteresis_sims (up_in, down_in1) ;
            [p3{looper}(b,2), polygon_area{looper}(b,2)] = hysteresis_sims (up_in, down_in2) ;

            % calculate phase of signals
            % first for the pre-condition
            stage_file_pre = [root_dir 'Data\Run_pre_a.stage'];
            opts = detectImportOptions(stage_file_pre, 'FileType', 'delimitedtext' );
            jj = readtable(stage_file_pre,opts);
            stage_pre = table2array(jj(4:end,2));

            % specify points
            wave_duration{looper} = [40, 40, 50];
            starting_point{looper} = [25, 80, 140];

            [val, idx]      = nanmax(stage_pre(ranger((b*2)-1):ranger(b*2))); % pull out the peak stage (pre)
            peak_stage_idx  = (idx + ranger((b*2)-1)-1)./12 - (289./12);
            peak_q_idx      = (q_idx{looper}(b,1)+0) ./ 12 - (289./12) ;
            phi{looper}(b,1) = calculate_phase(wave_duration{looper} (b), peak_stage_idx(1), peak_q_idx);
            
            % then for the post-condition
            stage_file_post = [root_dir 'Data\Run_post_a.stage'];
            opts = detectImportOptions(stage_file_post, 'FileType', 'delimitedtext' );
            kk = readtable(stage_file_post,opts);
            stage_post = table2array(kk(4:end,2));

            [val, idx]      = nanmax(stage_post(ranger((b*2)-1):ranger(b*2))); % pull out the peak stage (post)
            peak_stage_idx  = (idx + ranger((b*2)-1))./12 - (289./12);
            peak_q_idx      = (q_idx{looper}(b,2)+1) ./ 12 - (289./12) ;
            phi{looper}(b,2) = calculate_phase(wave_duration{looper}(b), peak_stage_idx(1), peak_q_idx);


            % pull out the epoch needed to extract distributed depths and velocities at the mean stage
            if b < 3

                t1 = ranger((b*2)-1):ranger(b*2);
                idx_epoch(b) = t1(1) + min(find(abs(stage_post(t1) - nanmean(stage_post(t1))) < 0.01)); % find within 1cm

            end


        end

    else
        for b = 1
            [Qp{looper}(b,1), q_idx{looper}(b,1)]  = nanmax(sum(pre_floodplain_values(:,:)));
            [Qp{looper}(b,2), q_idx{looper}(b,2)]  = nanmax(sum(post_floodplain_values(:,:)));
            q_idx{looper} = q_idx{looper};

            % do xcorr analysis
            up_in   = q_in(:,:);
            down_in1 = sum(pre_floodplain_values(:,:))';
            down_in2 = sum(post_floodplain_values(:,:))';
            [maxVar{looper}(b,1), peak_stage{looper}(b,1)] = cross_correlation_sims (up_in, down_in1);
            [maxVar{looper}(b,2), peak_stage{looper}(b,2)] = cross_correlation_sims (up_in, down_in2);


            % celerity calcs
            maxVar = cellfun(@(x) x*5,maxVar,'un',0);  % convert to mins
            % calculate the celerity
            maxVar_s = cellfun(@(x) x*60,maxVar,'un',0);  % convert to seconds
            temp = cellfun(@(x) 1350./x, maxVar_s, 'un',0);  % pre calc
            temp2 = cellfun(@(x) 1460./x, maxVar_s, 'un',0);  % post calc
            celerity{1,2}(:,1) = temp{1,2}(:,1) ; % pre sim B
            celerity{1,2}(:,2) = temp2{1,2}(:,2) ; % post sim B

            % perform hysteresis analysis
            [p3{looper}(b,1), polygon_area{looper}(b,1)] = hysteresis_sims (up_in, down_in1) ;
            [p3{looper}(b,2), polygon_area{looper}(b,2)] = hysteresis_sims (up_in, down_in2) ;

            % calculate phase of signals
            % first for the pre-condition
            stage_file_pre = [root_dir 'Data\Run_pre_b.stage'];
            opts = detectImportOptions(stage_file_pre, 'FileType', 'delimitedtext' );
            jj = readtable(stage_file_pre,opts);
            stage_pre = table2array(jj(4:end,2));

            % specify points
            wave_duration{looper} = [50];
            starting_point{looper} = [15];

            [val, idx]      = nanmax(stage_pre); % pull out the peak stage (pre)
            peak_stage_idx  = (idx./12) - (289./12) - starting_point{looper};
            peak_q_idx      = (q_idx{looper}(b,1)+0) ./ 12 - (289./12) - starting_point{looper};
            phi{looper}(b,1) = calculate_phase(wave_duration{looper} (b), peak_stage_idx(1), peak_q_idx);

            % then for the post-condition
            stage_file_post = [root_dir 'Data\Run_post_b.stage'];
            opts = detectImportOptions(stage_file_post, 'FileType', 'delimitedtext' );
            kk = readtable(stage_file_post,opts);
            stage_post = table2array(kk(4:end,2));

            [val, idx]      = nanmax(stage_post); % pull out the peak stage (post)
            peak_stage_idx  = (idx./12) - (289./12) - starting_point{looper};
            peak_q_idx      = (q_idx{looper}(b,2)+0) ./ 12 - (289./12) - starting_point{looper};
            phi{looper}(b,2) = calculate_phase(wave_duration{looper} (b), peak_stage_idx(1), peak_q_idx);




        end

    end

end

% add sub-plot annotations to finish
annotation(f0,'textbox',...
    [0.0617142857142855 0.907142857142861 0.0615 0.0476190476190478],...
    'String',{'A)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);


annotation(f0,'textbox',...
    [0.0617142857142855 0.426190477538687 0.0678571413670268 0.0690476176994188],...
    'String',{'B)'},...
    'LineStyle', 'none',...
    'FitBoxToText','on',...
    'fontsize',11);

% export the figure
exportgraphics(f0,['simA-B_out.png'],'Resolution',600)


