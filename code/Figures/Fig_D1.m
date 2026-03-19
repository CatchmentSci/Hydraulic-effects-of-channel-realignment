%% Plot selected storm hydrographs (Site 1 vs Site 2)
% This script is called from within Fig_3_4_5_D1.m
% As default, the call in the aforementioned script is commented out.
% To enable this script and generate Figure D1, you must uncommnt line 189
% within Fig_3_4_5_D1.m

purpose = 2; % 1 =  review 2 = supp info

% --- Ensure times are datetime
t1 = datetime(s1_date_matched);   % upstream time
h1 = s1_level_matched;            % upstream stage
t2 = datetime(s2_date_matched);   % downstream time
h2 = s2_level_matched;            % downstream stage

% --- Plot the six events where downstream stage exceeds 1.86 m
event_nums = [7, 76, 89, 90, 112, 135];

% Build events struct from event_period (row = event number)
events = struct([]);
for i = 1:numel(event_nums)
    n = event_nums(i);

    events(i).name  = sprintf('Event %d', n);
    events(i).start = event_period(n,1);
    events(i).end   = event_period(n,2);
end

% --- Optional plot window offsets (hours) for specific events
plotOffsetHours = containers.Map('KeyType','double','ValueType','double');
plotOffsetHours(89) = 48;   % shift start forward by 48 hours
plotOffsetHours(90) = 12;   % shift start forward by 12 hours


% --- Create figure: half screen width, full height (left side)
f = figure('Color','w', 'Units','normalized') ...%, 'Position',[0 0 0.5 1]);

if purpose == 1
    tl = tiledlayout(f, 3, 2, 'TileSpacing','compact', 'Padding','compact');
elseif purpose == 2
    tl = tiledlayout(f, 2, 2, 'TileSpacing','compact', 'Padding','compact');
end

% sample interval (minutes)
dt_minutes = minutes(t1(2) - t1(1));

% store handles for legend (grab from first subplot with data)
hLegendSite1 = [];
hLegendSite2 = [];

if purpose == 1
    include = [1,2,3, 4,5,6];
elseif purpose == 2
    include = [1,2,4,5];
end

for k = include

    ax = nexttile(tl);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    set(ax,'DefaultTextFontName','Arial')
    set(ax,'fontsize', 12);

    % Default plot window
    plotStart = events(k).start;
    plotEnd   = events(k).end;

    % Apply special offsets if required
    eventNum = sscanf(events(k).name, 'Event %d');
    if isKey(plotOffsetHours, eventNum)
        plotStart = plotStart + hours(plotOffsetHours(eventNum));
    end

    mask1 = (t1 >= plotStart) & (t1 <= plotEnd) & ~isnan(h1);
    mask2 = (t2 >= plotStart) & (t2 <= plotEnd) & ~isnan(h2);


    % Plot Site 1 (upstream)
    if any(mask1)
        hS1 = plot(ax, t1(mask1), h1(mask1), '-',...
            'LineWidth', 1.2,...
            'Color',  [215,25,28]./255 ...
            );
    else
        hS1 = plot(ax, nan, nan, '-',...
            'LineWidth', 1.2,...
            'Color',  [215,25,28]./255 ...
            );
    end

    % Plot Site 2 (downstream)
    if any(mask2)
        hS2 = plot(ax, t2(mask2), h2(mask2), '-', ...
            'LineWidth', 1.2,...
            'Color', [44,123,182]./255 ...
            );
    else
        hS2 = plot(ax, nan, nan, '-', ...
            'LineWidth', 1.2,...
            'Color',   [44,123,182]./255 ...
            );
    end

    % Save handles for the legend (once)
    if isempty(hLegendSite1) && any(mask1)
        hLegendSite1 = hS1;
    end
    if isempty(hLegendSite2) && any(mask2)
        hLegendSite2 = hS2;
    end

    if purpose == 1
        title(ax, events(k).name, 'Interpreter','none');
    end
    ylabel(ax, 'Stage [m]');
    xlabel(ax, 'Time');

    % X limits from data
    if any(mask1) || any(mask2)
        tmin = min([t1(mask1); t2(mask2)]);
        tmax = max([t1(mask1); t2(mask2)]);
        if tmin < tmax
            xlim(ax, [tmin tmax]);
        end
    end

    xtickangle(ax,45);

    % ---- Xcorr transmission time textbox ----
    eventStart = events(k).start;
    eventEnd   = events(k).end;

    overlap = min(event_period(:,2), eventEnd) - max(event_period(:,1), eventStart);
    [~, aa] = max(overlap);

    if aa >= 1 && aa <= numel(maxVar) && ~isnan(maxVar(aa))
        ttMin = maxVar(aa) * dt_minutes;
        bestR = max(R(:,aa), [], 'omitnan');
        txt = sprintf('Transmission time: %.0f min\nmax corr: %.2f', ttMin, bestR);
    else
        txt = 'Transmission time: n/a';
    end

    text(ax, 0.02, 0.98, txt, ...
        'Units','normalized', ...
        'VerticalAlignment','top', ...
        'BackgroundColor','w', ...
        'Margin',6);

    % ---- Vertical lines at first exceedance and last return below 1.86 m ----
    thr = 1.86;

    if any(mask2)
        t2w = t2(mask2);
        h2w = h2(mask2);

        exceed = h2w > thr;
        if any(exceed)
            d = diff([false; exceed(:); false]);

            % First exceedance
            iUp = find(d == 1, 1, 'first');
            if ~isempty(iUp)
                xline(ax, t2w(iUp), '--', ...
                    'LineWidth', 1.2,...
                    'Color',  [0.5 0.5 0.5] ...
                    );
            end

            % Last return below threshold
            iDownAll = find(d == -1);
            if ~isempty(iDownAll)
                iDown = iDownAll(end) - 1;
                if iDown + 1 <= numel(t2w)
                    xline(ax, t2w(iDown+1), '--', ...
                    'LineWidth', 1.2,...
                    'Color',  [0.5 0.5 0.5] ...
                    );
                end
            end
        end
    end

    legend(ax, {'Site 1 (upstream)','Site 2 (downstream)'}, ...
        'Location','southeast', ...
        'Box','off');

end


% --- Export ---
%exportgraphics(f, 'Fig_storm_hydrographs.pdf', 'ContentType','vector');
%exportgraphics(f, 'Fig_storm_hydrographs.png', 'Resolution',600);
