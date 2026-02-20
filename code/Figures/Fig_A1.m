
%XS A_DIFF_CALCS Compare bankfull cross-sections (DTM vs MS50) for 3 sections
% Produces 3 subplots (one per section), y-axis inverted (depth downward),
% and prints Δ bankfull area on each subplot.

% ---- INPUT ----
fileIn = 'D:\OneDrive - Newcastle University\Documents - Goldrill Beck Research\General\_shared\Lisflood-FP\Topo Comparison\DTM_MS50_comp.xlsx';
T      = readtable(fileIn,'Sheet','Summary');

% Each row: [x_dtm z_dtm x_ms50 z_ms50]
colGroups = {
    [ 1  2  4  5]   % S1
    [ 7  8 10 11]   % S2
    [13 14 16 17]   % S3
};

titles = {'Upper','Lower','Redundancy'};

% ---- FIGURE ----
figure('Color','w');
tlo = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

difference = nan(1,3);

% ---- LOOP ----
for i = 1:3
    ax = nexttile(tlo,i); hold(ax,'on'); box(ax,'on');
    set(ax, 'DefaultAxesFontName', 'Arial'); 
    set(ax,'fontsize', 12);

    s = T{:, colGroups{i}};

    xD = s(:,1); zD = s(:,2);
    xM = s(:,3); zM = s(:,4);

    % ---- DTM ----
    idxD = ~isnan(xD) & ~isnan(zD);
    yD   = -zD + abs(min(-zD(idxD)));   % your original "invert + shift" method
    areaD = trapz(xD(idxD), yD(idxD));

    % ---- MS50 ----
    idxM = ~isnan(xM) & ~isnan(zM);
    yM   = -zM + abs(min(-zM(idxM)));
    areaM = trapz(xM(idxM), yM(idxM));

    % ---- Differences ----
    dA = areaD - areaM;                          % absolute (m^2)
    difference(i) = (dA ./ areaM) * 100;         % relative (%), vs MS50

    % ---- Plot ----
    plot(ax, xD(idxD), yD(idxD), 'LineWidth',1.5, 'Color',[215,25,28]./255);
    plot(ax, xM(idxM), yM(idxM), 'LineWidth',1.5, 'Color',[44,123,182]./255);

    set(ax,'YDir','reverse'); % deeper = larger y

    xlabel(ax,'Distance (m)');
    ylabel(ax,'Depth below bankfull (m)');
    title(ax, titles{i});
    if i == 2
        legend(ax, {'DTM','MS50'}, 'Location','best');
    end

    % ---- Annotation (stays in corner regardless of axis limits) ----
    txt = sprintf('\\DeltaA_{bf} = %.3f m^2 (%.1f%%)', dA, difference(i));
    text(ax, 0.02, 0.95, txt, ...
        'Units','normalized', 'FontWeight','bold', ...
        'VerticalAlignment','top');
end

exportgraphics(gcf, 'xs_comp.pdf', 'ContentType','vector');
