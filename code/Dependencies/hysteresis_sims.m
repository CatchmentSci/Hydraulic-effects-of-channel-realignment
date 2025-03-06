function [p3, polygon_area] = hysteresis_sims (up_in, down_in) 

% dirIn  = 'D:\OneDrive - Newcastle University\Documents - Goldrill Beck Research\General\_shared\Hydro Data\PT_processed\';
% dirOut = 'D:\OneDrive - Newcastle University\Documents - Goldrill Beck Research\General\_shared\Hydro Data\PT_processed\';
% hysteresis(dirIn, dirOut)



%% undertake hysteresis analysis
array1{1,1}     = down_in;
half_way        = floor(length(array1{1,1})./2);
array1_1half    = array1{1,1}(1:half_way);
array1_2half    = array1{1,1}(half_way+1:half_way.*2);
low_value       = max(min([array1_1half,array1_2half])); % only need the low value for one array
idx2            = [min(find(array1{1,1} > low_value == 1)); max(find(array1{1,1} > low_value == 1))]; % min and max idx for data

array2{1,1}     = up_in;
array2_1half    = array2{1,1}(1:half_way);
array2_2half    = array2{1,1}(half_way+1:half_way.*2);

% clip the data so that the miniumum values are present on the rising
% and falling limb of the hydrograph
h{1,1}         = array2{1,1}(idx2(1):idx2(2));
h2{1,1}        = array1{1,1}(idx2(1):idx2(2));

% normalise
Hmin{1,1}      = min(h{1,1});
Hmax{1,1}      = max(h{1,1});
h{1,1}         = (h{1,1}-min(h{1,1})) ./ (Hmax{1,1} - Hmin{1,1});

% normalise
Hmin2{1,1}      = min(h2{1,1});
Hmax2{1,1}      = max(h2{1,1});
h2{1,1}         = (h2{1,1}-min(h2{1,1})) ./ (Hmax2{1,1} - Hmin2{1,1});


% calculate the hysteresis
for b = 1:length(h2{1,1})
    if b == 1
        p0_1(b,1) = h{1,1}(b,1) - 0;
        p0_2(b,1) = h2{1,1}(b,1) - 0;
        p1(b,1) = (h2{1,1}(b,1).*(p0_1(b,1)./300)) - (h{1,1}(b,1) .* (p0_2(b,1)./300));

    else
        p0_1(b,1) = h{1,1}(b,1) - h{1,1}(b-1,1);
        p0_2(b,1) = h2{1,1}(b,1) - h2{1,1}(b-1,1);
        p1(b,1) = h2{1,1}(b,1).*(p0_1(b,1)./300) - h{1,1}(b,1) .* (p0_2(b,1)./300);
    end
    p2(b,1) = p1(b,1) .* 300;
    if b == length(h2{1,1})
        p3(1) = 0.5.*sum(p2);
    end
end
clear p1 p0_1 p0_2 p1 p2

% alternative method
vertices = [h{1,1},h2{1,1}];  % Coordinates of vertices of a rectangle
polygon_area(1,1) = calculate_polygon_area(vertices);

