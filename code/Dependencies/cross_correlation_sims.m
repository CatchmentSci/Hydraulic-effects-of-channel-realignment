

function [maxVar, peak_stage] = cross_correlation_sims (up_in, down_in) 
%% undertake some xcorr analysis
clear R maxVar peak_stage idx


    array1{1,1}    = down_in;
    array2{1,1}    = up_in;

    % modification to pull out peak values for xcorr
    peak_idx        = find(max(array2{1,1}) == array2{1,1});
    ranger = 60;
    if peak_idx < ranger + 1
        array2{1,1}    = array2{1,1}(1:ranger);
        array1{1,1}    = array1{1,1}(1:ranger);
    elseif peak_idx + ranger > length(array2{1,1})
        array2{1,1}    = array2{1,1}(peak_idx-ranger:end);
        array1{1,1}    = array1{1,1}(peak_idx-ranger:end);
    else
        array2{1,1}    = array2{1,1}(peak_idx-ranger:peak_idx+ranger);
        array1{1,1}    = array1{1,1}(peak_idx-ranger:peak_idx+ranger);
    end

    % remove start data for s1
    padder1 = repmat(NaN,12,1);
    array2{1,1} = [array2{1,1}; padder1];
    array1{1,1} = [padder1; array1{1,1}];

    for a = 1:60
        shifter = a;
        prepData = zeros(1,shifter)';
        prepData = replace_num(prepData,0,NaN);
        xIn = [array1{1,1};prepData];
        yIn = [prepData;array2{1,1}];
        [xIn,yIn] = prepareCurveData(xIn,yIn);
        if ~isempty(xIn)
            R(a,1) = corr(xIn,yIn);
        else
            R(1:100,1) = NaN;
        end
        %max_corr(1,1) = max(R(:,1))
    end
    if ~isnan(R(1,1)) && max(R(:,1)) >0.90 && find ( R(:,1) == max(R(:,1))) > 11 % need an r2 of >0.90 to be accepted
        peak_stage(1,1) = nanmax(array1{1,1}); % peak flow at site 2
        maxVar(1,1) = find ( R(:,1) == max(R(:,1)));
        maxVar(1,1) = maxVar(1,1) - length(padder1);
    else
        maxVar(1,1) = NaN;
        peak_stage(1,1) = NaN;
    end

end
