

for a = 1:length(peak_stage_all)
    if ~isnan(peak_stage_all(a,1))
        idx(a) = find(min(abs(peak_stage_all(a,1) - rating(:,1))) == abs(peak_stage_all(a,1) - rating(:,1))) ;
        peak_q_all(a,1) = rating(idx(a),2);
    else
        peak_q_all(a,1) = NaN;
    end
end




for aa = 1:228

    neg(aa,1) = find ( R(:,aa) == max(R(:,aa))) < 12;
    if neg(aa,1) == true
        q(aa,1) = peak_q(aa,1);
    else 
         q(aa,1) = NaN;
    end

    poor_corr(aa,1) = max(R(:,aa)) < 0.90;

    if poor_corr(aa,1) == true
        q(aa,2) = peak_q(aa,1);
    else 
         q(aa,2) = NaN;
    end
end