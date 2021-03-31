%Funciton for calculating segment medians
function region_median = calc_region_median(segs,bound1,bound2,frac)
    [~, ix] = sort(segs.Segment_Mean);
    %sort from lowest to highest on segment mean
    segs = trimStruct(segs,ix);
    %correct length for parts of segment that do not overlap the region of
    %interest
    segs.length = max((min(segs.gend,bound2) - max(segs.gstart,bound1)),0); %+ sum(counts.total_length);
    %find midpoint of all parts of segment that overlap the region, assuming start of region is 0
    segs=trimStruct(segs,segs.length ~=0);
    medianix = (sum(segs.length))/frac;
    pos = 0;
    for i = 1:length(segs.length);
        %must update position before checking if we've passed median
        %also must be outside if statement
        pos = pos+segs.length(i);
        if pos >= medianix;
            region_median = segs.Segment_Mean(i);
            return;
        end
    end
end