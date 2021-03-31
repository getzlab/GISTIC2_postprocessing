function  segs = apply_cnv_blacklist(segs,cnv_blacklist,AF,DF,ARMS)
    cnv_blacklist.gstart = xhg19(cnv_blacklist.Chromosome,cnv_blacklist.Start);
    cnv_blacklist.gend   = xhg19(cnv_blacklist.Chromosome,cnv_blacklist.End  );
    cnv_blacklist.keep_event = zeros(size(cnv_blacklist.Chromosome));
    %arm_level_significance = trimStruct(arm_level_significance,arm_level_significance.significant_amplification | arm_level_significance.significant_deletion);
    AF.length = AF.gend - AF.gstart;
    DF.length = DF.gend - DF.gstart;
    cnv_blacklist.length = cnv_blacklist.gend - cnv_blacklist.gstart;
    ARMS.length = ARMS.gend - ARMS.gstart;
    for i = 1:length(AF.gstart)
        overlapix = find(((cnv_blacklist.gstart <= AF.gstart(i) & cnv_blacklist.gend >= AF.gstart(i)) | (cnv_blacklist.gstart <= AF.gend(i) & cnv_blacklist.gend >= AF.gend(i)) | (cnv_blacklist.gstart >= AF.gstart(i) & cnv_blacklist.gend <= AF.gend(i))) & (cnv_blacklist.length > AF.length(i)*0.01));
        cnv_blacklist.keep_event(overlapix,1) = 1;
    end
   for i = 1:length(DF.gstart)
        overlapix = find(((cnv_blacklist.gstart <= DF.gstart(i) & cnv_blacklist.gend >= DF.gstart(i)) | (cnv_blacklist.gstart <= DF.gend(i) & cnv_blacklist.gend >= DF.gend(i)) | (cnv_blacklist.gstart >= DF.gstart(i) & cnv_blacklist.gend <= DF.gend(i))) & (cnv_blacklist.length > DF.length(i)*0.01));
        cnv_blacklist.keep_event(overlapix,1) = 1;
    end
    for i = 1:length(ARMS.gstart)
        overlapix = find(((cnv_blacklist.gstart <= ARMS.gstart(i) & cnv_blacklist.gend >= ARMS.gend(i)) | (cnv_blacklist.gstart <= ARMS.gstart(i) & cnv_blacklist.gend >= ARMS.gend(i)) | (cnv_blacklist.gstart >= ARMS.gstart(i) & cnv_blacklist.gend <= ARMS.gend(i))) & ((cnv_blacklist.length > ARMS.length(i)*0.01)));
        cnv_blacklist.keep_event(overlapix,1) = 1;
    end
    cnv_blacklist = trimStruct(cnv_blacklist,cnv_blacklist.keep_event);
    segs.remove_seg = zeros(size(segs.Chromosome));
    
    for i = 1:length(cnv_blacklist.Chromosome)
        segix = find((segs.gstart <= cnv_blacklist.gstart(i) & segs.gend >= cnv_blacklist.gstart(i)) | (segs.gstart <= cnv_blacklist.gend(i) & segs.gend >= cnv_blacklist.gend(i)) | (segs.gstart >= cnv_blacklist.gstart(i) & segs.gend <= cnv_blacklist.gend(i)) & ~segs.remove_seg) ;
        newsegs = trimStruct(segs,zeros(size(segs.Chromosome)));
        for j = 1:length(segix)
            %if segment spans entire blacklisted region, remove if
            %blacklisted region is >80% of segment length
            %segments
            if (segs.gstart(segix(j)) < cnv_blacklist.gstart(i) & segs.gend(segix(j)) > cnv_blacklist.gend(i)) 
                cnv_length = cnv_blacklist.gend(i) - cnv_blacklist.gstart(i);
                seg_length = segs.length(segix(j));
                if (cnv_length/seg_length > 0.8)
                    segs.remove_seg(segix(j)) = 1;
                end
            %if segments spans just start of blacklisted region, set end of
            %segment to start of blacklisted region
            elseif (segs.gstart(segix(j)) < cnv_blacklist.gstart(i) & segs.gend(segix(j)) > cnv_blacklist.gstart(i))
                    segs.gend(segix(j))   = cnv_blacklist.gstart(i);
                    segs.End(segix(j)) = cnv_blacklist.Start(i);
            %if segment spans just end of blacklisted region, set start to
            %end of blacklisted region
            elseif (segs.gstart(segix(j))   < cnv_blacklist.gend(i) & segs.gend(segix(j)) > cnv_blacklist.gend(i)) 
                    segs.gstart(segix(j))   = cnv_blacklist.gend(i);
                    segs.Start(segix(j)) = cnv_blacklist.End(i);
            %if segment falls entirely within blacklisted region, remove
            %segment entirely
            elseif (segs.gstart(segix(j)) >= cnv_blacklist.gstart(i) & segs.gend(segix(j)) <= cnv_blacklist.gend(i))
                    segs.remove_seg(segix(j)) = 1;
            end
        end
        segs = trimStruct(segs,~segs.remove_seg);
        segs = mergeStruct(newsegs,segs);
    end
end