function [A, F]=GISTIC2_force_calling(segfile,amp_focal_file,del_focal_file,cnv_blacklist_file, arm_coordinates_file, PAR) 
if nargin<6
    PAR=[]
    PAR.single_amp_threshold = 0.1;
    PAR.single_del_threshold = -0.1;
    PAR.double_amp_threshold = 0.9;
    PAR.double_del_threshold = -0.9;
    PAR.arm_length_fraction_threshold = 2;
end
if nargin<5
    % assume all arms 
    arm_coordinates_file=''
end

segs=load_tsv(segfile);
%AL=load_tsv(focal_file);

AF=load_tsv(amp_focal_file);
AF=trimStruct(AF,1:3); AF=rmfield(AF,'N'); AF=flipStruct(AF); N=length(AF.cytoband); AF=trimStruct(AF,1:(N-1));
AF.chr=str2double(cellfun(@(x) x(1),regexp(regexprep(AF.widePeakBoundaries,'chr',''),':','split')));
x=cellfun(@(x) x(2),regexp(AF.widePeakBoundaries,':','split'));
AF.Start=str2double(cellfun(@(x) x(1),regexp(x,'-','split')));
AF.End=str2double(cellfun(@(x) x(2),regexp(x,'-','split')));
AF.gstart=xhg19(AF.chr,AF.Start);
AF.gend=xhg19(AF.chr,AF.End);
AF.cytoband=erase(AF.cytoband,"x");
AF.cytoband=strcat(AF.cytoband,':AMP');
AF.Descriptor=AF.cytoband;
%F.gene = AF.cytoband;
%AF.Descriptor=


DF=load_tsv(del_focal_file);
DF=trimStruct(DF,1:3); DF=rmfield(DF,'N'); DF=flipStruct(DF); N=length(DF.cytoband); DF=trimStruct(DF,1:(N-1));
DF.chr=str2double(cellfun(@(x) x(1),regexp(regexprep(DF.widePeakBoundaries,'chr',''),':','split')));
x=cellfun(@(x) x(2),regexp(DF.widePeakBoundaries,':','split'));
DF.Start=str2double(cellfun(@(x) x(1),regexp(x,'-','split')));
DF.End=str2double(cellfun(@(x) x(2),regexp(x,'-','split')));
DF.gstart=xhg19(DF.chr,DF.Start);
DF.gend=xhg19(DF.chr,DF.End);
DF.cytoband=erase(DF.cytoband,"x");
DF.cytoband=strcat(DF.cytoband,':DEL');
DF.Descriptor=DF.cytoband;



ARMS=load_tsv(arm_coordinates_file);
ARMS.gstart=xhg19(ARMS.chrn,ARMS.start);
ARMS.gend=xhg19(ARMS.chrn,ARMS.xEnd);

%load seg file, trim unnecessarily long names etc.
segs.gstart = xhg19(segs.Chromosome,segs.Start);
segs.gend = xhg19(segs.Chromosome,segs.End);
segs.length = segs.End - segs.Start;

%remove segments overlapping with blacklisted regions
cnv_blacklist = load_tsv(cnv_blacklist_file);
segs = apply_cnv_blacklist(segs,cnv_blacklist,AF,DF,ARMS);
segs.length = segs.End - segs.Start;

%correct segment_means for each sample by sample median
old_seg_means = segs.Segment_Mean;
for samp = unique(segs.Sample)'
    six = find(ismember(segs.Sample,samp));
    sampseg = trimStruct(segs,six);
    sample_median = calc_region_median(sampseg,min(ARMS.gstart),max(ARMS.gend),2);
    segs.Segment_Mean(six) = segs.Segment_Mean(six) - sample_median;
end


segs.log_segment_mean = segs.Segment_Mean;
%segs.Segment_Mean = 2.^(segs.Segment_Mean+1)-2;


%generate arm level calls first so we can check focals against them
ARMS=trimStruct(ARMS,ismember(ARMS.GISTICARM,'Y'));

A = [];
A.pair_id=unique(segs.Sample);
A.arm=ARMS.arm';
A1=A;
j=0;
check=0;
for pair=A.pair_id'
    j=j+1;
    %A.(pair{1}) = zeros(size(A.arm));
    for i=1:length(ARMS.arm)
        arm = ARMS.arm{i};
  
        %select segments that overlap arm, don't select by segment mean yet
        segix = ((segs.gstart <= ARMS.x1(i) & segs.gend >= ARMS.x1(i)) | (segs.gstart < ARMS.x2(i) & segs.gend > ARMS.x2(i)) | (segs.gstart >= ARMS.x1(i)& segs.gend <= ARMS.x2(i)  )) & ismember(segs.Sample,pair); %& segs.Segment_Mean <= single_del_threshold;
        seg1 = trimStruct(segs,segix);
        %skip event if patient has no supporting segs
        if length(seg1.Chromosome) == 0 | sum(seg1.length) < (ARMS.x2(i) - ARMS.x1(i))/PAR.arm_length_fraction_threshold
            check=check+1;
            continue
        end
        region_median = calc_region_median(seg1,ARMS.x1(i),ARMS.x2(i),2);
        A1.(['median_' arm ])(j,1) = region_median;           
        A.(['del_' arm ])(j,1) = 0;
        if region_median < PAR.single_del_threshold & region_median >= PAR.double_del_threshold
            A.(['del_' arm ])(j,1) = 1;
        elseif region_median < PAR.double_del_threshold
            A.(['del_' arm ])(j,1) = 2;
        end
        
        A.(['gain_' arm ])(j,1) = 0;
        if region_median > PAR.single_amp_threshold & region_median <= PAR.double_amp_threshold;
            A.(['gain_' arm ])(j,1) = 1;
        elseif region_median > PAR.double_amp_threshold;
            A.(['gain_' arm ])(j,1) = 2;
        end

        
    end
end


printStruct(A1,-1,'outputs/Arm_median_intensity.tsv');
printStruct(A,-1,'outputs/Arm_presence.tsv');


 
F=[];
F.pair_id=unique(segs.Sample);
F.arm=ARMS.arm';
%AF.cytoband=erase(AF.cytoband,"x");
%AF.cytoband=strcat(AF.cytoband,':AMP');
%AF.Descriptor=AF.cytoband;
F.gene = AF.cytoband;

%F.gene = AF.Descriptor;
%F.cohort = AF.cohort;

F1=F;
j=0;
check=0;
check3=0;
for pair=F.pair_id'
    j=j+1;
    %A.(pair{1}) = zeros(size(A.arm));
    for alix=1:length(AF.Descriptor)
        peak = AF.Descriptor{alix};
        if ismember(peak(2),{'q';'p'})
            arm = peak(1:2);
        else
            arm = peak(1:3);
        end
 
        %arm = AL.arm{i};
  
        %select segments that overlap arm, don't select by segment mean yet
        %segix = ((segs.gstart <= ARMS.x1(i) & segs.gend >= ARMS.x1(i)) | (segs.gstart < ARMS.x2(i) & segs.gend > ARMS.x2(i)) | (segs.gstart >= ARMS.x1(i)& segs.gend <= ARMS.x2(i)  )) & ismember(segs.Sample,pair); %& segs.Segment_Mean <= single_del_threshold;
        %seg1 = trimStruct(segs,segix);
        
        %Select segments overlapping with peak
        segix = ((segs.gstart < AF.gstart(alix) & segs.gend > AF.gstart(alix)) | (segs.gstart < AF.gend(alix) & segs.gend > AF.gend(alix)) | (segs.gstart > AF.gstart(alix) & segs.gend < AF.gend(alix))) & ismember(segs.Sample,pair);
        seg1 = trimStruct(segs,segix);
        %only select segments fully enclosed by the peak
        segix2 = (segs.gstart > AF.gstart(alix) & segs.gend < AF.gend(alix)) & ismember(segs.Sample,pair);
        seg2 = trimStruct(segs,segix2);
        
        %length(seg1.gstart)
        if length(seg1.gstart) ==1
            check=check+1;
        end
        if length(seg1.gstart) > 0
        %Subtract out arm level events ONLY IF PEAK OCCURS IN SIGNIFICANTLY
        %COPY-VARIED ARM
        % replace BA -> ARMS
        if (ismember(arm,ARMS.arm) & strfindk(F.gene(alix),'AMP'))
            check3=check3+1;
            armix = find(ismember(ARMS.arm,arm));
            bound1 = ARMS.gstart(armix);
            bound2 = ARMS.gend(armix);
            %arm_median = calc_region_median(seg1,bound1,bound2,2);
            arm_median = A1.(['median_' arm ])(j);
            seg1.Segment_Mean = seg1.Segment_Mean - arm_median;
        end
        if length(seg1.Chromosome) == 0
            continue
        end
        if strfindk(F.gene(alix),'AMP')
            seg1.Segment_Mean = -seg1.Segment_Mean;
        end
        region_median = calc_region_median(seg1,AF.gstart(alix),AF.gend(alix),7);
        %correct median back to the sign it originally had
        if strfindk(F.gene(alix),'AMP')
            region_median = -region_median;
        end
        end
        
        if length(seg1.gstart) == 0
            region_median='no_coverage';
            j
            pair
            peak
        end
        %create a valid field for the event (amp/del).
        key1 = extractBefore(F.gene{alix},":");
        %key2 = extractAfter(F.gene{alix},":");
        %key3 = extractBefore(key2, ":");
        
        if region_median ~='no_coverage'
        F1.(['median_' key1])(j,1) = region_median;           
        F.(['del_' key1])(j,1) = 0;
        if region_median < PAR.single_del_threshold & region_median >= PAR.double_del_threshold
            F.(['del_' key1])(j,1) = 1;
        elseif region_median < PAR.double_del_threshold
            F.(['del_' key1])(j,1) = 2;
        end
        %add another or condition (if a segment with (#probe>10, segment_mean>0.1) exists, also set it to be 1)
        if length(seg2.Chromosome) > 0
            %if length>=10000b.
            segix3 = (seg2.length >= 10000) & (seg2.Segment_Mean < PAR.single_del_threshold);
        if (sum(segix3) > 0) & (F.(['del_' key1])(j,1) ==0)
            F.(['del_' key1])(j,1) = 1;
        end
        end
        F.(['gain_' key1])(j,1) = 0;
        if region_median > PAR.single_amp_threshold & region_median <= PAR.double_amp_threshold;
            F.(['gain_' key1])(j,1) = 1;
        elseif region_median > PAR.double_amp_threshold;
            F.(['gain_' key1])(j,1) = 2;
        end
        %add another or condition (if a segment with (#probe>10, segment_mean>0.1) exists, also set it to be 1)
        if length(seg2.Chromosome) > 0
            segix3 = (seg2.length >= 10000) & (seg2.Segment_Mean > PAR.single_amp_threshold);
        if (sum(segix3) > 0) & (F.(['gain_' key1])(j,1) ==0)
            F.(['gain_' key1])(j,1) = 1;
        end
        end
        %if (strfindk(F.gene(alix),'AMP') & region_median > PAR.double_amp_threshold) | (strfindk(F.gene(alix),'DEL') & region_median < PAR.double_del_threshold)
        %    F.([F.gene{alix}])(j, 1) = 2;
        %elseif (strfindk(F.gene(alix),'AMP') & (region_median <= PAR.double_amp_threshold & region_median > PAR.single_amp_threshold))| (strfindk(F.gene(alix),'DEL') & (region_median >= PAR.double_del_threshold & region_median < PAR.single_del_threshold));
        %    F.([F.gene{alix}])(j, 1) = 1;
        %end
        end
        
        if region_median=='no_coverage'
        F1.(['median_' key1])(j,1) = -999999999999999;           
        F.(['del_' key1])(j,1) = -999999999999999;
        F.(['gain_' key1])(j,1) = -999999999999999;
         %if (strfindk(F.gene(alix),'AMP') & region_median > PAR.double_amp_threshold) | (strfindk(F.gene(alix),'DEL') & region_median < PAR.double_del_threshold)
        %    F.([F.gene{alix}])(j, 1) = 2;
        %elseif (strfindk(F.gene(alix),'AMP') & (region_median <= PAR.double_amp_threshold & region_median > PAR.single_amp_threshold))| (strfindk(F.gene(alix),'DEL') & (region_median >= PAR.double_del_threshold & region_median < PAR.single_del_threshold));
        %    F.([F.gene{alix}])(j, 1) = 1;
        %end
        end
    end
end   



printStruct(F1,-1,'outputs/Focal_amp_corrected_intensity.tsv');
printStruct(F,-1,'outputs/Focal_amp_presence.tsv');


F=[];
F.pair_id=unique(segs.Sample);
F.arm=ARMS.arm';
F.gene = DF.cytoband;


F1=F;
j=0;
check=0;
check2=0;
check3=0;
for pair=F.pair_id'
    j=j+1;
    %A.(pair{1}) = zeros(size(A.arm));
    for alix=1:length(DF.Descriptor)
        peak = DF.Descriptor{alix};
        if ismember(peak(2),{'q';'p'})
            arm = peak(1:2);
        else
            arm = peak(1:3);
        end
 
        %arm = AL.arm{i};
  
        %select segments that overlap arm, don't select by segment mean yet
        %segix = ((segs.gstart <= ARMS.x1(i) & segs.gend >= ARMS.x1(i)) | (segs.gstart < ARMS.x2(i) & segs.gend > ARMS.x2(i)) | (segs.gstart >= ARMS.x1(i)& segs.gend <= ARMS.x2(i)  )) & ismember(segs.Sample,pair); %& segs.Segment_Mean <= single_del_threshold;
        %seg1 = trimStruct(segs,segix);
        
        %Select segments overlapping with peak
        segix = ((segs.gstart < DF.gstart(alix) & segs.gend > DF.gstart(alix)) | (segs.gstart < DF.gend(alix) & segs.gend > DF.gend(alix)) | (segs.gstart > DF.gstart(alix) & segs.gend < DF.gend(alix))) & ismember(segs.Sample,pair);
        seg1 = trimStruct(segs,segix);
        %only select segments fully enclosed by the peak
        segix2 = (segs.gstart > DF.gstart(alix) & segs.gend < DF.gend(alix)) & ismember(segs.Sample,pair);
        seg2 = trimStruct(segs,segix2);

        if length(seg1.gstart) ==1
            check=check+1;
        end
        if length(seg1.gstart) > 0
        %Subtract out arm level events ONLY IF PEAK OCCURS IN SIGNIFICANTLY
        %COPY-VARIED ARM
        % replace BA -> ARMS
        if (ismember(arm,ARMS.arm) & strfindk(F.gene(alix),'DEL'))
            check3=check3+1;
            armix = find(ismember(ARMS.arm,arm));
            bound1 = ARMS.gstart(armix);
            bound2 = ARMS.gend(armix);
            %arm_median = calc_region_median(seg1,bound1,bound2,2);
            arm_median = A1.(['median_' arm ])(j);
            seg1.Segment_Mean = seg1.Segment_Mean - arm_median;
            if seg1.Segment_Mean ==0
                check2=check2+1
            end
            %if peak=='17q11_2'
            %    segix
            %end
        end
        if length(seg1.Chromosome) == 0
            continue
        end
        if strfindk(F.gene(alix),'AMP')
            seg1.Segment_Mean = -seg1.Segment_Mean;
        end
        %restrict length only within the peak.
        %subtract those that lay outside of the current peak.
        %seg1.gstart
        region_median = calc_region_median(seg1,DF.gstart(alix),DF.gend(alix),7);
        %correct median back to the sign it originally had
        if strfindk(F.gene(alix),'AMP')
            region_median = -region_median;
        end
        end
        
        if length(seg1.gstart) == 0
            region_median='no_coverage';
            j
            pair
            peak
        end
        %create a valid field for the event (amp/del).
        key1 = extractBefore(F.gene{alix},":");
        %key2 = extractAfter(F.gene{alix},":");
        %key3 = extractBefore(key2, ":");
        
        if region_median ~='no_coverage'
        F1.(['median_' key1])(j,1) = region_median;           
        F.(['del_' key1])(j,1) = 0;
        if region_median < PAR.single_del_threshold & region_median >= PAR.double_del_threshold
            F.(['del_' key1])(j,1) = 1;
        elseif region_median < PAR.double_del_threshold
            F.(['del_' key1])(j,1) = 2;
        end
        %add another or condition (if a segment with (#probe>10, segment_mean>0.1) exists, also set it to be 1)
        if length(seg2.Chromosome) > 0
            segix3 = (seg2.length >= 10000) & (seg2.Segment_Mean < PAR.single_del_threshold);
        if (sum(segix3) > 0) & (F.(['del_' key1])(j,1) ==0)
            F.(['del_' key1])(j,1) = 1;
        end
        end
        F.(['gain_' key1])(j,1) = 0;
        if region_median > PAR.single_amp_threshold & region_median <= PAR.double_amp_threshold
            F.(['gain_' key1])(j,1) = 1;
        elseif region_median > PAR.double_amp_threshold
            F.(['gain_' key1])(j,1) = 2;
        end
        %add another or condition (if a segment with (#probe>10, segment_mean>0.1) exists, also set it to be 1)
        if length(seg2.Chromosome) > 0       
            segix3 = (seg2.length >= 10000) & (seg2.Segment_Mean > PAR.single_amp_threshold);
        if (sum(segix3) > 0) & (F.(['gain_' key1])(j,1) ==0)
            F.(['gain_' key1])(j,1) = 1;
        end
        end
        %if (strfindk(F.gene(alix),'AMP') & region_median > PAR.double_amp_threshold) | (strfindk(F.gene(alix),'DEL') & region_median < PAR.double_del_threshold)
        %    F.([F.gene{alix}])(j, 1) = 2;
        %elseif (strfindk(F.gene(alix),'AMP') & (region_median <= PAR.double_amp_threshold & region_median > PAR.single_amp_threshold))| (strfindk(F.gene(alix),'DEL') & (region_median >= PAR.double_del_threshold & region_median < PAR.single_del_threshold));
        %    F.([F.gene{alix}])(j, 1) = 1;
        %end
        end
        
        if region_median=='no_coverage'
        F1.(['median_' key1])(j,1) = -999999999999999;           
        F.(['del_' key1])(j,1) = -999999999999999;
        F.(['gain_' key1])(j,1) = -999999999999999;
         %if (strfindk(F.gene(alix),'AMP') & region_median > PAR.double_amp_threshold) | (strfindk(F.gene(alix),'DEL') & region_median < PAR.double_del_threshold)
        %    F.([F.gene{alix}])(j, 1) = 2;
        %elseif (strfindk(F.gene(alix),'AMP') & (region_median <= PAR.double_amp_threshold & region_median > PAR.single_amp_threshold))| (strfindk(F.gene(alix),'DEL') & (region_median >= PAR.double_del_threshold & region_median < PAR.single_del_threshold));
        %    F.([F.gene{alix}])(j, 1) = 1;
        %end
        end

        

    end
end   



printStruct(F1,-1,'outputs/Focal_del_corrected_intensity.tsv');
printStruct(F,-1,'outputs/Focal_del_presence.tsv');

end


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
            if (segs.gstart(segix(j)) <= cnv_blacklist.gstart(i) & segs.gend(segix(j)) >= cnv_blacklist.gend(i)) 
                cnv_length = cnv_blacklist.gend(i) - cnv_blacklist.gstart(i);
                seg_length = segs.length(segix(j));
                if (cnv_length/seg_length > 0.8)
                    segs.remove_seg(segix(j)) = 1;
                end
            %if segments spans just start of blacklisted region, set end of
            %segment to start of blacklisted region
            elseif (segs.gstart(segix(j)) <= cnv_blacklist.gstart(i) & segs.gend(segix(j)) >= cnv_blacklist.gstart(i))
                    segs.gend(segix(j))   = cnv_blacklist.gstart(i);
                    segs.End(segix(j)) = cnv_blacklist.Start(i);
            %if segment spans just end of blacklisted region, set start to
            %end of blacklisted region
            elseif (segs.gstart(segix(j))   <= cnv_blacklist.gend(i) & segs.gend(segix(j)) >= cnv_blacklist.gend(i)) 
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



