function c=chrom2num(chr1)
if (iscell(chr1))
    chr=regexprep(chr1,'chr','');
    chr=regexprep(chr,'X','23');
    chr=regexprep(chr,'Y','24');
    chr=regexprep(chr,'MT','25');
    k=find(cellfun(@length,regexp(chr,'[a-c|h-z|A-C|H-Z]','start'))>0);
    chr(k)={'0'};
    c=str2num(char(chr));
else
    k=strmatch('chr',chr1);
    if (~isempty(k))
        for i=k(:)'
            chr(i,:)= regexprep(chr1(i,:), 'chr', '');
        end
    else
        chr=char(chr1);
    end

    n=size(chr);
    if (n(2)<2)
        chr=[chr repmat(' ',n(1),1)];
    end

    k=strmatch('X',chr);
    if (~isempty(k)), chr(k,:)=repmat('23',size(k),1); end
    k=strmatch('Y',chr);
    if (~isempty(k)), chr(k,:)=repmat('24',size(k),1); end
    k=strmatch('MT',chr);
    if (~isempty(k)), chr(k,:)=repmat('25',size(k),1); end
    k=strmatch('N',chr);
    if (~isempty(k)), chr(k,:)=repmat('00',size(k),1); end
    c=str2num(chr);
end
