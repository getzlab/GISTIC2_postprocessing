function X = load_tsv(tsv_file,comment,delimiter,first_field)
% load_tsv(tsv_file,tsv_file,comment,delimiter,first_field) skips header lines 
if nargin<2
    comment='#';
end
if nargin<3
    delimiter='\t';
end

cmd =['grep "^' comment '" ' '"' tsv_file '"' ' | wc -l'];
[o nh]=unix(cmd);
nh=round(abs(str2num(nh)));

% if nargin==1
%     cmd =['grep -n "^\w" ' tsv_file ' | cut -d":" -f1'];
%     [o nh]=unix(cmd);
%     n=regexp(nh,'\n','split');
%     nh=round(abs(str2num(n{1})))-1;    
% end

if nargin>3
    cmd =['grep -n "^' first_field '" "' tsv_file '"' ];
    [o q]=unix(cmd);
    nh=str2double(cellfun(@(x) x(1),regexp({q},':','split')))-1
end

if isempty(nh), nh=0;end
X0 = readtable(tsv_file,'Delimiter',delimiter,'FileType','text','HeaderLines',nh,'ReadVariableNames',1);
X = table2struct(X0,'ToScalar',true);

