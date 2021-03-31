function [X1]=flipStruct(X)
if nargin<1
    return
end
ff=fieldnames(X);
row=X.(ff{1});
rowname=matlab.lang.makeValidName(row);

X1.(ff{1})=ff(2:end);
for i=1:length(row)
    q=X.(ff{2})(i);
    if iscell(q)
      q1={''};
    else
      q1=NaN;
    end
    X1.(rowname{i})=repmat(q1,size(X1.(ff{1})));
end
for i=1:length(X1.(ff{1}))
    for j=1:length(row)
        X1.(rowname{j})(i)=X.(ff{i+1})(j);
    end
end
