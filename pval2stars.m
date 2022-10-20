function [stat_str,cutoff_str] = pval2stars(pval,cut_off)

if isempty(cut_off) || (cut_off ==0)
    cut_off = [0.05 0.01 0.001 0.0001];
end

cutoff_str = cell(length(cut_off),1);
stat_str = cell(size(pval));
stat_str(:) = {'N.S.'};
star = '';
for i =1:length(cut_off)
    indx = (pval <= cut_off(i));
    if ~any(indx)
        i=i-1;
        break
    end
    star = [star,'*'];
    
    stat_str(indx) =   {star};
    cutoff_str{i} = sprintf('%s p < %g',star,cut_off(i));
end
cutoff_str=cutoff_str(1:i);