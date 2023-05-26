function [DATA_out, nVal, rVal, RangeVal] = EPICv2_Compress(DATA)

id_orig = DATA.ColId;
id_simple = extractBefore(id_orig,"_");
id_simple_multiple = id_simple;
[id_simple_unique,ia,ic] = unique(id_simple,'stable');
id_simple_multiple(ia) = [];
counter = 0;

nVal = zeros(size(id_simple_multiple));

rVal = ones(size(id_orig)) * NaN;
RangeVal = ones(size(id_orig)) * NaN;


DATA_out = DATA;
DATA_out.ColId = id_simple_unique;
DATA_out.nCol = length(id_simple_unique);
DATA_out.X = DATA.X(:,ia);

for i=1:length(id_simple_multiple)
   
    indx = find(strcmp(id_simple_multiple(i),id_simple));
    nVal(i) = length(indx);
    x_mean = mean(DATA.X(:,indx),2,'omitnan');

    indx_x = strcmp(id_simple_multiple(i),id_simple_unique);
    DATA_out.X(:,indx_x) = x_mean;
    for j=1:length(indx)-1
        for k=2:length(indx)
            counter = counter + 1;
            rVal(counter) = corr( DATA.X(:,indx(j)),DATA.X(:,indx(k)),'Type','Spearman','Rows','pairwise');
            RangeVal(counter) = range(DATA.X(:,indx([j,k])),'all');

        end
    end

end

rVal = rVal(1:counter);
RangeVal = RangeVal(1:counter);
