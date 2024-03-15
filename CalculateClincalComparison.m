function CalculateClincalComparison(DATA_A,DATA_B,ClinicalIds,ClinicalType)

x1=cell(DATA_A.nRow,1);
x2=cell(DATA_B.nRow,1);
x1(:)={'A'};
x2(:)={'B'};
Group=[x1;x2];
[ClinVar,ia,ib] = intersect(DATA_A.RowAnnotationFields,DATA_B.RowAnnotationFields,'Stable');
Clinical=[DATA_A.RowAnnotation(:,ia);DATA_B.RowAnnotation(:,ib)];

nA=DATA_A.nRow;
nB=DATA_B.nRow;

fprintf('\n')
fprintf('\tA(n=%u)\tB(n=%u)\tp-value\n',nA,nB)


for i=1:numel(ClinicalIds)
    indx = strcmpi(ClinicalIds{i},ClinVar);
    if any(indx)
        switch lower(ClinicalType{i})

            case {'continues','continuous','cont'}
                x_A=DATA_A.RowAnnotation(:,indx);
                x_A=strrep(x_A,'[Not Available]','NaN');
                x_A=cellfun(@(x) str2num(x), x_A, 'UniformOutput', 0);
                x_A=cell2mat(x_A);
                x_B=DATA_B.RowAnnotation(:,indx);
                x_B=strrep(x_B,'[Not Available]','NaN');
                x_B=cellfun(@(x) str2num(x), x_B, 'UniformOutput', 0);
                x_B=cell2mat(x_B);
                [~,p_tt] = ttest2(x_A,x_B,0.05,'both','unequal');
                if p_tt < 0.0001
                    p_txt = sprintf('<0.0001');
                else
                    p_txt = sprintf('%.4f',p_tt);
                end
                fprintf('%s\t%.1f ± %.1f\t%.1f ± %.1f\t%s\n',ClinicalIds{i},mean(x_A,"omitnan"),std(x_A,"omitnan"),mean(x_B,"omitnan"),std(x_B,"omitnan"),p_txt)
                fprintf('\n')
            case {'categorical','cat'}

                fprintf('%s\n',ClinicalIds{i})
                [tbl,chi2,p,labels] = crosstab(Clinical(:,indx),Group);
                [lab,sort_indx ]= sort(labels(:,1));
                tbl=tbl(sort_indx,:);
                if length(sort_indx) == 2
                    [~,p,~] = fishertest(tbl);
                end
                if p < 0.0001
                    p_txt = sprintf('<0.0001');
                else
                    p_txt = sprintf('%.4f',p);
                end
                fprintf('%s\t%u\t%u\t%s\n',lab{1},tbl(1,1),tbl(1,2),p_txt)
                for i=2:length(sort_indx)
                    fprintf('%s\t%u\t%u\n',lab{i},tbl(i,1),tbl(i,2))
                end
                fprintf('\n')
        end


    end


end
