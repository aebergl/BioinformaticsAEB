function fh = Plot_Manhattan(DATA,ChrPos,Chrfield,Y_variable,varargin)
% USAGE:
%   fh = Plot_Manhattan(DATA,ChrPos,Chrfield,Y_variable)
%   Chreates a Manhattan plot across all chromosomes
%
% INPUTS:
% * DATA: DATA structure with results
% * ChrPos: Name of chr position field in RowAnnotationFields 'CpG_beg'
% * Chrfield: Name of Chr field in RowAnnotationFields 'CpG_chrm'
% * Y_variable: Name variable to be used on Y-axis from ColId 'HR coxreg DSS'
%
% OUTPUTS:
% * fh: Figure handle to Chromosome figure
%
%   options ---------------------------------------
%
%   'MarkSelected'  Highligh selected points
%   'Print'         Do not transpose the input file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2020 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sex=0;
FontSize = 7;
LineWidth = 0.5;
GridLines = 'on';
RightMargin = 0.1;
LeftMargin = 0.1;

FigSize = [7 2.5];

Cmap=lines;


% Select and convert base pair position
ChrPosColumn = strcmpi(ChrPos,DATA.RowAnnotationFields);
if ~any(ChrPosColumn)
    error('Given value for ChrPos not found, %s not found in RowAnnotationFields',ChrPos{1})
end
ChrPos = DATA.RowAnnotation(:,ChrPosColumn);
ChrPos = cellfun(@(x) str2num(x), ChrPos, 'UniformOutput', 0);
ChrPos = cell2mat(ChrPos);
ChrPos = double(ChrPos);


% Select Chromosome number
ChrColumn = strcmpi(Chrfield,DATA.RowAnnotationFields);
if ~any(ChrColumn)
    error('Given value for Chrfield not found, %s not found in RowAnnotationFields',ChrColumn{1})
end
Chr = DATA.RowAnnotation(:,ChrColumn);

if any(strcmp('chr1',Chr))
    ChrList = {'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM','NA'};
    Chr = categorical(Chr,ChrList);
    Chr = double(Chr);
else


end


indx_Y_Val = strcmpi(Y_variable,DATA.ColId);
if any(indx_Y_Val)
    Y_Val = DATA.X(:,indx_Y_Val);
else
    error('%s not found',Y_Type)
end


tab2=[Chr,ChrPos,Y_Val];

if sex==0
    tab2=tab2(tab2(:,1)<23,:); %remove chr23 if not wanting sex chromosomes
end

fh=figure('Name','GSEA Plot','Color','w','Tag','GSEA Plot figure','Units','inches');
fh.Position(3:4) = FigSize;
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize,'Units','inches',...
    'PositionConstraint','outerposition','Clipping','off');

ah.LineWidth = LineWidth;
ah.XGrid = GridLines;
ah.YGrid = GridLines;
ah.Colormap = Cmap;
ah.ColorOrder = Cmap;



bptrack=0; %variable to track base pairs, this helps the plotting function to know where to start plotting the next chromosome
tab2(tab2(:,3)==0,:)=[];
for i=1:22+sex
    hold on
    plot(ah,tab2(tab2(:,1)==i,2)+max(bptrack),tab2(tab2(:,1)==i,3),'.'); %a scatterplot. On the x axis, the base pair number + bptrack, which starts at 0. On the y axis, - log 10 p
    bptrack=[bptrack,max(bptrack)+max(max(tab2(tab2(:,1)==i,:)))]; %this updates bptrack by adding the highest base pair number of the chromosome. At the end, we should get the total number of base pairs in the human genome. All values of bptrack are stored in a vector. They're useful later
end

xlim([0,max(bptrack)])

M=movmean(bptrack,2); %this calculates the moving average of bptrack and uses them as chromosome label markers. This puts them in the middle of each chromosome.
xticks(M(2:end));
xlabel('Chromosome')
xticklabels( 1:23 );
