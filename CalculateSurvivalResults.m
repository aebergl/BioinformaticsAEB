function RESULTS_DATA = CalculateSurvivalResults(DATA,Group1,GroupId,Surv_Type,DataType,varargin)


MinNumSampleSize = 10;

if isempty(DataType)
    DataType = 'Beta-value';
end

if ~iscell(Group1)
    Group1 = {Group1};
end


if ~isempty(GroupId)
    GroupIdColumn = (strcmp(GroupId,DATA.RowAnnotationFields));
    GroupIdColumn = GroupIdColumn & (cumsum(GroupIdColumn) == 1); % Picks the first occurence

    SampleIndxMatrix = false(DATA.nRow,numel(Group1));
    for i=1:numel(Group1)
        SampleIndxMatrix(:,i) = strcmp(Group1{i},DATA.RowAnnotation(:,GroupIdColumn));
    end
    SampleIndx_x = any(SampleIndxMatrix,2);
    SampleId_x = DATA.RowId(SampleIndx_x);
else
    SampleIndx_x = true(DATA.nRow,1);
    SampleId_x = DATA.RowId;
    GroupId = 'All samples';
    Group1 = 'All samples';
end

RESULTS_DATA = CreateDataStructure(DATA.nCol,8,[],[]);

% Add Info
RESULTS_DATA.Title = 'Survival analysis results';
RESULTS_DATA.Info.Source = inputname(1);
RESULTS_DATA.Info.OS_Type = Surv_Type;
RESULTS_DATA.Info.Platform = DataType;
RESULTS_DATA.Info.GroupID = GroupId;
RESULTS_DATA.Info.Group1 = Group1;
RESULTS_DATA.Info.MinNumSampleSize = MinNumSampleSize;
RESULTS_DATA.Info.SampleId_x = SampleId_x;
RESULTS_DATA.Info.numSamples_x = length(SampleId_x);
RESULTS_DATA.RowId = DATA.ColId;
RESULTS_DATA.RowAnnotation = DATA.ColAnnotation;
RESULTS_DATA.RowAnnotationFields = DATA.ColAnnotationFields;

VarNames = {'p logrank','q logrank','fdr logrank','HR logrank','p coxreg','q coxreg','fdr coxreg','HR coxreg'}';
VarNames = strcat(VarNames,{' '},Surv_Type);
VarNames = [VarNames; {'Range','Num Values'}'];
RESULTS_DATA.ColId=VarNames;

RESULTS_DATA.ColAnnotationFields = {'VarId'};
RESULTS_DATA.ColAnnotation = RESULTS_DATA.ColId;

X_x = DATA.X(SampleIndx_x,:);
if strcmpi('M-value',DataType)
    X_x = B2M(X_x);
end

SurvColId = strcmp(Surv_Type,DATA.SURVIVAL.SurvivalTypes);
TimeVar = DATA.SURVIVAL.SurvTime(SampleIndx_x,SurvColId);
EventVar = DATA.SURVIVAL.SurvEvent(SampleIndx_x,SurvColId);

% Check time variable for missing data and timepoints < TimeMin
rem_indx_time = ( isnan(TimeVar) | (TimeVar < 0) );

% Check Event variable for missing data and/or empty cells and cells with NA
if isnumeric(EventVar) || islogical(EventVar)
    rem_indx_event = isnan(EventVar);
elseif iscell(EventVar)
    rem_indx_event = (cellfun('isempty',EventVar) | strcmpi('[Not Available]',EventVar) | strcmpi('NA',EventVar));
end
rem_indx = (rem_indx_time | rem_indx_event);
TimeVar(rem_indx) = [];
EventVar(rem_indx) = [];
X_x(rem_indx,:) = [];

EventVarBin = zeros(size(EventVar));

indx_Event = strcmpi('Dead',EventVar) | strcmpi('Deceased',EventVar) | strcmpi('Relapsed',EventVar)...
    |  strcmpi('Yes',EventVar) | strcmpi('Event',EventVar) | strcmpi('Progression',EventVar)...
    | strcmpi('Progressed',EventVar) | strcmpi('TRUE',EventVar);

indx_NoEvent = strcmpi('Alive',EventVar) | strcmpi('Living',EventVar) | strcmpi('NotRelapsed',EventVar)...
    | strcmpi('DiseaseFree',EventVar) | strcmpi('No',EventVar) | strcmpi('Censored',EventVar)...
    | strcmpi('NoProgression',EventVar) | strcmpi('NoEvent',EventVar) | strcmpi('FALSE',EventVar);
EventVarBin(indx_Event) = 1;
EventVarBin(indx_NoEvent) = 0;

coxphopt = statset('coxphfit');
coxphopt.MaxFunEvals = 500;
coxphopt.MaxIter = 500;
coxphopt.TolBnd = 1.0000e-08;

p_logrank = ones(DATA.nCol,1) * NaN;
HR_logrank = ones(DATA.nCol,1) * NaN;
p_coxreg = ones(DATA.nCol,1) * NaN;
HR_coxreg = ones(DATA.nCol,1) * NaN;
RangeVal  = ones(DATA.nCol,1) * NaN;
nVal = ones(DATA.nCol,1) * NaN;


if DATA.nCol > 100
    parfor i=1:DATA.nCol
        x = X_x(:,i);
        nVal(i) = sum(~isnan(x));
        x_tmp = x(isfinite(x));
        [~,n]=GroupCount(x_tmp);

        if nVal(i) > MinNumSampleSize && nVal(i) - max(n) + 1 > MinNumSampleSize
            RangeVal(i) = range(x);
            try
                [p_logrank(i),~,stats]= MatSurv(TimeVar,EventVar,x,'CutPoint','median','NoPlot',true,'Print',false,'NoWarnings',true);
                HR_logrank(i) = stats.HR_logrank;
            catch
            end
            try
                [~,~,~,stats_cox] = coxphfit(normalize(x),TimeVar,'censoring',~EventVarBin,'Options',coxphopt,'Baseline',0);
                p_coxreg(i)=stats_cox.p;
                HR_coxreg(i) = exp(stats_cox.beta);
            catch
            end
        end
    end
else
    for i=1:DATA.nCol
        x = X_x(:,i);
        nVal(i) = sum(~isnan(x));
        x_tmp = x(isfinite(x));
        [~,n]=GroupCount(x_tmp);

        if nVal(i) > MinNumSampleSize && nVal(i) - max(n) + 1 > MinNumSampleSize
            RangeVal(i) = range(x);
            try
                [p_logrank(i),~,stats]= MatSurv(TimeVar,EventVar,x,'CutPoint','median','NoPlot',true,'Print',false,'NoWarnings',true);
                HR_logrank(i) = stats.HR_logrank;
            catch
            end
            try
                [~,~,~,stats_cox] = coxphfit(normalize(x),TimeVar,'censoring',~EventVarBin,'Options',coxphopt,'Baseline',0);
                p_coxreg(i)=stats_cox.p;
                HR_coxreg(i) = exp(stats_cox.beta);
            catch
            end
        end
    end
end
RESULTS_DATA.X(:,1) = p_logrank;
RESULTS_DATA.X(:,4) = HR_logrank;
RESULTS_DATA.X(:,5) = p_coxreg;
RESULTS_DATA.X(:,8) = HR_coxreg;
RESULTS_DATA.X(:,9) = RangeVal;
RESULTS_DATA.X(:,10) = nVal;
try
    [~, RESULTS_DATA.X(:,2),~] = mafdr(p_logrank);
catch
    RESULTS_DATA.X(:,2) = ones(DATA.nCol,1);
end

try
    [RESULTS_DATA.X(:,3)] = mafdr(p_logrank,'BHFDR',true);
catch
    RESULTS_DATA.X(:,3) = ones(DATA.nCol,1);
end


try
    [~, RESULTS_DATA.X(:,6),~] = mafdr(p_coxreg);
catch
    RESULTS_DATA.X(:,6) = ones(DATA.nCol,1);
end

try
    [RESULTS_DATA.X(:,7)] = mafdr(p_coxreg,'BHFDR',true);
catch
    RESULTS_DATA.X(:,7) = ones(DATA.nCol,1);
end
