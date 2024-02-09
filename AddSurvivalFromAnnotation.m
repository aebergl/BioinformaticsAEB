function DATA = AddSurvivalFromAnnotation(DATA,EventsColumnId,TimeColumnId,SurvivalTypes,Units)


if ~iscell(EventsColumnId)
    EventsColumnId = {EventsColumnId};
end
[EventColumns]  = ismember(DATA.RowAnnotationFields,EventsColumnId);
if sum(EventColumns) < length(EventsColumnId)
    error('Not all Event column found')
end
Events = DATA.RowAnnotation(:,EventColumns);

if ~iscell(TimeColumnId)
    TimeColumnId = {TimeColumnId};
end
[TimeColumn]  = ismember(DATA.RowAnnotationFields,TimeColumnId);
if sum(TimeColumn) < length(TimeColumnId)
    error('Not all Event column found')
end
Time = DATA.RowAnnotation(:,TimeColumn)


Time=cellfun(@(x) str2num(x), Time, 'UniformOutput', 0)
Time = cell2mat(Time)


DATA.SURVIVAL.RowId = DATA.RowId;
DATA.SURVIVAL.SurvivalTypes = SurvivalTypes;
DATA.SURVIVAL.Units = Units;

DATA.SURVIVAL.SurvEvent = Events;
DATA.SURVIVAL.SurvTime = Time;

