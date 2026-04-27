function S = CreateSignature(Name,Description,Type,Identifiers,IDs,Coeff)

if isempty(Name)
    S.Name = "";
else
    S.Name = ConvertStr(Name,'string');
end

if isempty(Description)
    S.Description = "";
else
    S.Description = ConvertStr(Description,'string');
end

if isempty(Type)
    S.Type = "";
else
    S.Type = ConvertStr(Type,'string');
end

if isempty(Identifiers)
    S.Identifiers = "";
else
    S.Identifiers = ConvertStr(Identifiers,'string');
end

if isempty(IDs)
    S.IDs = "";
else
    S.IDs = ConvertStr(IDs,'string');
end

if isempty(Name)
    S.Coeff = zeros(0,0);
else
    S.Coeff = Coeff;
end

