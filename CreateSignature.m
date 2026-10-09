function S = CreateSignature(Name, Description, Source, Type, Identifiers, IDs, Coeff)
%CREATESIGNATURE Construct one signature definition with file provenance.
%   S = CreateSignature(Name, Description, Source, Type, Identifiers, IDs, Coeff)
%
%   Name         Signature name, such as a worksheet name or gene-set name.
%   Description  Human-readable signature description and reference metadata.
%   Source       Originating file or another caller-supplied provenance label.
%                CreateSignaturesFromExcel supplies the filename + extension.
%   Type         Scoring-method label, such as "PCA", "linear model", "mean".
%                It is NOT the workbook layout or an assay label.
%   Identifiers  Labels describing the columns of IDs.
%   IDs          Member identifiers/annotations, one row per signature member.
%   Coeff        Member coefficients. This function does not calculate scores.
%
%   This is the supplied seven-argument constructor with one functional fix:
%       S.Source = ConvertStr(Source, 'string');
%   The uploaded version incorrectly converted Type in that assignment.
%   All other conversion and empty-input behavior is intentionally retained
%   for compatibility with code that already calls this constructor.
%
%   Empty IDs historically become scalar "" here. The Excel importer restores
%   0-by-N IDs / 0-by-1 Coeff for its own empty-member signatures. This avoids
%   changing the behavior of other code using this constructor directly.
%
%   Dependency: the existing ConvertStr.m utility must be on the MATLAB path.
%   See also CreateSignaturesFromExcel.

% Convert textual metadata with the same project utility used previously.
if isempty(Name)
    S.Name = "";
else
    S.Name = ConvertStr(Name, 'string');
end

if isempty(Description)
    S.Description = "";
else
    S.Description = ConvertStr(Description, 'string');
end

% File provenance and scoring method are independent fields.
if isempty(Source)
    S.Source = "";
else
    S.Source = ConvertStr(Source, 'string');
end

if isempty(Type)
    S.Type = "";
else
    S.Type = ConvertStr(Type, 'string');
end

if isempty(Identifiers)
    S.Identifiers = "";
else
    S.Identifiers = ConvertStr(Identifiers, 'string');
end

if isempty(IDs)
    S.IDs = "";
else
    S.IDs = ConvertStr(IDs, 'string');
end

% Preserve the original empty-Name template behavior. The importer validates
% and aligns member coefficients before passing them to this constructor.
if isempty(Name)
    S.Coeff = zeros(0, 0);
else
    S.Coeff = Coeff;
end
end
