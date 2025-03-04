function [RESULTS] = GetCellosariousID(IdList,CelloDATA,ManualCorrection)

HumanOnly = true;

CleanChar = {' ','-', '_'};
ReplaceChar = "";

n = length(IdList);

PrimaryId = [CelloDATA.CellLine(:).ID]';
PrimaryIdClean = replace(PrimaryId,CleanChar,ReplaceChar);

RESULTS = strings(n,5);

% Start checking

if HumanOnly
    indx_human = contains( [CelloDATA.CellLine(:).OX]','NCBI_TaxID=9606');
else
    indx_human = true(CelloDATA.nCellLines,1);
end
 EntryNumHuman = find(indx_human);


for i = 1:n
    Id = IdList(i);
    Id_Clean = replace(Id,CleanChar,ReplaceChar); 
    RESULTS(i,1) = Id;
    % Start to check against Primary Id, case sensative
    
    if ~isempty(ManualCorrection)
        indx_MC = strcmp(Id,ManualCorrection(:,1));
        if any(indx_MC)
            Id = ManualCorrection(indx_MC,2);
        end
    end
    indx = strcmp(Id,PrimaryId) & indx_human;
   

    if any(indx)
        RESULTS(i,2) = PrimaryId(indx);
        RESULTS(i,3) = CelloDATA.CellLine(indx).AC;
        RESULTS(i,4) = CelloDATA.CellLine(indx).OX;
        RESULTS(i,5) = "Primary Id match";
    else

        % Check agains Synonyms, case sensative
        indx = strcmp(Id,CelloDATA.Synonyms);
        EntryNum = CelloDATA.SynonymsIndx(indx);
        if (sum(indx) == 1 || isscalar(unique((CelloDATA.SynonymsIndx(indx))))) && ismember(EntryNum,EntryNumHuman)% Only one unique Hit
            
            RESULTS(i,2) = CelloDATA.CellLine(EntryNum).ID;
            RESULTS(i,3) = CelloDATA.CellLine(EntryNum).AC;
            RESULTS(i,4) = CelloDATA.CellLine(EntryNum).OX;
            RESULTS(i,5) = "Unique match to synonym";

        elseif sum(indx) > 1 % Multiple hits
            % Check if only one is human
            OX = [CelloDATA.CellLine(EntryNum).OX]';
            indx_hum = strcmp('NCBI_TaxID=9606; ! Homo sapiens (Human)',OX);
            if sum(indx_hum) == 1
                RESULTS(i,2) = CelloDATA.CellLine(EntryNum(indx_hum)).ID;
                RESULTS(i,3) = CelloDATA.CellLine(EntryNum(indx_hum)).AC;
                RESULTS(i,4) = CelloDATA.CellLine(EntryNum(indx_hum)).OX;
                RESULTS(i,5) = "Multiple match to synonym, only one human";

            else
                % Check if ID works with a cleaned Primary ID
                ID_clean_tmp = [CelloDATA.CellLine(EntryNum).ID]';
                ID_clean_tmp = replace(ID_clean_tmp,CleanChar,ReplaceChar);
                indx_id = strcmp(Id_Clean,ID_clean_tmp);
                if sum(indx_id) == 1
                    RESULTS(i,2) = CelloDATA.CellLine(EntryNum(indx_id)).ID;
                    RESULTS(i,3) = CelloDATA.CellLine(EntryNum(indx_id)).AC;
                    RESULTS(i,4) = CelloDATA.CellLine(EntryNum(indx_id)).OX;
                    RESULTS(i,5) = "Multiple match to synonym, one match cleaned primary Id";
                else
                    RESULTS(i,5) = "MULTIPLE MATCHES!";

                end

            end
        else
            indx = strcmpi(Id,PrimaryId) & indx_human;
            if sum(indx) == 1 % Only one unique Hit
                RESULTS(i,2) = PrimaryId(indx);
                RESULTS(i,3) = CelloDATA.CellLine(indx).AC;
                RESULTS(i,4) = CelloDATA.CellLine(indx).OX;
                RESULTS(i,5) = "Unique match to Primary Id, ignoring case";
            else

                indx = strcmpi(Id,CelloDATA.Synonyms); % Case insensative
                if sum(indx) == 1 || isscalar(unique((CelloDATA.SynonymsIndx(indx))))% Only one unique Hit
                    EntryNum = CelloDATA.SynonymsIndx(indx);
                    RESULTS(i,2) = CelloDATA.CellLine(EntryNum).ID;
                    RESULTS(i,3) = CelloDATA.CellLine(EntryNum).AC;
                    RESULTS(i,4) = CelloDATA.CellLine(EntryNum).OX;
                    RESULTS(i,5) = "Unique match to synonym, ignoring case";
                else
                    indx = strcmpi(Id_Clean,PrimaryIdClean) & indx_human ;

                    if sum(indx) == 1 % Only one unique Hit
                        RESULTS(i,2) = PrimaryId(indx);
                        RESULTS(i,3) = CelloDATA.CellLine(indx).AC;
                        RESULTS(i,4) = CelloDATA.CellLine(indx).OX;
                        RESULTS(i,5) = "Primary Id matched with cleaned Id";
                    else

                        RESULTS(i,5) = "NO MATCH FOUND";
                    end


                end
            end
        end


    end
end
