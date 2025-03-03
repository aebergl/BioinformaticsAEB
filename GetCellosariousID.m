function [RESULTS] = GetCellosariousID(IdList,CelloDATA);
n = length(IdList);

PrimaryId = [CelloDATA.CellLine(:).ID]';

RESULTS = strings(n,5);

% Start checking

for i = 1:n
    Id = IdList(i);
    % Start to check against Primary Id
    indx = strcmp(Id,PrimaryId);
    RESULTS(i,1) = Id;
    if any(indx)
        RESULTS(i,2) = PrimaryId(indx);
        RESULTS(i,3) = CelloDATA.CellLine(indx).AC;
        RESULTS(i,4) = CelloDATA.CellLine(indx).OX;
        RESULTS(i,5) = "Primary Id match";
    else
        % Check agains Synonyms
        indx = strcmp(Id,CelloDATA.Synonyms);
        if sum(indx) == 1 % Only one unique Hit
            EntryNum = CelloDATA.SynonymsIndx(indx);
            RESULTS(i,2) = CelloDATA.CellLine(EntryNum).ID;
            RESULTS(i,3) = CelloDATA.CellLine(EntryNum).AC;
            RESULTS(i,4) = CelloDATA.CellLine(EntryNum).OX;
            RESULTS(i,5) = "Unique match to synonym";
        elseif sum(indx) > 1
            EntryNum = CelloDATA.SynonymsIndx(indx);
            % Check if only one is human
            OX = [CelloDATA.CellLine(EntryNum).OX]';
            indx_hum = (strcmp('NCBI_TaxID=9606; ! Homo sapiens (Human)',OX));
            if sum(indx_hum) == 1
                RESULTS(i,2) = CelloDATA.CellLine(EntryNum(indx_hum)).ID;
                RESULTS(i,3) = CelloDATA.CellLine(EntryNum(indx_hum)).AC;
                RESULTS(i,4) = CelloDATA.CellLine(EntryNum(indx_hum)).OX;
                RESULTS(i,5) = "Multiple match to synonym, only one human";
            else
                % Check if ID works with a cleaned is
                ID_tmp = [CelloDATA.CellLine(EntryNum).ID]';
                ID_tmp = strrep(ID_tmp,'-','');
                ID_tmp = strrep(ID_tmp,' ','');
                indx_id = strcmp(Id,ID_tmp);
                if sum(indx_id) == 1
                    RESULTS(i,2) = CelloDATA.CellLine(EntryNum(indx_id)).ID;
                    RESULTS(i,3) = CelloDATA.CellLine(EntryNum(indx_id)).AC;
                    RESULTS(i,4) = CelloDATA.CellLine(EntryNum(indx_id)).OX;
                    RESULTS(i,5) = "Multiple match to synonym, One matche cleaned primary Id";
                else
                    RESULTS(i,5) = "MULTIPLE MATCHES!";

                end

            end
        else
            indx = strcmpi(Id,CelloDATA.Synonyms);        if sum(indx) == 1 % Only one unique Hit
                EntryNum = CelloDATA.SynonymsIndx(indx);
                RESULTS(i,2) = CelloDATA.CellLine(EntryNum).ID;
                RESULTS(i,3) = CelloDATA.CellLine(EntryNum).AC;
                RESULTS(i,4) = CelloDATA.CellLine(EntryNum).OX;
                RESULTS(i,5) = "----------";
            end

        end


    end


end

