function [GeneIdOut, GeneSymbolOut, FullAnnotationOut] = GeneSymbol2GeneId(GeneSymbols,RefseqGeneInfoFile)

if isempty(RefseqGeneInfoFile)
    ftpobj = ftp('ftp.ncbi.nlm.nih.gov');
    cd(ftpobj,'refseq/H_sapiens/');
    mget(ftpobj,'Homo_sapiens.gene_info.gz');
    gunzip('Homo_sapiens.gene_info.gz');
    [FidInputFile,message] = fopen('Homo_sapiens.gene_info','r');
    DATA.Info.GeneAnnotation = 'ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info';
    DATA.Info.GeneAnnotationData = date;
    if  FidInputFile == -1
        disp(InputFile)
        disp(message)
        return
    end
else
    [FidInputFile,message] = fopen(RefseqGeneInfoFile,'r');
    DATA.Info.GeneAnnotation = 'file_name';
    DATA.Info.GeneAnnotationData = date;
    
    if  FidInputFile == -1
        disp(RefseqGeneInfoFile)
        disp(message)
        return
    end
    
end

%These are the column iDs
% 1     {'#tax_id'                              }
% 2     {'GeneID'                               }
% 3     {'Symbol'                               }
% 4     {'LocusTag'                             }
% 5     {'Synonyms'                             }
% 6     {'dbXrefs'                              }
% 7     {'chromosome'                           }
% 8     {'map_location'                         }
% 9     {'description'                          }
% 10    {'type_of_gene'                         }
% 11    {'Symbol_from_nomenclature_authority'   }
% 12    {'Full_name_from_nomenclature_authority'}
% 13    {'Nomenclature_status'                  }
% 14    {'Other_designations'                   }
% 15    {'Modification_date'                    }
% 16    {'Feature_type'                         }

GeneColumnsToUse = [2 3 11 9 5 7 8 10 6];

tline = fgetl(FidInputFile);
tline =  textscan(tline,'%s','delimiter','\t');
ColumnIds = tline{1}
numColumns = numel(ColumnIds)

S = textscan(FidInputFile,repmat('%s',1,numColumns),'delimiter','\t');

GeneId = S{2};
Symbol = S{3};
Synonyms = S{5};
FullAnnotation = cat(2,S{1:numColumns});
Synonyms_Cell = cellfun(@(x) strsplit(x,'|'), Synonyms, 'UniformOutput', false);

nInput = numel(GeneSymbols);

GeneIdOut =  cell(nInput,1);
GeneSymbolOut  =  cell(nInput,1);
FullAnnotationOut = cell(nInput,numColumns);
GeneIdOut(:) = {'---'};
GeneSymbolOut(:) = {'---'};
FullAnnotationOut(:) = {'---'};

[~,indx1,indx2] = intersect(GeneSymbols,Symbol,'Stable');

GeneSymbolOut(indx1,:) = Symbol(indx2,:);
GeneSymbolOut(cellfun('isempty',GeneSymbolOut)) = {'---'};
    
GeneIdOut(indx1,:) = GeneId(indx2,:);
GeneIdOut(cellfun('isempty',GeneIdOut)) = {'---'};

FullAnnotationOut(indx1,:) = FullAnnotation(indx2,:);
FullAnnotationOut(cellfun('isempty',FullAnnotation)) = {'---'};


% get gene ids with no maching annotion
[MissingIds, indxMissing] = setdiff(GeneSymbols,Symbol,'stable');
length(MissingIds)

for i=1:length(MissingIds)
    indx_tmp = cellfun(@(x) strcmp(MissingIds{i},x),Synonyms_Cell,'UniformOutput',false);
    indx_tmp = cellfun(@(x) any(x),indx_tmp,'UniformOutput',false);
    indx_tmp = cell2mat(indx_tmp);
    
    if sum(indx_tmp) ==  1
            GeneSymbolOut(indxMissing(i)) = Symbol(indx_tmp);
            GeneIdOut(indxMissing(i)) = GeneId(indx_tmp);
            FullAnnotationOut(indxMissing(i),:) = FullAnnotationOut(indx_tmp,:);
            
    elseif sum(indx_tmp) > 1
        fprintf('Multiple synonyms found for %s\n',MissingIds{i});
        indx_tmp = find(indx_tmp);
        GeneSymbolOut(indxMissing(i)) = Symbol(indx_tmp(1));
        GeneIdOut(indxMissing(i)) = GeneId(indx_tmp(1));
        FullAnnotationOut(indxMissing(i),:) = FullAnnotationOut(indx_tmp(1),:);
        
    else
        fprintf('No synonyms found for %s\n',MissingIds{i}) 
    end
end
%NewGeneId =  GetReplacedGeneId(MissingIds);











