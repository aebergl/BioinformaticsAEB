function DATA = AddREFseqGeneInfoU133Plus2(DATA,U133_BA_Path,RefSeqPath,ReplaceAppend)


U133AnnotationFile = fullfile(U133_BA_Path,'HGU133Plus2_Hs_ENTREZG_desc.annot');
if isfile(U133AnnotationFile)
    DATA = AddColAnnotationFromFile(DATA,U133AnnotationFile);
else
    error('Could not load %s',U133AnnotationFile)
end

RefSeqFile = fullfile(RefSeqPath,'Homo_sapiens.gene_info');


if isempty(RefSeqFile)
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
    [FidInputFile,message] = fopen(RefSeqFile,'r');
    DATA.Info.GeneAnnotation = 'file_name';
    DATA.Info.GeneAnnotationData = date;
    
    if  FidInputFile == -1
        disp(RefSeqFile)
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
GeneColumnsToUse = [2 3 10 11 9 5 7 8 6];

tline = fgetl(FidInputFile);
tline =  textscan(tline,'%s','delimiter','\t');
ColumnIds = tline{1};
numColumns = numel(ColumnIds);
Annotation = cell(DATA.nCol,length(GeneColumnsToUse));
Annotation(:) = {'---'};

S = textscan(FidInputFile,repmat('%s',1,numColumns),'delimiter','\t');

GeneId = S{2};

GeneInfo = cat(2,S{GeneColumnsToUse});

[~,indx1,indx2] = intersect(DATA.ColAnnotation(:,2),GeneId,'Stable');

Annotation(indx1,:) = GeneInfo(indx2,:);
Annotation(cellfun('isempty',Annotation)) = {''};
Annotation(indx1,:) = GeneInfo(indx2,:);

% get gene ids with no maching annotion
[MissingIds, indxMissing] = setdiff(DATA.ColAnnotation(:,2),GeneId,'stable');
length(MissingIds)

% for i=1:length(MissingIds)
%     
%     if strcmpi(NewGeneId{i},'withdrawn') % Should not be used
%         Annotation(indxMissing(i),:) = {'Withdrawn'};
%         ExtraInfo{indxMissing(i)} = sprintf('%s is withdrawn',MissingIds{i});
%     elseif isempty(NewGeneId{i})
%         Annotation(indxMissing(i),:) = {'Missing'}; % No annotation was found
%         ExtraInfo{indxMissing(i)} = sprintf('%s could not be find',MissingIds{i});
%     else        
%         indx = strcmp(NewGeneId{i},DATA.ColId); % Check of this gene is already in the dataset
%         if any(indx)
%             Annotation(indxMissing(i),:) = {sprintf('%s has been replaced',MissingIds{i}, NewGeneId{i})};
%             ExtraInfo{indxMissing(i)} = sprintf('%s is replaced with %s and it already exist',MissingIds{i}, NewGeneId{i});
%         else
%             indx_tmp = strcmp(NewGeneId{i},GeneId);
%             Annotation(indxMissing(i),:) = GeneInfo(indx_tmp,:);
%             ExtraInfo{indxMissing(i)} = sprintf('%s replaced with %s',MissingIds{i}, NewGeneId{i});
%         end
%     end
%     
%     
% end

switch lower(ReplaceAppend)
    case 'replace'
        DATA.ColAnnotationFields = [ColumnIds(GeneColumnsToUse)];
        DATA.ColAnnotation = [Annotation];
    case 'append'
        DATA.ColAnnotationFields = [DATA.ColAnnotationFields; ColumnIds(GeneColumnsToUse)];
        DATA.ColAnnotation = [DATA.ColAnnotation Annotation];
end

end















