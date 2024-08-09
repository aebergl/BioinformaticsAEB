function DATA = Methylation_ExpandExtendedImmune(DATA)

% Check that all cell types are present
CellTypes = {'Bas','Bmem','Bnv','CD4mem','CD4nv','CD8mem','CD8nv','Eos','Mono','Neu','NK','Treg'};

indx =  ismember(DATA.ColId,CellTypes);
if sum(indx) ~= 12
    error('Not all celltypes found in DATA')
end

New_ColId = {'Gran','Bcells','CD4T','CD8T','Tcells','Lymphocytes'}';
X_new = zeros(DATA.nRow,length(New_ColId));
 
% 1. Granulocytes
indx = ismember(CellTypes,{'Bas','Eos','Neu','Mono'});
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,1) = Xsum;

% 2. B cells
indx = ismember(CellTypes,{'Bmem','Bnv'});
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,2) = Xsum;

% 3. CD4 T cells
indx = ismember(CellTypes,{'CD4mem','CD4nv','Treg'});
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,3) = Xsum;

% 4. CD8 T cells
indx = ismember(CellTypes,{'CD8mem','CD8nv'});
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,4) = Xsum;

% 5. T cells
indx = ismember(CellTypes,{'CD4mem','CD4nv','Treg','CD8mem','CD8nv'});
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,5) = Xsum;

% 6. Lymphocytes
indx = ismember(CellTypes,{'Bmem','Bnv','CD4mem','CD4nv','Treg','CD8mem','CD8nv','NK'});
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,6) = Xsum;

DATA.X = [DATA.X, X_new];
DATA.ColId = [DATA.ColId; New_ColId];
DATA.nCol = size(DATA.X,2);