function DATA = Methylation_ExpandExtendedImmune(DATA)

% ExtendedBlood
CellTypesExtendedBlood = {'Bas','Bmem','Bnv','CD4mem','CD4nv','CD8mem','CD8nv','Eos','Mono','Neu','NK','Treg'};

%EpiDish
CellTypesEpiDish = {'Baso','Bmem','Bnv','CD4Tmem','CD4Tnv','CD8Tmem','CD8Tnv','Eos','Mono','Neu','NK','Treg'};



if sum(ismember(DATA.ColId,CellTypesExtendedBlood)) == 12
    CellTypes = CellTypesExtendedBlood;
elseif  sum(ismember(DATA.ColId,CellTypesEpiDish)) == 12
    CellTypes = CellTypesEpiDish;
else
    error('Not all celltypes found in DATA')
end


New_ColId = {'Gran','Bcells','CD4T','CD8T','Tcells','Lymphocytes'}';
X_new = zeros(DATA.nRow,length(New_ColId));

% 1. Granulocytes
indx = ismember(DATA.ColId,CellTypes([1 8 9 10]));
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,1) = Xsum;

% 2. B cells
indx = ismember(DATA.ColId,CellTypes([2 3]));
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,2) = Xsum;

% 3. CD4 T cells
indx = ismember(DATA.ColId,CellTypes([4 5 12]));
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,3) = Xsum;

% 4. CD8 T cells
indx = ismember(DATA.ColId,CellTypes([6 7]));
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,4) = Xsum;

% 5. T cells
indx = ismember(DATA.ColId,CellTypes([4 5 6 7 12]));
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,5) = Xsum;

% 6. Lymphocytes
indx = ismember(DATA.ColId,CellTypes([2 3 4 5 6 7 11 12]));
Xtmp = DATA.X(:,indx);
Xsum = sum(Xtmp,2,"omitmissing");
X_new(:,6) = Xsum;

DATA.X = [DATA.X, X_new];
DATA.ColId = [DATA.ColId; New_ColId];
DATA.nCol = size(DATA.X,2);