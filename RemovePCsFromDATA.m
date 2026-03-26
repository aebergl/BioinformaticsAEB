function DATA = RemovePCsFromDATA(DATA,PCsToRemove,varargin)

ScaleX = false;
MaxPCAComp = max(PCsToRemove);


[PCA_Model,X] = NIPALS_PCA(DATA.X,'NumComp',MaxPCAComp,'ScaleX',ScaleX,'MVAverage',false);

AllComps = 1:MaxPCAComp;
CompsToAddBack = setdiff(AllComps,PCsToRemove);

if ~isempty(CompsToAddBack)
    X = X + PCA_Model.T(:,CompsToAddBack)*PCA_Model.P(:,CompsToAddBack)';
end
if ScaleX
    PCA_Model.x_weight(PCA_Model.x_weight < 1e-8) = 1e-8;
    %X = bsxfun(@rdivide, X, PCAmodel.x_weight); %same as X = X .* (ones(N,1) * x_weight);
    X = X ./ PCA_Model.x_weight;
end

%X = bsxfun(@plus, X, PCA_Model.x_mean); %same as X = X - (ones(N,1) * x_mean);
X = X + PCA_Model.x_mean;
DATA.X = X;
