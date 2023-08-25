function [px,py] = GetPatch(nRows,nCols)
%Calculate Patch coordinates
%https://stackoverflow.com/questions/6614207/how-to-export-non-blurry-eps-images
px = bsxfun(@plus, [-0.5; 0.5; 0.5; -0.5], reshape(1:nCols, [1 1 nCols]));
py = bsxfun(@plus, [-0.5; -0.5; 0.5; 0.5], 1:nRows);

px = reshape(repmat(px, [1 nRows 1]), 4, nCols*nRows);
py = reshape(repmat(py, [1 1 nCols]), 4, nCols*nRows);

end