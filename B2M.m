function M = B2M(B)

MaxVal = 8;

M = log2(B ./ (1-B));

M(M < -MaxVal) = -MaxVal;
M(M >  MaxVal) =  MaxVal;