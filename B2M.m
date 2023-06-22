function M = B2M(B)

MaxVal = 18;

M = log2(B ./ (1-B));

M(M < -MaxVal) = -MaxVal;
M(M >  MaxVal) =  MaxVal;