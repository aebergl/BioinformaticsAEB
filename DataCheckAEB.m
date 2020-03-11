function DATA = DataCheckAEB(DATA,varargin)

NaN_Check = true;
Zero_Check = true;
Std_Check = true;

[N,K] = size(DATA.X);

X = ones(N,K,'uint8');

X(isnan(DATA.X)) = 1;
X(DATA.X == 0) = 2;

% Create Ficure
fh=figure('Name','Data Check plot','Color','w','Tag','Data Check plot');
ah = axes(fh,'NextPlot','add','tag','MV and 0 figure','box','on','Linewidth',1,...
    'YAxisLocation','origin','Layer','top','position',[0.3 0.3 0.4 0.4]);

pcolor(ah,X)