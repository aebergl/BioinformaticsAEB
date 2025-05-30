function SavePDF_AEB(FileName)
fh=gcf;
set(fh, 'Units', 'inches');
set(fh, 'PaperUnits', 'inches');
set(fh, 'PaperSize', [fh.Position(3), fh.Position(4)]);
set(fh, 'PaperPosition', [0, 0, fh.Position(3), fh.Position(4)]);
set(fh, 'PaperPositionMode', 'auto');
print(fh, FileName, '-dpdf', '-vector');

end

