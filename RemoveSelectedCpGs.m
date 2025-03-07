function [DATA] = RemoveSelectedCpGs(DATA,varargin)

KeepCgOnly      = true;
RemoveMasked    = false;
RemoveX         = false;
RemoveY         = false;

ChrColumnId     = "CHR";
ChrY            = "Y";
ChrX            = "X";