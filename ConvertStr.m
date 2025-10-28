function StrOut = ConvertStr(StrIn,OutPutType)


switch class(StrIn)
    case 'string'
        switch lower(OutPutType)
            case 'cell'
                StrOut = cellstr(StrIn);
            case 'char'
                StrOut = char(StrIn);
            case 'string'
                StrOut = StrIn;
        end

    case 'cell'
        switch lower(OutPutType)
            case 'cell'
                StrOut = StrIn;
            case 'char'
                StrOut = char(StrIn);
            case 'string'
                StrOut = string(StrIn);
        end

    case 'char'
        switch lower(OutPutType)
            case 'cell'
                StrOut = cellstr(StrIn);
            case 'char'
                StrOut = StrIn;
            case 'string'
                StrOut = string(StrIn);
        end

    otherwise
        error('Wrong input type')
end

