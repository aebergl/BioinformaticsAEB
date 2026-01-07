function [DirNames] = GetDirs(Dir)

    files = dir(Dir);
    isDir = [files.isdir];
    subDirs = files(isDir); % A structure with extra info.
    DirNames = {subDirs(3:end).name};
end