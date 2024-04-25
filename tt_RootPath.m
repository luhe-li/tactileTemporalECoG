function rootPath = tt_RootPath()

rootPath = which('tt_RootPath');
rootPath = fileparts(rootPath);
addpath(genpath(rootPath))

end