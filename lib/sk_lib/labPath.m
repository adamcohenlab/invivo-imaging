function [labpath] = labPath()
%Generate lab root directory path for the native environment
%   Detailed explanation goes here
labpathwin = 'X:\Lab\';
labpathunix = '/n/cohen_lab/Lab';
if ispc
    labpath = labpathwin;
else
    labpath = labpathunix;
end
end

