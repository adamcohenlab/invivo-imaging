function [out] = parsavetiff(data,filename)
%PARSAVETIFF wrapper for saveastiff for use in parfor loop
%   Detailed explanation goes here

out = saveastiff(data, filename);

end

