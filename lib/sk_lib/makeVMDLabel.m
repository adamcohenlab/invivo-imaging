function [lb] = makeVMDLabel(pt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ivqn = str2num(cell2mat(regexp(pt,'IVQ(\d+)','tokens','once')));
daten = cell2mat(regexp(pt,'(\d+-\d+-\d+)','tokens','once'));
session = str2num(cell2mat(regexp(pt,'IVQ\d+-S(\d+)','tokens','once')));
slice = str2num(cell2mat(regexp(pt,'slice(\d+)','tokens','once')));
fov = str2num(cell2mat(regexp(pt,'FOV(\d+)','tokens','once')));
dirname = cell2mat(regexp(pt,'[\\\/](\d{6}.*)[\\\/]','tokens','once'));
[filepath,name,ext] = fileparts(pt);
timen = dirname(1:6);
datevector = (cellfun(@str2num,({daten(1:4); daten(6:7); daten(9:10); timen(1:2); timen(3:4); timen(5:6)})))';
t = datestr(datevector,'mmm dd yyyy,HH:MMAM');
lb = sprintf('IVQ%iS%i, slice%i, FOV%i | %s/%s | %s',...
    ivqn,session,slice,fov,dirname(1:end),[name ext],t);

end

