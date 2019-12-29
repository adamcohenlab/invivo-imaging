function [str_out] = vmd_toString(obj)
%vmd_toString generates a string array that contains all the settings and
%properties of the object.

if isa(obj,'VMD')
    s1 = strsplit(string(evalc('obj')),newline);
    s2 = strsplit(string(evalc('obj.set')),newline);

    s3 = ([s2(1,4:end) s1(1,5:end)]);
elseif isa(obj,'FVM')
    s1 = strsplit(string(evalc('obj')),newline);
    s2 = strsplit(string(evalc('obj.set')),newline);
    s3 = ([s2(1,4:end) s1(1,6:end)]);
elseif isa(obj,'PCAobj')
    s1 = strsplit(string(evalc('obj')),newline);
    s2 = strsplit(string(evalc('obj.set')),newline);
    s3 = ([s2(1,4:end) s1(1,4:end)]);
else
    s1 = strsplit(string(evalc('obj')),newline);
    s2 = strsplit(string(evalc('obj.set')),newline);
    s3 = ([s2(1,4:end) s1(1,4:end)]);
end

str_out = s3;
% th = text(1,1,s3);
% set(th,'FontName','FixedWidth');