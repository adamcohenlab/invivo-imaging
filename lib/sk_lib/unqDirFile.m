function [thisfilename] = unqDirFile(savedir,str,ext)

if ~exist(savedir,'dir')==7
    mkdir(savedir);
end



%generate unique filename
stopflag = 0;
i=1;
while ~stopflag
    thisfilename = fullfile(savedir,sprintf([str '_%i.' ext],i));
    if isa(thisfilename,'file')==2
        i=i+1;
    else
        stopflag=1;
    end
end

end

