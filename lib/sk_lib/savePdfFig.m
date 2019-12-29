function thisfilename = savePdfFig(h,vmd,str)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
savedir = fullfile(vmd.set.dir, 'quivvia');

thisfilename = unqDirFile(savedir,str,'pdf');
figure(h.fh);
print('-dpdf','-painters',thisfilename);
thisfilename2 = unqDirFile(savedir,str,'fig');
savefig(thisfilename2);

%save to appropriate folder
%do .ps first (to make multiple pages) and then convert to pdf
%determine path from vmd.set.dir
%or save as both ps and pdf, so that ps can be collated for summaries, but
%pdf is available locally within the file structure.
end

