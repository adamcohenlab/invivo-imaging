function [vUnits] = ic_makeUnits(q)
%IC_MAKEUNITS create 'vUnit' objects (cells) from ica components defined by
%q.icells
%   Detailed explanation goes here
% todo:
%   -how to deal with running makeUnits twice... overrwrite?

if isempty(q.icells)
    %or instead plot the ICA components and ask to select cells by hand
    error('Please define indices of cell-like ICA components via obj.icells');
end

ncells = length(q.icells);

for i = 1:ncells
    iu = q.icells(i);
    v = VUNIT();
    v.parent = q;
    v.calcmethod = 'ICA';
    v.index = iu;
    v.footprint = q.footprints.frame(iu);
    v.filter = q.filters.frame(iu);
    v.timetrace = q.traces(iu,:);
    v.timetracelp = q.more.traceslp(iu,:);
    v.timetracenp = q.more.tracesnp(iu,:);
    v.com = []; %calculate COM position of footprint
    v.sigma = []; %rough 'half-width' estimate
    
    
    %calculate stats
    v.stats = [];
    vUnits(i)=v; %add cell to array of vUnits.
end

if isempty(q.vunits) %add vUnits array to ICAobj
    q.vunits = vUnits;
else
    q.vunits((end+1):(end+ncells))=vUnits;
end

end
    
    




