function [p] = fvm_doPCA(f,s)
%FVM_DOPCA do pca on a filtered movie
%   make settings s using makePCAsettings.m
p=PCAobj;
p.timestamp = datestr(now);
%% do PCA on HP filtered movie
[p.uvm,p.d,p.v,p.stats] = vmPCA(f.movHP,s.npc,f.tvec,0); %do PCA

%% apply same filters to raw movie
[m.dnp,m.vnp,m.statsnp] = applyFilters(p.uvm,f.dmov,f.tvec); %apply filters to another version of the movie

%% apply same filters to LP movie
[m.dlp,m.vlp,m.statslp] = applyFilters(p.uvm,f.movLP,f.tvec); %apply filters to another version of the movie

p.parent = f;
p.set = s;
p.more = m;
p.npc = s.npc;
if isempty(f.PCAobj)
    f.PCAobj = p;
else
    f.PCAobj(end+1) = p;
end

end