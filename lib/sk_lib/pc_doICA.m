function [q] = pc_doICA(p,s)
%FVM_DOICA do pca on a filtered movie
%   make settings s using makeICAsettings.m
q=ICAobj;
q.timestamp = datestr(now);
%% do ICA on HP filtered movie
if isempty(s.pkeep)
    p_keep = 1:p.npc;
else
    p_keep = s.pkeep;
end

 %Aplha is the ratio between spatial and Temporal ICA (0 = only space, 0.99 = only time):

[icsFilter,icsFootprint,icsTime,nICs,mixmat,sepmat,extra] = vmICA(p.uvm(:,:,p_keep),p.d(p_keep),p.v(:,p_keep),s.nic,s.alpha,p.parent.tvec);
VT = (p.more.vnp)';
tracesnp = sepmat*VT(p_keep,:);

VT = (p.more.vlp)';
traceslp = sepmat*VT(p_keep,:);

%%
q.parent = p;
q.set = s;

q.filters = icsFilter;
q.footprints = icsFootprint;
q.traces = icsTime;
q.more.traceslp = traceslp;
q.more.tracesnp = tracesnp;
q.stats = extra;
q.nics = nICs;
q.mixmat = mixmat;
q.sepmat = sepmat;


if isempty(p.ICAobj)
    p.ICAobj = q;
else
    p.ICAobj(end+1) = q;
end

end