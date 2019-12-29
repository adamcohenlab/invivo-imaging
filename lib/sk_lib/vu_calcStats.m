function vus = vu_calcStats(vus)
%VU_CALCSTATS take a vu object with raw data (footprint, timetrace) and
%calculate some useful parameters/statistics...
for i = 1:length(vus)
    v = vus(i);
    fp = v.footprint;
    gvx = 1:size(fp,2);
    gvy = 1:size(fp,1);
    [X,Y] = meshgrid(gvx,gvy);
    comx = X.*fp/sum(fp(:));
    comx = sum(comx(:));
    comy = Y.*fp/sum(fp(:));
    comy = sum(comy(:));
    
    
    v.com = [comy comx];
    
end