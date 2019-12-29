function [h] = vmd_basicPlot(v,varargin)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
%   Simon kheifets 6/3/2008
% add...
%   PSD and autocorrelation of mean trace...
%   statistics abt image...
h.parent = v;

p = inputParser;
p.KeepUnmatched = 1;
addParameter(p,'topdf',0);
addParameter(p,'fixylims',1); %adjust the y limits for mean nad variancce
parse(p,varargin{:});
v2struct(p.Results);
if topdf
    h.fh = figure('Units','inches',...
        'Position',[17 1.5 8.5 11],...
        'Color','w',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperSize',[8.5 11]);
else
    h.fh = figure('Position',[403   100   837   651],'Color','w');
end

nplots = 4;
nrows = 5;
ncols = 1;
titlesize = 0.15;
hh = 0.9/nplots;
hhi = 0.7*hh;
ww = 0.8;
ll =(1-ww)/2;
postitle = [0,0.9,1,0.1];
pos1 = [ll,0.95-hh,ww,hhi];
pos2 = [ll,0.95-2*hh,ww,hhi];
pos3 = [ll,0.95-3*hh,ww,hhi];
pos4 = [ll,0.95-4*hh,ww,hhi];




axcount = 1;
titletext = ['Raw voltage movie summary:' newline v.label];
h.ah(axcount) = subplot('Position',postitle);
%set(gca,'box','on');
h.th = text(0.5, 0.5, titletext, 'Units', 'Normalized',...
    'HorizontalAlignment','Center',...
    'FontUnits','Normalized','FontSize',0.3,...
    'Interpreter','none');
if topdf
    set(h.th,'FontUnits','inches','FontSize',titlesize);
end
set(gca,'Visible','off')


axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos1);
imagesc(v.meanimg);
daspect([1 1 1]);
colorbar;
colormap(jet);
title('Mean of each pixel');

axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos2);
imagesc(log10(v.varimg));
daspect([1 1 1]);
colorbar;
colormap(jet);
title('Var of each pixel (log scale)');

axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos3);
plot(v.tvec,v.meantrace,'LineWidth',0.05);
title('Mean of each frame');
xlabel('Time (s)');
ylabel('mean counts');
if fixylims
    ylim([prctile(v.meantrace,.1) max(v.meantrace)]);
end
    

axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos4);
plot(v.tvec,v.vartrace,'LineWidth',0.05);
title('Variance of each frame');
xlabel('time (s)');
ylabel('var (counts^2)');
if fixylims
    ylim([prctile(v.vartrace,.1) max(v.vartrace)]);
end



end

