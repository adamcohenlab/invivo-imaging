function [h] = ic_quickPlot(q,varargin)
%IC_QUICKPLOT plot summary of ICAobj object, aka results of ICA done on the
%output of PCA analysis
%Simon Kheifets 6/20/2018
%   Detailed explanation goes here
% todo:
% -fix title
h.parent = 1;

pp = inputParser;
pp.KeepUnmatched = 1;
addParameter(pp,'topdf',0);
parse(pp,varargin{:});
v2struct(pp.Results);
titlesize=0.15; %font size for title

if topdf
    h.fh = figure('Units','inches',...
        'Position',[17 1.5 8.5 11],...
        'Color','w',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperSize',[8.5 11]);
else
    h.fh = figure('Position',[324   145   973   643],'Color','w');
end

if topdf
    htit = 0.1;
    htop = 0.55;
    wimg = 0.7;
    hfrac = 0.25;
    %position of title
    postitle = [0,1-htit,1,htit];
    
    %position of images
    hh = htop;
    ww = wimg;
    hhi = htop;
    wwi = wimg;
    ll = 0.5*(ww-wwi);
    bb = 1-htit-hh+0.5*(hh-hhi);
    pos1 = [ll bb wwi hhi];
    
    %position of timetraces
    hh = 1-htop-htit;
    ww = 1;
    hhi = 0.8*hh;
    wwi = 0.8*ww;
    bb = 0.5*(hh-hhi);
    ll = 0.5*(ww-wwi);
    pos2 = [ll bb wwi hhi];
    
    %position of fractional variance
    hh = hfrac;
    ww = 1-wimg;
    hhi = 0.8*hh;
    wwi = 0.8*ww;
    ll = wimg+0.5*(ww-wwi);
    bb = 1-htit-htop+0.5*(hh-hhi);
    pos3 = [ll bb wwi hhi];
    
    %position of info string
    hh = htop-hfrac;
    ww = 1-wimg;
    h0 = 0.8*hh;
    w0 = 0.8*ww;
    ll = wimg+0.5*(ww-w0);
    bb = 1-htit-hh+0.5*(hh-h0);
    pos4 = [ll bb w0 h0];
else
    pos1 = [0.0 0 0.6 1];
    pos2 = [0.625 0.45 0.35 0.45];
    pos3 = [0.625 0.05 0.35 0.3];
end
%% plot title
axcount = 1;
titletext = ['ICA summary:' newline q.parent.parent.parent.label];
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
%% plot images
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos1);
set(gca,'Visible','off');
h.imh = stackImshow(q.footprints,'subtitles',sprintfc('N=%i',1:q.footprints.frames),...
    'hframe',h.ah(axcount),'addcolorbar',0);
if topdf
%     set(h.ah(axcount),'FontUnits','inches',...
%         'FontSize',titlesize,'Interpreter','none');
end

%% plot time traces
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos2);
stackplot(q.traces',q.parent.parent.tvec);
xlabel('time');
title('component time traces');

ylims = get(gca,'YLim');
ylim([ylims(1) size(q.traces,1)+1]);
%%
%% plot fractional variance
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos3);

D = q.parent.d;
tv = q.parent.stats.totvar;
Dc = [0; cumsum(D)];
semilogy(0:1:length(D),(tv-Dc)/tv,'o-');

title('Fractional Variance of pca....');
xlabel('Component number');
ylabel('fraction of variance remaining');

%% plot info string
axcount = axcount+1;
h.ah(axcount)=subplot('Position',pos4);
infostring = vmd_toString(q);
h.th2 = text(-0.5,1,infostring,...
    'Interpreter','none',...
    'VerticalAlignment','top',...
    'FontName','FixedWidth');
set(gca,'Visible','off')



end

