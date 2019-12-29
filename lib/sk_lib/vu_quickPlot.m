function [h] = vu_quickPlot(vus)
%VU_QUICKPLOT plot summary of array of vUnit objects, aka putative cells
%footprints and COM (for now)
%Simon Kheifets 6/20/2018
%   Detailed explanation goes here
% todo:
%   -plot time traces on one axes so they can be zoomed together...
%       (or write some callback function to update the zoom on the others if one 
%       of them is zoomed.

h.fh = figure('Position',[324   145   973   643],'Color','w');

nrows = length(vus);
ncols = 2;
postitle = [0,0.9,1,0.1];
w0(1) = 0.4;
w0(2) = 0.6;
h0 = 0.9/nrows;
ww = w0*0.8;
hh = h0*0.65;
lls = cumsum(w0);
lls = [0 lls];
for i = 1:nrows
    b = 0.9-i*h0+0.5*(h0-hh); %bottom
    for j=1:ncols
        l = lls(j)+0.5*(w0(j)-ww(j));
        pos{i,j} = [l b ww(j) hh];
    end
end

titletext = ['Putative Cells'];

axcount = 1;
h.ah(axcount) = subplot('Position',postitle);
h.th = text(0.5, 0.5, titletext, 'Units', 'Normalized',...
    'HorizontalAlignment','Center',...
    'FontUnits','Normalized','FontSize',0.3);
set(gca,'Visible','off')


for i = 1:nrows
    axcount = axcount+1;
    h.ah(axcount) = subplot('Position',pos{i,1});
    
    f = vus(i).footprint;
%     tt = vus(i).timetrace;
    tt = vus(i).timetracenp;
    tvec = vus(i).parent.parent.parent.tvec;
    
    imshow(f,[]);
    colorbar();
    set(gca,'CLim', max(abs(f(:)))*[-1 1]);
    colormap(bluewhitered);
    
    axcount = axcount+1;
    h.ah(axcount) = subplot('Position',pos{i,2});
    plot(tvec,tt);
end
    

    
    
    
end