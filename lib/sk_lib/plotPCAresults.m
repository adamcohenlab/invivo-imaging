function [handles] = plotPCAresults(obj,traces,D,totvar,title)

ah = subplot('Position',[0.0 0 0.6 1]);
set(ah,'Visible','off');
handles.imhandles = stackImshow(obj,'maintitle',...
    title,'subtitles',sprintfc('N=%i',1:obj.frames),...
    'hframe',ah,'addcolorbar',1);
%
handles.tracehandle = subplot('position',[0.625 0.45 0.35 0.5]);
stackplot(traces);

handles.varhandle = subplot('position', [0.625 0.05 0.35 0.3]);
Dc = [0; cumsum(D)];
semilogy(0:1:length(D),(totvar-Dc)/totvar,'o-');
%ylim([1e-5 1]);
xlabel('Component #');
ylabel('fraction of variariance');
