function [handles] = plotComponents(obj,traces,title)
ah = subplot('Position',[0.0 0 0.6 1]);
set(ah,'Visible','off');
handles.imhandles = stackImshow(obj,'maintitle',title,'subtitles',sprintfc('N=%i',1:obj.frames),'hframe',ah);
%
handles.tracehandle = subplot('position',[0.625 0.1 0.35 0.8]);
stackplot(traces);