classdef VUNIT < handle
    %vUnit defines a putative cell, as identified from a vi movie via ICA, clicky,
    %or some other method. it has an image, COM, and a time trace

    
    properties
        parent = [];
        calcmethod = [];
        index = [];
        timestamp = [];
        footprint = [];
        filter = [];
        timetrace = [];
        timetracelp = [];
        timetracenp = [];
        com = []; %center of mass (pixels; relative to filtered movie...)
        sigma = []; %second moment
        stats = [];
    end
    
    methods
        function s = saveobj(obj) %function to save object to struct
            fnames = fieldnames(obj);
            for i = 1:length(fnames)
                fname = fnames{i};
                if strcmp(fname,'parent')% dont save parent, otherwise it gets into an infinite loop
                    s.(fname) = [];
                else
                    s.(fname) = obj.(fname);
                end
            end
        end
    end
    
            methods (Static)
        function obj = loadobj(s) %load vunit and restore references correctly
            newObj = ICAobj();
            fnames = fieldnames(newObj);
            
            for i = 1:length(fnames)
                fname =fnames{i};
                newObj.(fname)=s.(fname);
            end
 
            obj = newObj;
        end
    end
end