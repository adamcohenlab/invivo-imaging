classdef FVM < handle
    %FVM defines the filtered voltage movie class as
    %simply as possible, just so it inherits from handle, and has well
    %defined properties which reflect the type of data/analysis that it
    %stores.
    
    properties
        parent = [];
        set = [];
        timestamp = [];
        mov = [];
        tvec = [];
        spectrogram = [];
        meanimg = [];
        meantrace = [];
        duration = [];
        label = [];
        varimg = [];
        dmov = [];
        movHP = [];
        movLP = [];
        meanimgHP = [];
        meanimgLP = [];
        varimgHP = [];
        varimgLP = [];
        meantraceHP = [];
        meantraceLP = [];
        PCAobj = [];
        ICAobj = [];
    end
    
    methods
        function s = saveobj(obj) %function to save without running out of memory
            fnames = fieldnames(obj);
            for i = 1:length(fnames)
                fname = fnames{i};
                if strcmp(fname,'parent')% || strcmp(fname,'fvm')
                    s.(fname) = [];
                else
                    s.(fname) = obj.(fname);
                end
            end
        end
    end
    
        methods (Static)
        function obj = loadobj(s) %load fvm and restore references correctly
            newObj = FVM();
            fnames = fieldnames(newObj);
            for i = 1:length(fnames)
                fname =fnames{i};
                newObj.(fname)=s.(fname);
            end
            
            if ~isempty(newObj.PCAobj)
                newObj.PCAobj.parent = newObj;
            end
            
            if ~isempty(newObj.ICAobj)
                newObj.ICAobj.parent = newObj;
            end
            
                        
            obj = newObj;
        end
    end
    
    
end

        