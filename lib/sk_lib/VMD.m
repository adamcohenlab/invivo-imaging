classdef VMD < handle
    %VMD defines a voltage movie data class as
    %simply as possible, just so it inherits from handle, and has well
    %defined properties which reflect the type of data/analysis that it
    %stores.
    
    properties
        set = []; %import settings
        timestamp = []; %time of creation
        mov = [];
        label = [];
        tvec = [];
        dt = [];
        duration = [];
        meanimg = [];
        varimg = [];
        meantrace = [];
        vartrace = [];
        sgrm = []; %spectrogram information
%        info = [];
%        clicky = [];        %struct containing size and results from clicky
        fvm = [];
    end
    
    methods
        function s = saveobj(obj) %function to save without running out of memory
            fnames = fieldnames(obj);
            for i = 1:length(fnames)
                fname = fnames{i};
                if strcmp(fname,'mov')%
                    s.(fname) = [];
                else
                    s.(fname) = obj.(fname);
                end
            end
        end
    end
    methods (Static)
        function obj = loadobj(s) %load object, but load the movie based on its original location (so as not to save multiple times)
            newObj = VMD();
            fnames = fieldnames(newObj);
            for i = 1:length(fnames)
                fname =fnames{i};
                newObj.(fname)=s.(fname);
            end    
            [mov0, nrow, ncol]=readBinMov4(newObj.set.dir,newObj.set.file,newObj.set.transpose);
            newObj.mov=vm(mov0);
            
            if ~isempty(newObj.fvm)
                newObj.fvm.parent = newObj;
            end
            obj = newObj;
        end
    end
end
