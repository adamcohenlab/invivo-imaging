function [v,f] = vmd_makeFVM(v,s)
%VMD_MAKEFVM makes a filtered voltage movie out of a raw movie.
% use makeFVMsettings() to generate settings struct 's';
%   Simon kheifets 6/3/2018

if ~isempty(s.tlim)
    mov = (v.mov(:,:,v.tvec>s.tlim(1)&v.tvec<s.tlim(2)));
    tvec = v.tvec(v.tvec>s.tlim(1)&v.tvec<s.tlim(2));
else
    mov = v.mov;
    tvec = v.tvec;
end
    
mov = mov-400; %correct hammamatsu offset

if s.cropmotion % crop motion edges
    %this does not account for the cropped time yet
    [mov, xShifts, yShifts] = fixMotionEdges(mov,v.set.dir);
end

%mov = v.mov(:,3:end,:); %get rid of the bad edge pixel

if s.dopblc
    mov = mov.pblc; %photobleach correction
end

if s.removerownoise
    mov = rmRowNoise(mov,tvec); %remove row-symmetric noise
end

%%
if~isempty(s.nbin)&&s.nbin>1
    mov = do_binning(mov,s.nbin); %do binning in x and y
end
%mov=mov./mov.mean;
%mov = mov.*sqrt(mov.mean); %normalize by stdev so shot noise is similar magnitude in each pixel)

%% fix flashing DMD frames
if s.fixdmdframes
    intensS = smooth(mov.frameAverage,100);% plot(tvec,[mov.frameAverage intensS]);
    intensHP = mov.frameAverage-intensS;
    badFrames = find(abs(intensHP)> 4);
    mov(:,:,badFrames) = mov(:,:,badFrames-1);
end


%% do low/high pass
dmov = mov - mean(mov);
% movLP = imfilter(dmov, ones(1,1,20)/20, 'replicate');
movLP = vm((filter(ones(1,s.hpsamples)/s.hpsamples,1,(dmov(:,:).data)'))',[mov.rows,mov.cols]);
movHP = dmov - movLP;

indkeep = s.hpsamples:movHP.frames;
movHP = movHP(:,:,indkeep);
movLP = movLP(:,:,indkeep);
mov = mov(:,:,indkeep);
dmov = dmov(:,:,indkeep);
tvec = tvec(indkeep);

%% remove bloodflow-type pixles
if s.rmbloodflow %try to remove bloodflow-type pixels
    difmov = vm(diff(movLP.data,1,3));
    nfilt = 500;
    difmovLP = vm((filter(ones(1,nfilt)/nfilt,1,(difmov(:,:).data)'))',[difmov.rows,difmov.cols]);
    difmovHP=difmov-difmovLP;

    difmovHP = difmovHP(:,:,nfilt:end);

    vardif = var(difmovHP.data,0,3);
    %imagesc(vardif);daspect([1 1 1]);colorbar
    indreject = find(vardif(:)>s.tblood);
    movHP(indreject,:) = 0;
    movLP(indreject,:) = 0;
    mov(indreject,:) = 0;
    dmov(indreject,:) = 0;
end

%% save data in struct
f = FVM();
f.timestamp = datestr(now);
f.parent = v;
f.set=s;
f.mov=mov;
f.tvec = tvec;

f.meanimg = f.mov.mean;
f.meantrace = f.mov.frameAverage;
f.duration = f.tvec(end);
if isempty(s.label)
    f.label = v.label;
else
    f.label = s.label;
end
    
f.varimg = var(single(f.mov.data),0,3);
f.dmov = dmov;
f.movHP = movHP;
f.movLP = movLP;
f.meanimgHP = f.movHP.mean;
f.meanimgLP = f.movLP.mean;
f.varimgHP = var(movHP);
f.varimgLP = var(movLP);
f.meantraceHP = f.movHP.frameAverage;
f.meantraceLP = f.movLP.frameAverage;

if isempty(v.fvm)
    v.fvm = f;
else
    v.fvm(end+1)=f;
end

end

