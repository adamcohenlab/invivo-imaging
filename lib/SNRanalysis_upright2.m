%% default ROI = 1;
function [SNR, dFF, Tconv, C,thresh] = SNRanalysis_upright2(clicktrace, Frequency,threshscale, ROI, ChR);

if nargin< 2;
    ROI = 1;
    ChR = 0;
end

if nargin>2 & nargin<3;
    ChR =0;
end

 min_width = 20;
for ri = ROI; % pick the ROI to analyze
    figure(881);clf;
    trace1 = clicktrace(:,ri); % avoid the first several frames of photobleach artifact
    time = (1:length(trace1))/Frequency;
    plot(time,trace1);
    %     trace1 = (trace0-min(trace0))/(max(trace0)-min(trace0));
    %% % get rid of blue artifacts and fluctuation
    t = time;
    %background
    bkgd_interp_width = 10;     % number of points to smooth over when calculating background
    smpix = 9;                  % the number of pixels to use in the smoothing function
    sigma = 0.8;                % the width of the gaussian in the smoothing function in frames
    min_interspike_time = 3;    % the minimum number of frames between spikes
    dt = 1000/Frequency; % ms
    %smooth function % smooth with peak
    x = -ceil((smpix-1)/2):ceil((smpix-1)/2);
    f_sm = exp(-x.^2/(2*sigma^2));
    f_sm = f_sm/sum(f_sm);
    %     figure(1)
    %     plot(x,f_sm,'-*'); title('smoothing function');
    lensm = (length(f_sm)-1)/2;
    %
    FramesPerStep=length(trace1);
    T = trace1;
%     [T, pbleach] = rem_pbleach(trace1,ceil(FramesPerStep*3/5));     % remove long-timescale photobleaching.  Interp value shoule be longer than blue pulse duration
    bkgd = medfilt2(T, [bkgd_interp_width 1],'symmetric');  % use median filter to calculate background
%     bkgd = smooth(T, 20);
    Tmed = T - bkgd;
%     a = 1/5;
%     Tmed = filter([1-a a-1], [1 a-1], T);
%     bkgd = T - Tmed;
    
    Tconv = conv([ones(1,lensm)*Tmed(1), Tmed', ones(1,lensm)*Tmed(end)],f_sm,'valid');  % smooth using calulated smooth function
    Sm = round(FramesPerStep/8);
    %Sm = bkgd_interp_width;
    Tconv_bkgd = smooth(ordfilt2(Tconv, round(0.1*Sm), ones(1,Sm),'symmetric'),Sm); % raise up the regions that stick down too far during spiking
    
    figure(882); clf;
    plot(t,T,t,bkgd,t,Tconv,t,Tconv_bkgd)
    xlabel('time (sec)'); ylabel('fluorescence intensity');
    legend('raw','bkgd','conv','conv bkgd')
    hold on
%% try to define a new threshold
%     Tsort = sort(T);
%     figure; hold all; plot(1:length(Tsort), Tsort,'-o');
%     Tturn = Tsort(2:end)-Tsort(1:(end-1));
%     figure; plot(Tturn);
%     Tyaxismin = Tsort (ceil(length(Tsort)*0.005));
% Tyaxismax = Tsort (ceil(length (Tsort) *0.995));
% figure; subplot (3,1,1); plot (T);
% ylim ([Tyaxismin, Tyaxismax]);
% title (['Original Traces']);

    %%
    pulsei = 1
    Tstep = Tconv(FramesPerStep*(pulsei-1) + (1:FramesPerStep));
    if threshscale ==0;
        guide = mean(Tstep) + 1/2*(max(Tstep)-mean(Tstep));
        figure; plot(Tstep);
        title('Right-click to indicate threshold')
        [x, y] = getpts(gca);
        thresh = y(end);
    else
        thresh = threshscale;
    end
   
    [sT, ns] = spikefind2_sl(Tstep,min_interspike_time,thresh,thresh,min_width);
    
    % fit each spike to a Gaussian to extract the width and fine-tune
    % the timing
    sT2 = zeros(1,ns);
    sW = zeros(1,ns);
    Amp = zeros(1,ns);
    % Set up a kernel for the average spike waveform
    nback = 40/dt;
    nfront = 60/dt;
    Lk = nback + nfront + 1;
    kernel = [];
    %         ChR= 0;
    %%
    if ChR
        kernel1 = [];
        %             [lighton, lightoff] = ChR_timing(bkgd,interval,threshscale)
         % try to find the light on and off timing by directly extracting the trace (if there is a leak through)
%         [lighton, lightoff] = ChR_timing(bkgd,20,20);
%         %             lighton = 1:500:length(Tstep);  % hard coding the times the light turns on.  Need to confirm for your data
%         %             lightoff = 250:500:length(Tstep);  % hard coding the times the light turns off.  Need to confirm for your data
%         if isempty(lighton);
        display('did not find the edges of stimulation on and off by ChR_time');
        cyclelength = 1 %second
        stimon = 0.5 %second
        stimoff = cyclelength-stimon; %second
        cycn = floor(length(clicktrace)/cyclelength/Frequency);
        for icycn = 1:cycn;
            lighton(icycn) = round((cyclelength*(icycn-1))*Frequency)+1;
            lightoff(icycn) = round((cyclelength*(icycn-1)+ stimon)*Frequency)+1;
        end
        figure(995);clf; hold all; plot(clicktrace); plot(lighton, clicktrace(lighton),'r*');
        plot(lightoff, clicktrace(lightoff),'g*');
%         end
        firstspikedelay = 10;  % maximum delay (in frames) for the first spike after illumination to count as a "first spike"
        c1 = 1;
        
        
    end;
  
    %%
      c = 1; % counter of good spikes;%%
    if ~isempty(ns) && ns~=0;
        figure(886); clf; plot(Tconv);
        for k = 1:ns
            Frame = max([1,sT(k)-lensm]):min([length(Tstep),sT(k)+lensm]);    % pull out indices a few points on either side of spike
            Peak = Tstep(Frame);    % corresponding intensity values
            Time = Frame - sT(k);
            ffun = fittype('1 + A*exp(-4*log(2)*(x-t0)^2/w^2)'); % Fit function.  A is amplitude, t0 is time shift, w is full width half max
            fopt = fitoptions('Method','NonlinearLeastSquares','Lower',[0, -lensm, 1],'Upper',[1.5*(max(Peak)-1), lensm, 2*lensm+1],'StartPoint',[max(Peak)-1, 0, lensm/2]);   % specify min, max, and guess for each parameter
            peak_fit = fit(Time',Peak',ffun,fopt);
            Frame_interp = interp1(1:length(Frame),Frame,1:.1:length(Frame));    % plot fit using more finely spaced points
            figure(886); hold on;
            plot(Frame_interp,1 + peak_fit.A*exp(-4*log(2)*(Frame_interp - sT(k) - peak_fit.t0).^2/peak_fit.w^2),'c')
            text(sT(k) + peak_fit.t0,1+peak_fit.A,int2str(k));
            sT2(k) = sT(k) + peak_fit.t0;
            %Amp(k) = peak_fit.A;
            Amp(k) = max(Peak)-min(Peak);
            %                 figure(999); subplot(ceil(ns^0.5),ceil(ns^0.5),k); plot(Peak); hold all;
            % get half amplitude AP width
            for iP = 1: (length(Peak)-1);
                if Peak(iP+1)>= min(Peak)+ Amp(k)*0.5 && Peak(iP)<= min(Peak)+Amp(k)*0.5
                    Time1 = iP;
                end
                if Peak(iP+1)<= min(Peak)+Amp(k)*0.5 && Peak(iP)>= min(Peak)+Amp(k)*0.5
                    Time2= iP; %in the unit of frames
                end
            end
            
            try
                sW(k) = (Time2-Time1)*dt; % width of the action potential in the unit of time ms
            catch
                display('error finding half amplitude width time');
            end
            
            Traces = trace1';
            %             j = ROI;
            if ~ChR
                if (sT(k) > nback) & ((sT(k) + nfront) <= length(Tstep));
                    Tsingle = Traces((sT(k)-nback):(sT(k)+nfront));
                    Tnorm = (Tsingle-min(Tsingle))/(max(Tsingle)-min(Tsingle));
                    kernel = [kernel; Tnorm];
                    c = c + 1;
                end;
                %                 kernel = kernel/c;
            else  % if ChR
                if (sT(k) > nback) & ((sT(k) + nfront) <= length(Tstep));  % Check that the spike isn't too close to front or back of data
                    if find((sT(k) - lighton > 0) & (sT(k) - lighton < firstspikedelay))  % Check whether a "first spike"
                        Tsingle = Traces((sT(k)-nback):(sT(k)+nfront));
                        Tnorm = (Tsingle-min(Tsingle))/(max(Tsingle)-min(Tsingle));
                        kernel1 = [kernel1; Tnorm];
                        c1 = c1 + 1;
                        hold all;
                        plot(sT(k), Tconv(sT(k)), 'ro');
                    elseif ~sum((sT(k) < lightoff) & (sT(k) + nfront > lightoff))
                        Tsingle = Traces((sT(k)-nback):(sT(k)+nfront));
                        Tnorm = (Tsingle-min(Tsingle))/(max(Tsingle)-min(Tsingle));
                        kernel = [kernel; Tnorm];
                        c = c + 1;
                    end;
                end;
            end;
            %calculate noise with signal
            Signal(k) = T(sT(k))-bkgd(sT(k));
            %              spikedelete = T;
            %              spikedelete((sT(k)-nback):(sT(k)+nfront)) = 1;
        end
        %         kernel = kernel/c;
        %         if ChR
        %         kernel1 = kernel1/c1;
        %         end
        title({['cell ',int2str(1),', pulse ',int2str(ri)];'Press space to input bad indices. To keep all peaks press any other key.'});
        hold off;
        %look at the kernels
        figure(883);clf;
        tau = (-nback:nfront)*dt;
        if ~ChR
            try
                plot(tau, kernel,'k','linewidth',1); hold on;
                plot(tau, mean(kernel),'r','linewidth',2); hold on;
                %             legend('individual kernel','average kernel');
                xlabel('Time (ms)')
                ylabel('Fluorescence (A.U.)')
                title(['Average spike waveform (n = ' num2str(ns) ' spikes)'])
            catch
                display('No spikes found');
            end
        else
            if ~isempty(kernel1);
                subplot(1,2,1);
                plot(tau, kernel1,'k','linewidth',1); hold on;
                plot(tau, mean(kernel1),'r','linewidth',2); hold on;
                xlabel('Time (ms)')
                ylabel('Fluorescence (A.U.)')
                title(['Average first spike waveform(n = ' num2str(c1) ' spikes)'])
            else
                display('No first spike in stepstim');
            end
            if ~isempty(kernel)
                subplot(1,2,2);
                plot(tau, kernel,'k','linewidth',1); hold on;
                plot(tau, mean(kernel),'r','linewidth',2); hold on;
                xlabel('Time (ms)')
                ylabel('Fluorescence (A.U.)')
                title(['Average not-first spike waveform(n = ' num2str(c1) ' spikes)'])
            end
        end;
      
        dFF = mean(Amp)*100; %percentage of dF/F
        C(ri).nspike(1) = ns;
        C(ri).spikeT{1} = sT2;
        C(ri).spikeW{1} = sW;
        C(ri).spikeA{1} = Amp;
        C(ri).Tconv{1} = Tconv;
        try
            C(1).kernel{1,2} = kernel';
            C(1).kernel{1,1} = tau';
        end
        C(ri).C=c;
        if ChR
            try
                C(ri).kernel1{1,2} = kernel1';
                C(ri).kernel1{1,1} = tau';
                C(ri).C1 = c1;
            end
        end;
        % calculate singal to noise
        sig = mean(Signal);
        Bkgd = T - bkgd;
        
        for sTi = 1:ns;
            try
                Bkgd((sT(sTi)-nback):(sT(sTi)+nfront)) = 0;
            end
        end
        figure(887);
        plot(Bkgd);
        title('Bkgd trace without real trace');
        noise =std(Bkgd);
        SNR = sig/noise;
        display(['SNR = ', num2str(SNR)]);
        
        
        figure(884); clf;
        % plot original trace
        subplot(3,8,1:8);
        plot(time,trace1,'b','linewidth',2);
        xlabel('Time (s)','Fontsize',14);
        ylabel('Raw fluorescence (A.U.)', 'Fontsize', 14);
%         xlim([0,10]);
        try
            title(['SNR = ', num2str(SNR), '  dF/F = ', num2str(dFF),'%'],'Fontsize',14);
        catch
            title(['No spikes found'],'Fontsize',14);
        end
        legend('Raw trace');
        % plot the convoluted trace
        subplot(3,8,9:16);
        plot(time,Tconv,'k','linewidth',2);
        xlabel('Time (s)','Fontsize',14);
        ylabel('dF/F', 'Fontsize', 14);
%         xlim([0,10]);
        title(['Convoluted trace ', num2str(ri)],'Fontsize',14);
        legend('Convoluted trace');
        % plot kernels
        if ~ChR
            try
                subplot(3,8,17:18);
                plot(tau, kernel,'k','linewidth',1); hold on;
                plot(tau, mean(kernel),'r','linewidth',3); hold on;
                %             legend('individual kernel','average kernel');
                xlabel('Time (ms)')
                ylabel('Fluorescence (A.U.)')
                title(['Average spike waveform (n = ' num2str(ns) ' spikes)'])
            end
        else
            subplot(3,8,17:18);
            try
                plot(tau, kernel1,'k','linewidth',1); hold on;
                plot(tau, mean(kernel1),'r','linewidth',3); hold on;
                %             legend('individual kernel','average kernel');
                xlabel('Time (ms)')
                ylabel('Fluorescence (A.U.)')
                title(['Stepstimulation, first spike (n = ' num2str(c1) ' spikes)'])
            end
            if ~isempty(kernel)
                subplot(3,8,20:21);
                plot(tau, kernel,'k','linewidth',1); hold on;
                plot(tau, mean(kernel),'r','linewidth',3); hold on;
                %             legend('individual kernel','average kernel');
                xlabel('Time (ms)')
                ylabel('Fluorescence (A.U.)')
                title(['Stepstimulation, non-first spike (n = ' num2str(ns) ' spikes)'])
            end
            %     legend(['dF/F = ', num2str(dFF),'%']);
        end
        
    else
        SNR = [];
        dFF = [];
        C = [];
        display('no spike was found');
    end
end