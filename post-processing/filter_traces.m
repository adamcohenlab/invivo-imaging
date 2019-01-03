% addpath '/Volumes/cohen_lab/Lab/Labmembers/Shan Lou/Lab updates/2017 Ai155 trace analysis/subfunction/';
% addpath '/Volumes/cohen_lab/Lab/Labmembers/Shan Lou/code/ultrawidefield_analysis/SNR function package/';

%% get spikes
clicktrace = traces(3,:)' * 10;
Frequency = 1000;
threshscale = 0;
ROI = 1;
ChR = 0;
[SNR, dFF, Tconv, C2,thresh] = SNRanalysis_upright2(clicktrace, Frequency,threshscale, ROI, ChR);
spikes2 = uint16(C2(2).spikeT{1,1});

%% plot spikes
figure;
hold on
plot(pure(:,ROI),'Color','Black')
scatter(spikes,35*ones(size(spikes)),'filled','Red')
% scatter(C2.spikeT{1,1},37*ones(size(C2.spikeT{1,1})),'filled','Red')

%% filter
filtered = zeros(size(pure(:,ROI)));
filtered(spikes) = pure(spikes,ROI);
saveastiff(sub,'sub.tif');
%--run TrendFiltering in python--
denoised = double(vm('sub.tif'));
filtered(filtered == 0) = denoised(filtered == 0);

%% plot
figure;
plot(pure(:,ROI),'LineWidth',1.5)
hold on
plot(filtered,'LineWidth',0.75,'Color','Red')
scatter(spikes,15*ones(size(spikes)),'filled','Blue')
legend('Original Trace','Filtered Trace','Spikes')

