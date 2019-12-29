function [Mr]=NoRMCorre2(mov,output)
gcp;
%% read data and convert to double
Yf = single(mov);
[d1,d2,T] = size(Yf);

%% perform some sort of deblurring/high pass filtering

if (0)    
    hLarge = fspecial('average', 40);
    hSmall = fspecial('average', 2); 
    for t = 1:T
        Y(:,:,t) = filter2(hSmall,Yf(:,:,t)) - filter2(hLarge, Yf(:,:,t));
    end
    %Ypc = Yf - Y;
    bound = size(hLarge,1);
else
    gSig = 7; 
    gSiz = 17; 
    psf = fspecial('gaussian', round(gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    Y = imfilter(Yf,psf,'same');
    %Ypc = Yf - Y;
    bound = 2*ceil(gSiz/2);
end
%% first try out rigid motion correction
    % exclude boundaries due to high pass filtering effects
options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50,'max_shift',20,'iter',1,'correct_bidir',false);

%% register data and apply shifts to removed percentile
tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
    % exclude boundaries due to high pass filtering effects
tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset
    % apply shifts on the whole movie
try
    %% compute metrics 
[cY,mY,vY] = motion_metrics(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
[cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);

[cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
[cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);

%% plot rigid shifts and metrics
shifts_r = squeeze(cat(3,shifts1(:).shifts));
figure;
    subplot(311); plot(shifts_r);
        title('Rigid shifts','fontsize',14,'fontweight','bold');
        legend('y-shifts','x-shifts');
    subplot(312); plot(1:T,cY,1:T,cM1);
        title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
        legend('raw','rigid');
    subplot(313); plot(1:T,cYf,1:T,cM1f);
        title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
        legend('raw','rigid');
        saveas(gca,fullfile(output,'RigidCorrectionMetrics.fig'));
catch
end
%         %% now apply non-rigid motion correction
% % non-rigid motion correction is likely to produce very similar results
% % since there is no raster scanning effect in wide field imaging
% 
% options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
%     'grid_size',[128,128],'mot_uf',4,'correct_bidir',false, ...
%     'overlap_pre',32,'overlap_post',32,'max_shift',20);
% 
% tic; [M2,shifts2,template2] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_nr,template1); toc % register filtered data
% tic; Mpr = apply_shifts(Yf,shifts2,options_nr,bound/2,bound/2); toc % apply the shifts to the removed percentile
% 
% %% compute metrics
% 
% [cM2,mM2,vM2] = motion_metrics(M2,options_nr.max_shift);
% [cM2f,mM2f,vM2f] = motion_metrics(Mpr,options_nr.max_shift);
% 
% %% plot non rigid shifts        
% 
% shifts_r = squeeze(cat(3,shifts1(:).shifts));
% shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
% shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
% shifts_x = squeeze(shifts_nr(:,2,:))';
% shifts_y = squeeze(shifts_nr(:,1,:))';
% 
% patch_id = 1:size(shifts_x,2);
% str = strtrim(cellstr(int2str(patch_id.')));
% str = cellfun(@(x) ['patch # ',x],str,'un',0);
% 
% figure;
%     ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients for filtered data','fontsize',14,'fontweight','bold')
%             set(gca,'Xtick',[],'XLim',[0,T-3])
%     ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
%             set(gca,'Xtick',[])
%     ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
%             xlabel('timestep','fontsize',14,'fontweight','bold')
%     linkaxes([ax1,ax2,ax3],'x')
% saveas(gca,'NonRigidCorrectionMetrics.fig')
end