% tcs_plot_all_subj
%
% plots each subj individually
%
% TCS (plotting data analyzed w/ SpatialEM.m)

root = '/Volumes/data/FosterEtAl2016/';

subj_to_plot = [1,2,3,5,7,8,9,10,12,13,14,15,16,18,20];

nrows = 5;
ncols = 3;

figure;
for ss = 1:length(subj_to_plot)
    
    this_CRFs = load(sprintf('%sData/%i_SpatialTF.mat',root,subj_to_plot(ss)));
    this_CRFs = this_CRFs.em; % get relevant field
    
    subplot(nrows,ncols,ss); hold on;
    thism = mean(mean(this_CRFs.tfs.total,2),4); % average over iters, blocks
    
    imagesc(this_CRFs.time,this_CRFs.cCenters-180,squeeze(thism).');
    title(sprintf('Subj %02.f',subj_to_plot(ss)));
    
end