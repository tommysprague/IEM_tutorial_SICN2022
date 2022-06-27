% tcs_make_fMRI_data.m
%
% takes raw fMRI data and saves out files used for SICN 2019 IEM tutorials
%
% TC Sprague

%% save out 'fundamentals' data

delay_window = [5.25 12]*1000;
nbins = 8;

fn = '../../data/orig/CC_MGSMap1_V3AB_surf_trialData.mat';

load(fn);


% get correct bins for positions (for fundamentals tutorial)
c_all = c_map;
c_all(:,2) = ceil(c_map(:,1)/(360/nbins));
c_all(:,1) = round(c_all(:,1));


r_all = r_map;

excl_all = zeros(size(r_all)); % won't exclude any trials

tpts = which_TRs * TR;

which_tpts = tpts >= delay_window(1)*0.001 & tpts < delay_window(2)*0.001;



data_all = mean(dt_mapz(:,:,which_tpts),3);

which_vox = rf.ve >= 0.1 & std(data_all)>0;

data_all = data_all(:,which_vox);

fn_fund = '../../data/fMRI/fMRI_fundamentals.mat';
save(fn_fund,'c_all','data_all','delay_window','excl_all','r_all');



%% now save out 'advanced' dataset

dt_all = dt_mapz(:,which_vox,:);
fn_adv = '../../data/fMRI/fMRI_advanced.mat';

tpts = tpts*1000; % convert to ms to match EEG

save(fn_adv,'c_all','dt_all','excl_all','r_all','tpts');