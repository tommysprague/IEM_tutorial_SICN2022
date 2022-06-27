% iem_sim_simple
% written: js (jserences@ucsd.edu) and ts (tsprague@ucsb.edu), 05.2019
%
% NOTE: this script is available at github.com/JohnSerences/iem_sim - it is
% not striclty a tutorial, but demonstrates some simple simulation
% principles, as well as some features of 'fixed' encoding model analyses.
% Moreover, it implements invertible linear transforms of channel models,
% as described in Gardner & Liu, 2019, eNeuro.
%
% implements a simple neuron-to-voxel-to-IEM simulation similar to the
% simulations in Gardner and Liu, 2019 eNeuro (GandL2019). 
% 
% To keep the linear transform part of this simulation simple and to improve readability, 
% the number of features (nf) is always fixed at 8, as is
% the number of 'channels' (nc) (after GandL2019). See in-line comments for more
% on this and see references below for code that can actually be used to analyze
% fMRI/EEG data where the constraint is only that nf>nc.
%
% see here for the more practical tutorial code: 
% https://github.com/tommysprague/IEM-tutorial (Matlab)
% https://github.com/tommysprague/IEM_tutorial_SICN2019 (Matlab)
% https://github.com/JohnSerences/encoding_models (Python)
% 
% In this simulation, each voxel contains a bunch of simulated neurons,
% with each neuron having a von Mises tuning function with randomly generated tuning preference mu and 
% concentration k (i.e. tuning width).
%
% After voxel responses are generated, the simulation will apply a forward
% model that consists either of raised cos basis functions, stick
% functions (delta), or linearly transformed (xformed) versions of these standard basis functions so
% that they become bimodal. The simulation will then fit the forward model
% (however specified), invert the model to recover channel responses, and
% then either multiply the resulting channel response profiles by inv(xform) or not...
% 
% Unlike Gardner and Liu, this simulation will implement a two-condition
% experiment in which there is attentional gain on the neural responses in
% condition 2 compared to condition 1 (note this could be contrast-induced
% gain, or whatever...its just a scalar applied to stimulus-driven 
% neural population response profile)
%
% At the end of the simulation, a simple linear regression will be
% performed to quantify the difference in response amplitude between
% conditions. This is not suggested as a real final analysis for this type of data,
% especially in the case of funky basis functions. Instead, this is just meant to
% provide a quick summary of the estimated amplitude of the
% difference between conditions. 
%
% sample settings...
% Setting 1: 
% use_cos_or_stick = 1;
% use_xform = 0; 
% inverse_before_shift = 0;
% will analyze the data with typical cos^x basis set
%
% Setting 2: 
% use_cos_or_stick = 1;
% use_xform = 1; 
% inverse_before_shift = 0;
% will analyze with bimodal basis in forward model, and after inversion, channel
% response profile will also tend to be bimodal. Importantly, the difference 
% in amplitude between condtions is well-preserved, even though a bimodal basis set was used.
% Note that with a lot of noise, quantifying the amplitude (and
% baseline/shape) of bi-modal channel response profiles is more difficult, and 
% some small differences in amp/baseline of the responses may appear. But
% comparisons between condition 1 and condition 2 (esp amplitude minus baseline) yield similar 
% results as Setting 1 across a wide range of noise values.  
% 
% Setting 3: 
% use_cos_or_stick = 1;
% use_xform = 1; 
% inverse_before_shift = 1;
% will analyze with bimodal basis in forward model, but will then multiply
% channel responses by inv(xform). Will produce same differences between
% condition 1 and condition 2 as Setting 1
% even though the forward model used bimodal as opposed to cos^x basis funcs. 
% Note also that because in both cases we're quantifying cos-like channel response profiles 
% when assessing the differences between condition 1 and condition 2, that 
% the results of this analysis will be the same as under Setting 1 across noise levels.
% 
% You can also set use_cos_or_stick to 0 to see similar patterns of
% condition-wise differences given xformed and non-xformed basis sets when
% stick basis functions are used in place of cos^x

clear
close all

%% number of stimulus features and number of channels in forward model
% note that features are currently spaced evenly over a 360 deg feature
% space - for the moment set these to 8 and 8 and in this code and leave these
% as fixed values so that the xform steps work in a straightforward manner. 
% To actually implement this analysis on real data, nc and nf do not have to be the same, so long 
% as nf>nc. For more practical guides to implementation on real data, see:
% https://github.com/tommysprague/IEM-tutorial (Matlab)
% https://github.com/JohnSerences/encoding_models (Python)
% This code simply meant to emulate a version of the GandL2019
% simulation with nf==nc==8. 

nf = 8;    % number of features in experiment (orientations, direction, locations, etc) - note: will define over 2*pi by default
nc = 8;    % number of channels in basis set

%% use cos (1) or stick functions (0) for forward model
use_cos_or_stick = 1;

%% apply linear transform to make a bimodal basis? 1=yes, 0=no
use_xform = 0; 

%% multiply channel responses by inv(xform) before shifting/averaging/fitting? this 
% will be applied to channel responses after model inversion and will undo the effects of 
% the xform on the original basis (akin to GandL2019 Figure 3)
inverse_before_shift = 0;

%% set the gain for condition differences here. In this simulation, only have 
% a scalar gain parameter. 
attn_gain = 1.8;

%% scalar for setting randn noise added to the simulated data set - added at 
% the level of simulated neurons - keep in mind that this noise should be
% considered in the context of each neurons min and max response that is
% set below.
n = 0;     

%% Mixing vec for generating bi-modal basis (xform matrix). The 
% nc x nc xform matrix is made by circularly shifting this vector. Note that 
% this imposes a constraint on the hard-coded nc=8...could replace this
% with a function to generate a mixing-vec of length(nc) if you want a more
% flexible simulation. 
mixing_vec = 0.8 * [0 1 0.95 0 0.05 0 0.95 1];

%% number of trials in simulated experiment (must be evenly divisible by nf)
% 1/2 of the trials will be in condition 1, and 1/2 in condition 2, and
% furthermore a balanced 1/2 will be used as the training set and the other
% 1/2 as the test set. 
ntrials = 256;       

%% how many voxels to simulate? 
nvox = 50;              

%% number of 'neurons' in each voxel - doesn't matter too much so long as it's 
% reasonable, but will slow things down quite a bit if this gets too large.
% 50-100 seems to work fine and finished pretty quickly on most machines. 
nn = 100;

%% max response of neurons in each voxel, and range of max responses. if 
% max_neural_resp_range == 0, then max of all neural tfs peak at max_neural_resp, if
% max_neural_resp_range > 0 then max neural tfs given by: 
% max_neural_resp + 2*(rand-.5)*max_neural_resp_range
% note that these values are not meant to be super meaningful, just to get
% responses off of 0...
max_neural_resp = 30;         % max of neural tfs
min_neural_resp = 5;          % min of neural tfs
max_neural_resp_range = 0;          

%% range of concentrations for neural tuning functions - if [x,x] then all 
% neurons will have concentration set by x, if [x,y], then random uniform
% distribution of concentration parameters over the range x:y.
% randomly within each voxel across the range x:y
k_range = [3,8];   % modest range, note that bigger k == narrower tuning
%k_range = [1 14]; % bigger range
%k_range = [4,4];  % shared value for all neurons
           
%% raised cos for making design matrix of forward
% model - specify xx and mu in radians.
make_cos_basis = @(xx,mu) max(0, cos(xx - mu)).^(nc-mod(nc,2));

%% circular gaussian (vm) - used for neural tuning 
% - use radians for xx and mu...use this as opposed to make_cos_basis just for 
% historical reasons...can adjust bandwidth with k (k->bigger == narrower
% tuning)
make_vm = @(xx,mu,k) exp(k.*(cos(xx-mu)-1));

%% main section for generating neural tuning and then voxel responses
% first define the tuning preference of 'neurons' in each simulated voxel
% Then compute summed neural response in each voxel and add noise to response on each trial, and
% gain per condition as appropriate. Note that gain is scalar 

% stim features - defined over a 2*pi circular space 
stim_x = linspace(0, 2*pi-(2*pi/nf), nf); 

% neural tuning function x-axis - 2*pi radian space to define neural tfs. 
neural_x = linspace(0, 2*pi-(2*pi/360), 360);

% make a list of the feature on each trial - keep balanced for now...
trial_features = repmat(stim_x', ntrials/nf, 1); 

% make a list of trial conditions - first nf trials in cond1, next nf
% trials in cond2, repeat...
trial_conds = repmat([ones(nf,1);ones(nf,1)*2],ntrials/(nf*2),1);

% set up scalar gain for each trial - vector of ones for cond1 trials, and
% attn_gain for cond2 trials. Since the trials are organized as 1:nf ==
% cond1, nf+1:nf*2 are cond2, etc, can stack nf ones on top of nf attn_gain
% and then repmat that to assign gain to each trial in the simulation
trial_gain = repmat([ones(nf, 1); repmat(attn_gain, nf, 1)], ntrials/(nf*2), 1);

% figure out neural tfs for each voxel.
% random tuning of neurons in each voxel...this is a 
% nvox x nn matrix of tuning preferences for the neurons in each voxel 
n_prefs = rand(nvox,nn) * (2*pi);

% compute the actual neural response on each trial and generate voxel responses...
% start loop over voxels...generate a #trial x #vox matrix, just like we normally use for analysis of real data. 
sim_vox_data = nan(ntrials,nvox);

for vox = 1:nvox

    % generate a vector k to govern concentration of neural tuning functions.
    % If k_range = [x,x] or k_range = [x], then k==x for every voxel
    k = rand(1,nn)*(k_range(end)-k_range(1)) + k_range(1); 

    % neural tuning functions - will return a nn x 360 matrix, with a
    % neural tuning function centered on n_prefs(vox,neuron) and with concentration 
    % k in each row...doing this without looping over #neurons
    n_tfs_temp = make_vm(repmat(neural_x,nn,1), repmat(n_prefs(vox,:)',1,360), repmat(k',1,360));
    
    % then scale the n_tfs...can all have same max resp (if max_neural_resp_range==0) 
    % or different if max_neural_resp_range>0...
    max_resp_in_each_neuron = repmat(max_neural_resp + 2*(rand(nn,1)-.5)*max_neural_resp_range, 1, 360);
    n_tfs = (n_tfs_temp .* (max_resp_in_each_neuron-min_neural_resp) + min_neural_resp);
    
    % now compute the response on every trial for every voxel
    % add gain and noise as specified...
    for trial = 1:ntrials

        % then find the index into the columns of 
        % the n_tfs matrix that most closely corresponds to the stimulus presented 
        % on the present trial (this will also support a full range of stimulus
        % features that take on random values in future versions of simulation 
        % if desired). 
        [~,feat_ind] = min(abs(neural_x-trial_features(trial)));
        
        % grab column out of n_tfs to generate stimulus-driven response on current 
        % trial in current voxel...then add noise, then 
        % apply gain to n_tfs...note that noise added here before gain! 
        n_tfs_w_noise = n_tfs(:, feat_ind)+(randn(nn,1)*n);
        
        % then multiply by scalar gain factor on current trial
        n_tfs_w_noise_gain = trial_gain(trial) * n_tfs_w_noise;        
        
        % get the vector of responses, take sum as voxel response on this
        % trial
        sim_vox_data(trial,vox) = sum(n_tfs_w_noise_gain);
        
    end 
end

%% use half the data for training
% this will be balanced across condition 1 and condition 2
trn  = sim_vox_data(  1:(ntrials/2),:);
trnf = trial_features(1:(ntrials/2),:);

%% the other half for testing
tst  = sim_vox_data(  (ntrials/2)+1:end, :);
tstf = trial_features((ntrials/2)+1:end, :);

% vector for labeling conditions in the tst data 
tst_cond = trial_conds((ntrials/2)+1:end);

%% make design matrix 
channel_x = linspace(0, 2*pi-(2*pi/nc), nc); 
X = zeros(size(trn,1), nc);

% use stick (delta) or cos basis functions
if use_cos_or_stick==1
    % cos basis
    for i = 1:size(trn,1)
        X(i,:) = make_cos_basis(channel_x, trnf(i));
    end

elseif use_cos_or_stick==0
    % stick functions
    for i = 1:size(trn,1)
        X(i,trnf(i)==stim_x)=1;
    end
    
else
    disp('error in dm construction')
    return
end

% normalize to have unit amp
X = X./max(X(:)); 

%% linear transform of X? 
if use_xform==1
    mix_mat = zeros(nc,nc);
    % circshift...
    for i = 1:nc
        mix_mat(i,:) = circshift(mixing_vec, [0,(i-1)]);
    end 
    X = X*mix_mat;
end

%% compute weights (forward model)
w = X\trn;

%% invert to compute channel responses
% cr = (inv(w*w')*w*tst')'; 
cr = (w.'\tst.').';

%% then shift the rows to recenter channel response profiles

% apply inverse transform before circ shifting? 
if use_xform==1 && inverse_before_shift==1
    cr = cr * inv(mix_mat); 
end

% figure out how much to shift each profile by to 'center' it (although
% with xformed basis, not really sure so much where the center is anymore,
% but can shift it anyway to a common point at least). 
shiftby = (nf/2)+1;

% do the shifting
for ii=1:size(cr,1)
   cr(ii,:) = circshift(cr(ii,:), shiftby-find(tstf(ii)==stim_x), 2);
end

%% fitting the difference between conditions, just as a check to see if the 
% magnitude of the difference is preserved across transforms. 

% compute difference 
d = mean(cr(tst_cond==2,:))-mean(cr(tst_cond==1,:));

% build a dm for regression to estimate amp and baseline offset of
% difference - grab a row out of the main design matrix X and use that,
% plus a column of ones

if use_xform==1 && inverse_before_shift==1 
    % invert the dm X so we can grab a row
    nX = X*inv(mix_mat);  % assign to nX for plotting (so we can preserve the original X)
else
    nX = X; 
end

% plots...
marksersz = 15; 
lw = 4; 
fsz = 30; 
px = -180:360/nc:180-(360/nc);
 
% make some plots
figure(1), clf, hold on
%title('Cond1 vs. Cond2')
plot(px, mean(cr(tst_cond==1,:)), 'o-', 'markersize', marksersz, 'linewidth', lw, 'MarkerFaceColor', 'auto')
plot(px, mean(cr(tst_cond==2,:)), 'o-',  'markersize', marksersz, 'linewidth', lw, 'MarkerFaceColor', 'auto')
set(gca, 'XTick', px,'TickDir','out')
set(gca, 'XLim', [-180, 180])
set(gca, 'XTickLabel', {'-180','','-90','','0','','90',''})
set(gca, 'YLim', [-.1, 1.5])
set(gca, 'YTick', 0:.5:1.5)
xlabel('Feature difference from true')
ylabel('Channel Response (a.u)')
%legend({'Cond1', 'Cond2'}, 'FontSize', fsz, 'Location', 'NorthEast')
set(gca, 'FontSize', fsz)

figure(2), clf
%title('Difference (Cond2 - Cond1)')
plot(px, mean(cr(tst_cond==2,:))-mean(cr(tst_cond==1,:)), 'o-', 'markersize', marksersz, 'MarkerFaceColor', 'auto', 'linewidth', lw)
set(gca, 'XTick', px,'TickDir','out')
set(gca, 'XLim', [-180, 180])
set(gca, 'XTickLabel', {'-180','','-90','','0','','90',''})
set(gca, 'YLim', [-.1, .8])
set(gca, 'YTick', 0:.4:.8)
xlabel('Feature difference from true')
ylabel('Channel Response (a.u)')
set(gca, 'FontSize', fsz)

figure(3), clf, hold on
%title('Neural tuning in a voxel')
plot(neural_x*180/pi, n_tfs(1:10,:)', '-', 'markersize', marksersz, 'MarkerFaceColor', 'auto', 'linewidth', lw)
%plot(neural_x*180/pi, n_tfs([2,21,44],:)', '-', 'linewidth', lw*4)
tmpx = linspace(0,360-360/nf,nf/2);
set(gca, 'XLim', [0, 360],'XTick',neural_x(1:90:end)*180/pi,'TickDir','out');
set(gca, 'YTick',[])
xlabel('Feature value')
ylabel('Firing Rate')
set(gca, 'FontSize', fsz)

% plot a continuous version of the basis set for making figure...

% if cosine basis...
if use_cos_or_stick == 1
    nth = 360; myth = linspace(360/nth,360,nth); chth = linspace(0,360-360/nf,nf);
    tmpX = make_cos_basis(deg2rad(myth.'),deg2rad(chth));
    
    if use_xform==1
        contX = tmpX*mix_mat;
    else
        contX = tmpX;
    end
    
    figure(4), clf, hold on
    %title('Basis set (continuous)')
    plot(myth,contX,'-','LineWidth',lw/2);
    plot(myth,contX(:,5),'-','LineWidth',lw*2);
    set(gca, 'XLim', [0, 360],'XTick',chth(1:2:end),'TickDir','out');
    xlabel('Feature value')
    ylabel('Channel Sensitivity')
    set(gca, 'FontSize', fsz)
    
end