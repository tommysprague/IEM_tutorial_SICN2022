% IEM_tutorial_advanced_full.m
%
% Originally presented at cuttingEEG 2018 in Paris, France
% Adapted for fMRI, presented at Summer Institute for Cognitive
% Neuroscience, UC Santa Barbara, July 4 2019
% Tommy Sprague (tsprague@ucsb.edu)
% Asst Prof, Dept Psychological & Brain Sciences, UC Santa Barbara
%
% Now that you've spent some time learning the fundamentals of how to run an IEM
% analysis on a small set of data, let's explore some more advanced aspects
% of the analysis you'll need to understand if you want to run an IEM-based
% experiment yourself.
%
% like the "fundamentals" tutorial, you'll have the option of working with
% fMRI or EEG data. For EEG, there will be a few extra components (marked
% as "EEG only") that cover best practices for selecting relevant signals.
% For Brain Camp, we'll only cover the fMRI parts, but feel free to use
% these materials to help design EEG experiments/analyses as well
%
% This tutorial will cover:
%
% - kinds of EEG data to use: focus on time course of alpha power (EEG
%   only)
%
% - run-wise cross-validation: compute reconstructions for each trial by
%   looping over runs (leave-one-run-out). We'll also discuss other
%   cross-validation schemes.
%
% - balancing trials within CV folds - due to artifact rejection, trial
%   censoring based on motion correction, etc, you sometimes end up with
%   lopsided #'s of trials over feature space. I'll show an example for how
%   to address this, and we'll see how much it impacts the data quality
%
% - Converting channel responses into 'reconstructions' - we'll weight the
%   channel sensitivity profile for each channel by its estimated
%   activation. This is a 'smoothed' version of the channel response
%   function (CRF), and can be 'exactly' coregistered across trials (not
%   just based on bin). This will also let us more fully explore how
%   different basis sets impact our analyses
%
% - Quantifying reconstructions based on their 'fidelity' (Sprague et al,
%   2016) and/or their 'slope' (Foster et al, 2016), the two most
%   commonly-used metrics right now. We'll also discuss when it's
%   appropriate to do curve-fitting to quantify 'amplitude' and 'size'
%
% - Computing reconstructions over time to uncover neural representations
%   at the timescale of the measurements - ms (EEG) or ~s (fMRI). We'll go
%   over a few ways to do this, and discuss how to compare results from
%   each method.
%
% - Understanding why it's important to compare reconstructions across
%   conditions (or timepoints) computed using a 'fixed' encoding model
%   (Sprague et al, 2018)
%
% - Evaluating the 'significance' of reconstructions: we'll generate a null
%   distribution of reconstructions against which we'll compare our actual
%   data. This is an important step to ensure any structure in the data is
%   not due to chance. We'll discuss these procedures at the group level as
%   well.
%
% Feel free to share this tutorial and the data associated with it (which
% comes from the OSF repositiory for Expt 1 of Foster et al, 2016:
% https://osf.io/q6sxh/). If you found it helpful, it can never hurt to
% cite either the GitHub repsository (TODO), and/or our recent commentary
% on these methods (Sprague et al, 2018, eNeuro:
% http://www.eneuro.org/content/5/3/ENEURO.0098-18.2018)
%
% If you have any questions, concerns, ideas, or want advice on these
% methods, I'm always available and happy to help - tsprague@ucsb.edu
%
% And stay tuned to my github (github.com/tommysprague) - I'm planning on
% building a set of IEM utilities to help make some common procedures
% quicker/easier to implement. This will *not* be a toolbox, as I believe
% it's important to understand the steps of your analyses at a deep level -
% but it will provide recipes and outlines for common analyses (similar to
% below), at both the individual subj and group level.

% do you want to use EEG ('EEG') or fMRI ('fMRI') data?
which_data = 'fMRI';  % for SICN: focus on fMRI (but will work for EEG too)

% IMPORTANT: set this if you are using an older version of matlab that does not
% have new features (if you get errors, try this first)
old_matlab_ver = 0;

% a few 'global' variables
pbins = 0:45:(360-45); % center of position bins


if strcmpi(which_data,'EEG')
    twindow = [-500 2000]; % the time window we'll consider for analyses
    sampling_rate = 250; % recorded eeg data at 250 Hz
    delay_window = [750 1250]; % for delay-period comparisons, which tpts?
elseif strcmpi(which_data,'fMRI')
    twindow = [-2250 12750]; % time window we'll consider for analyses (ms)
    sampling_rate = 1/.75; % TR = 750 ms (sampling rate in Hz)
    delay_window = [5250 12000]; % ms
end


%% first - load data
addpath util/;
load(sprintf('../data/%s/%s_advanced.mat',which_data,which_data));

% in your workspace, you'll have:
% - dt_all: full timeseries for each signal dimension (electrode or voxel) on each trial (n_trials x
%   n_signaldim x n_tpts)
% - c_all:  condition label for each trial. first column is the polar angle 
%   remembered on that trial, in degrees (Cartesian coords, 45 = upper right); second column is the position bin (1-8; beginning at 0 deg and moving CCW, so bin 2 is centered at 45, 3 at 90, etc); n_trials x 2
% - excl_all: logical vector, contains 1 for each trial that was marked by
%   Foster et al (2016) for exclusion based on artifacts, etc; n_trials x 1
%   (all 0's for fMRI data)
% - r_all: label for each trial as to which run it came from; n_trials x 1
% - chan_labels: EEG only - cell array of strings labeling each electrode; 
%   n_electrodes+2 x 1 cell array (we don't include EOG channels in data here)
% - tpts: time at each sample in ms (size(dt_all,3) x 1)

% let's look at the 'raw' data - pick a channel and sort trials by position
% bin and plot average over all trials; look at mean delay period response
% vs position bin

dim_to_plot = 7;  % which electrode/vox? see chan_labels

delay_tpts = tpts >= delay_window(1) & tpts < delay_window(2); % which tpts we'll use for 'delay' analyses

figure; subplot(1,2,1); hold on;
pu = unique(c_all(:,2)); % get all position indices
pos_colors = hsv(length(pu)); tmp_mu = nan(length(pu),1); tmp_std = nan(length(pu),1);
for pp = 1:length(pu)
    thisidx = c_all(:,2)==pu(pp) & ~excl_all; % which trials we average
    plot(tpts,squeeze(mean(dt_all(thisidx,dim_to_plot,:),1)),'Color',pos_colors(pp,:));
    tmp_mu(pp) = mean(mean(dt_all(thisidx,dim_to_plot,tpts>=delay_window(1)&tpts<=delay_window(2)),3),1);
    tmp_std(pp) = std(mean(dt_all(thisidx,dim_to_plot,tpts>=delay_window(1)&tpts<=delay_window(2)),3),[],1);
end
xlim(twindow);
xlabel('Time (ms)');
ylabel('Response');
title(sprintf('Dimension %i',dim_to_plot));
set(gca,'TickDir','out','FontSize',14);

subplot(1,2,2);
hold on; plot(pbins,tmp_mu,'k-');
for pp = 1:length(pu)
    plot(pbins(pp)*[1 1],tmp_mu(pp)+[-1 1]*tmp_std(pp),'-','Color',pos_colors(pp,:));
    plot(pbins(pp),tmp_mu(pp),'o','MarkerSize',10,'Color','w','MarkerFaceColor',pos_colors(pp,:))
end
xlabel('Position bin (center, \circ)');
ylabel('Mean delay response');
title(sprintf('%i to %i ms',delay_window(1),delay_window(2)));
xlim([-45 360]);
set(gca,'TickDir','out','FontSize',14);

% Try looking at a few signal dimensions (voxels/electrodes) to see if you
% can get a sense for the stability of the delay-period potential as a
% function of stimulus position.



% NOTE: we're skipping this step in the brain camp tutorial! if you are
% looking at EEG data, this will happen behind-the-scenes, but you can set
% force_plot=1 to see the results. This step is not required for fMRI
% processing

if strcmpi(which_data,'EEG')
    
    %% filter EEG data
    
    % Based on our quick look at a few example electrodes, it doesn't seem
    % promising to use 'raw' EEG data (though, if you like, feel free to try it
    % below!). We do know that the spatial distribution of alpha (8-12 Hz)
    % power on the scalp reflects changes in locus of spatial attention - so it
    % was deemed a good candidate band for looking at spatial WM
    % representations as well (see Foster et al, 2016 for analyses of other
    % bands as well; or play w/ the filter properties below).
    
    % set up our filter properties - we're going to use the simple 'filter,
    % then Hilbert' method to get timecourses of power at a given frequency
    % band. Wavelets, etc, work equivalently well. In my experience, power is
    % all that matters - though using the real & imag parts of Fourier/wavelet
    % coefficients can also be useful (as [real(coeffs) imag(coeffs)] ).
    % I'd stay away from cicular quantities like angle(coeffs) though.
    
    
    filt_band = [8 12]; % filter from 8-12 Hz
    
    % NOTE: you could set this up to look at different freqs - left as exercise
    % for reader...
    
    % use EEGFILT to band-pass the data, then get instantaneous power via
    % Hilbert transform
    
    dbp_all = nan(size(dt_all)); % band-passed data
    df_all = nan(size(dt_all));  % freq data
    for cc = 1:size(dt_all,2)
        % band-pass data using eeg_filt
        dbp_all(:,cc,:) = eegfilt(squeeze(dt_all(:,cc,:)),sampling_rate,filt_band(1,1),filt_band(1,2));
        
        % hilbert-transform and compute power (abs(H(X)).^2)
        df_all(:,cc,:)  = abs(hilbert(eegfilt(squeeze(dt_all(:,cc,:)),sampling_rate,filt_band(1,1),filt_band(1,2)).').').^2; % instantaneous power calculated here for induced activity.
    end
    
    % if we're using EEG data, our "data" (d_all) will be the filtered data
    % (df_all)
    d_all = df_all;
    
    % if you're curious about how data filtering works for EEG, set this to 1
    force_plot = 0;
    
    %% inspect filtered data
    if force_plot == 1
        
        % let's look at an example trial, channel:
        tnum = 51; chnum = 7;
        figure;
        hold on;
        % plot raw signal
        plot(tpts,squeeze(dt_all(tnum,chnum,:)),'LineWidth',1);
        
        xlabel('Time (ms)');
        ylabel('Signal');
        title(sprintf('Trial %i, Channel %s',tnum,chan_labels{chnum}));
        set(gca,'TickDir','out','FontSize',14);
        
        % plot band-pass signal
        plot(tpts,squeeze(dbp_all(tnum,chnum,:)),'LineWidth',1);
        % plot abs(Hilbert-transform) - envelope
        plot(tpts,squeeze(d_all(tnum,chnum,:)).^0.5,'LineWidth',1);
        % plot power from Hilbert
        plot(tpts,squeeze(d_all(tnum,chnum,:)),'LineWidth',1);
        
        legend({'Raw (\muV)','Filtered (\muV)','Envelope (\muV)','Power (\muV^2)'});
        
        
        %% Sort filtered data by position
        % Similar to above, let's look at the same channel sorted by position
        figure; subplot(1,2,1); hold on;
        tmp_mu = nan(length(pu),1); tmp_std = nan(length(pu),1);
        for pp = 1:length(pu)
            thisidx = c_all(:,2)==pu(pp) & ~excl_all; % which trials we average
            plot(tpts,squeeze(mean(d_all(thisidx,dim_to_plot,:),1)),'Color',pos_colors(pp,:),'LineWidth',1);
            tmp_mu(pp) = mean(mean(d_all(thisidx,dim_to_plot,tpts>=delay_window(1)&tpts<=delay_window(2)),3),1);
            tmp_std(pp) = std(mean(d_all(thisidx,dim_to_plot,tpts>=delay_window(1)&tpts<=delay_window(2)),3),[],1);
        end
        xlim(twindow);
        xlabel('Time (ms)');
        ylabel(sprintf('%i to %i Hz power (\\muV^2)',filt_band(1),filt_band(2)));
        title(sprintf('Electrode %s',chan_labels{dim_to_plot}));
        set(gca,'TickDir','out','FontSize',14);
        
        subplot(1,2,2);
        hold on; plot(pbins,tmp_mu,'k-');
        for pp = 1:length(pu)
            plot(pbins(pp)*[1 1],tmp_mu(pp)+[-1 1]*tmp_std(pp),'-','Color',pos_colors(pp,:));
            plot(pbins(pp),tmp_mu(pp),'o','MarkerSize',10,'Color','w','MarkerFaceColor',pos_colors(pp,:))
        end
        xlabel('Position bin (center, \circ)');
        ylabel('Mean delay power (\muV^2)');
        title(sprintf('%i to %i ms',delay_window(1),delay_window(2)));
        set(gca,'TickDir','out','FontSize',14);
        xlim([-45 360]);
    end
    
else
    
    % if we're using fMRI data, our data will be the original data source
    % (dt_all)
    
    d_all = dt_all;
    
end



%% Step 1: Build spatial IEM (same as in fundamentals tutorial!)
%
% Motivated by the spatial sensitivity of delay-period responses in many
% signal dimensions, let's build a spatial IEM, which projects measured
% responses from signal space into 'information space' (polar angle
% coordinates), so that we can directly average channel response profiles
% across all trials.
%
% This is overall the same procedure in fMRI and EEG, though SNR issues may
% require some adjustments to analysis. A few valid analysis choices
% involve:
%
% - selecting 'selective' signal dimensions based on training data, or using
%   reported EEG channels in the literature (e.g., Wolff et al, 2017)
%
% - ensuring perfectly equated representation of each location bin in
%   training set: due to artifact rejection, it may be necessary to perform
%   model estimation using average exemplars within position bins rather
%   than the actual positions. Out of an abundance of caution, it would
%   also be necessary to ensure the same # of trials contribute to each
%   position bin in the training set, and shuffle this a few times to be
%   sure results don't come from accidental trial groupings
%
% - training on 'average' delay period response, then test on each timepoint


% this will proceed just like the 'fundamentals' tutorial - so let's first
% define the properties of our basis set:

n_chan = 8; % # of channels, evenly spaced around the screen
chan_centers = linspace(360/n_chan,360,n_chan);

chan_width = 180; % default value - try making wider/thinner

% evaluate basis set at these
angs = linspace(-179,180,360);


% I moved the code used to create channel profiles into a separate
% function: build_basis_polar_mat (in util/). It also automatically changes
% the power the channels are raised to based on the number of channels,
% unless you specify your own (typically 8). This has roots in some old
% papers (Freeman & Adelson, 1992) which describe properties of 'steerable
% filters'. The basic idea is that a steerable filter is one that
% can be "steered" by linear weighting to any position. That is, each
% of these basis functions has a peak at a different polar angle.
% If you want to make the exact same shaped basis function to have a peak
% at a different polar angle, for steerable filters, there is a linear
% combination of these basis functions that will do so. A simple
% example of this property is sine and cosine functions. You can
% make a sinusoid of any phase by properly weighting a sine and cosine
% function. Turns out this is true for any power of sine and cosine
% as well. You need n+1 exponent raised sinusods where n is the 
% power you are raising to for a steerable basis set.
myb_orig = build_basis_polar_mat(angs,chan_centers,chan_width);

% now let's look at the basis set:
figure; plot(angs,myb_orig,'LineWidth',1.5);
xlabel('Polar angle (\circ)');ylabel('Channel sensitivity'); title('Basis set (information channels)');
xlim([-180 180]);set(gca,'XTick',-180:90:180,'TickDir','out','Box','off','FontSize',14);


%% Step 2a: Use basis set to compute channel responses


% just like in the 'fundamentals' tutorial, now we need to use our
% condition labels on each trial to compute channel responses, given the
% channel sensitivity profiles defined by our basis set (above).

% Use our c_all variable, which contains the exact angle viewed on
% each trial to compute channel responses
%
% start by making a matrix n_trials x 360, with 1's and 0's, which we can
% multiply by our basis_set (b_orig) to generate computed channel
% responses

% first, we'll need to wrap around angles greater than 180 to -179
c_all(c_all(:,1)>180,1) = c_all(c_all(:,1)>180,1)-360;

stim_mask = zeros(size(c_all,1),length(angs));
for tt = 1:size(c_all,1)
    % for now, we're going to assume the representation 'should be' perfect
    % - the feature mask we're projecting into channel space will be a
    % delta function. it's also possible to try a 'blob' here too.
    stim_mask(tt,angs==c_all(tt,1)) = 1;
    clear idx;
end

% HINT: here's a faster way to do this w/ newer versions of matlab:
% stim_mask = c_all(:,1)==angs;

% compute channel responses for each trial
% (when using a delta stimulus, this is equivalent to shifting a basis
% function and evaluting it at the channel centers - but sometimes your
% model may be more complicated than a delta stimulus, so it's good to be
% in the habit of writing it like this)
X_all = stim_mask * myb_orig;


% let's check out the computed channel responses:
figure;
% first, a single trial
whichtrial_C = 125;
subplot(1,2,1);
hold on;  chan_colors = lines(n_chan);
plot(chan_centers,X_all(whichtrial_C,:),'k-','LineWidth',1.5);
for cc = 1:n_chan
    plot(chan_centers(cc),X_all(whichtrial_C,cc),'ko','MarkerSize',10,'MarkerFaceColor',chan_colors(cc,:))
end
plot([1 1]*mod(c_all(whichtrial_C),360),[0 1],'k--','LineWidth',2);
xlabel('Channel center (\circ)');
ylabel('Predicted channel response');
title(sprintf('Trial %i: %i\\circ',whichtrial_C,mod(c_all(whichtrial_C,1),360)));
set(gca,'TickDir','out','XTick',0:90:360,'FontSize',14);


% all trials
subplot(1,2,2); hold on;
imagesc(chan_centers,1:size(c_all,1),X_all); axis ij;
plot([0 360],whichtrial_C+[-0.5 0.5],'r-');
xlim([22.5 360+22.5]);ylim([0.5 size(c_all,1)+0.5]);
xlabel('Channel center (\circ)');
ylabel('Trial');

% put the rank of the design matrix (X_all) in the title
title(sprintf('All trials (rank = %i)',rank(X_all)));
set(gca,'TickDir','out','XTick',0:90:360,'FontSize',14);


%% Step 2b & 3: Train/test IEM (full delay period) - leave-one-run-out
%
% With our encoding model and predicted channel responses in hand, now we
% can estimate channel weights and reconstruct WM representations from each
% trial. In the fundamentals tutorial, we just split data in half and
% trained using one half and tested with the other. But that of course left
% half the data non-reconstructed! Here we implement a n-fold
% cross-validation procedure (common in decoding analyses) in which we
% split the data into n parts, loop over the parts and hold one part out as
% a 'test' set and use the remaining data as a 'training' set for model
% estimation. For simplicity, we'll use the run number as the partition.
% (note: for fMRI data, if any within-run normalization, like detrending or
% Z-scoring is applied, to be safe it's best to avoid training/testing
% within the same run if possible)
%
% Since we only have one task here, we'll also perform
% leave-one-run-out cross-validation (to start with). We'll also explore a
% couple of other options for training/testing IEMs

chan_resp = nan(size(X_all)); % fill this in with estimated channel responses

if strcmpi(which_data,'EEG')
    IEM_delay = [750 1250];   % time window to use for delay analyses (750:1250)
elseif strcmpi(which_data,'fMRI')
    IEM_delay = [5250 12000]; % time window to use for delay analyses
end
    
IEM_delay_tpts = tpts >= IEM_delay(1) & tpts < IEM_delay(2);

% NOTE: here, you could change which data we're using! (voxels, normalization scheme,
% etc)
delay_data = mean(d_all(:,:,IEM_delay_tpts),3);

% variable used for cross-validation:
% since we're doing leave-one-run-out (LORO), we'll use the run label to
% sort trials into training/testing sets
cv_all = r_all; % NOTE: try other schemes, like: cv_all = ceil(rand(size(r_all))*n_folds); (below)

%cv_all = ceil(rand(size(r_all))*50);  % try this one too! change the
%number to change the # of folds
cv_folds = length(unique(cv_all));

cvu = unique(cv_all);
for rr = 1:length(cvu)
    trnidx = cv_all~=cvu(rr) & ~excl_all; % train using all 'included' trials except testing run
    tstidx = cv_all==cvu(rr); % for now, reconstruct with all trials (can drop excluded trials later)
    
    trndata = delay_data(trnidx,:);
    tstdata = delay_data(tstidx,:);
    trnX    = X_all(trnidx,:);
    
    % estimate channel weights using 'training' data and predicted channel
    % responses on those trials
    % This is the heart of the encoding model here. Recall that
    % we are trying to fit fMRI voxel/EEG electrode responses as a linear
    % combination of channels, so the equation looks like:
    % C(n x k) x W(k x v) = R(n x v)
    % Where C is the n trials x k channels matrix of channel responses (trnx)
    % W is the k channels x v voxels weight matrix that we are trying to find
    % R is the n trials x v voxel response matrix (trndata)
    % To solve, use your favorite method for solving
    W_hat = pinv(trnX)*trndata;
    
    % use estimated channel weights to reconstruct channel responses using
    % 'test' data
    % Now to invert to predict channel responses given the weights
    % and the voxels, you go back through the same equation
    % C(n x k) = R(n x v) * pinv(W(k x v)
    chan_resp(tstidx,:) = tstdata*pinv(W_hat);
    
    clear w_hat trndata trnX tstdata trnidx tstidx;
end

% Plot 'raw' reconstructed channel responses

% now we have channel responses computed on each trial
% note that these are 'raw' - different positions will result in different
% profiles of reconstructed channel response functions

% so let's average all trials of matching position bins and plot them, just
% to see if reconstructions track remembered position approximately
figure;
hold on;  pbin_plot = pbins; pbin_plot(pbin_plot==0)=360;
for pp = 1:length(pu)
    plot(chan_centers,mean(chan_resp(c_all(:,2)==pu(pp)&~excl_all,:),1),'-','LineWidth',2,'Color',pos_colors(pp,:));
    plot([1 1]*pbin_plot(pu(pp)),[0 1],'--','LineWidth',1.5,'Color',pos_colors(pp,:));
end
xlim([22.5 382.5]);
set(gca,'XTick',0:90:360,'TickDir','out','FontSize',14);
xlabel('Channel center (\circ)');
ylabel('Reconstructed channel response (a.u.)');
title(sprintf('Delay-period reconstructions: %i to %i ms, %i CV folds',IEM_delay(1),IEM_delay(2),cv_folds));

% This looks pretty good to me! Try other cross-validation schemes if you
% like. If you switch to the 'randomized' mode, what happens when you split
% trials into fewer CV folds? More?



%% Convert channel responses to 'reconstructions'
%
% In the fundamentals tutorial, we stayed in 'channel space' - on each
% trial, we had one number for each channel. But in this dataset, stimulus
% positions were chosen from a uniform distribution around the screen. We
% were potentially hurting ourselves by binning trials based on their
% approximate position. (that is, trials at 23 deg polar angle were
% averaged with those at 67 deg polar angle). One way to get around this is
% to compute a 'reconstruction' using the estimated channel responses
% (chan_resp) and the channel sensitivity profile (or 'basis set'). By
% weighting each channel's sensitivity profile by the estimated channel
% response, we have a smooth function that we can shift more precisely.
% (you could even render those channel responses at a finer resolution with
% build_basis_polar_mat if you like!).

% weight channel profile by channel response on each trial to get a smooth
% reconstruction, then shift & align all reconstructions together

% first, let's weight each basis function by the corresponding channel
% response. in matlab, we can write this succinctly like:
recons_raw = chan_resp * myb_orig.';

% recons_raw is n_trials x length(angs) - each row is the reconstruction
% for a single trial.

% Now, we need to circularly shift weighted channel profile to align all trials

% make an empty variable to put things in
recons_aligned = nan(size(recons_raw));

% loop over trials
for tt = 1:size(c_all,1)
    % we want to adjust so that each position is set to 0
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned(tt,:) = circshift(recons_raw(tt,:),-1*shift_by);
end

% NOTE: the above assumes angs is a set of integers, spaced by 1. I'll
% leave it as an exercise to re-implement this procedure for an arbitrary
% feature space. there are hints about how to do this in Sprague et al,
% 2016 (and the code is available online)


% plot aligned reconstructions: all indiv trials and average
% (we're only going to plot non-excluded trials here)
figure;
subplot(3,1,[1 2]);
imagesc(angs,1:sum(~excl_all),recons_aligned(~excl_all,:)); axis ij;
ylabel('Trial (after exclusion)');
title('Aligned delay-period reconstructions');
xlim([-180 180]);
set(gca,'TickDir','out','XTick',-180:90:180,'Box','off','FontSize',14);

% mean across trials
subplot(3,1,3);
hold on;
plot(angs,mean(recons_aligned(~excl_all,:),1),'k-','LineWidth',2);
% also plot +- SEM over trials (not best metric, but gives some sense of
% errorbars)
plot(angs,mean(recons_aligned(~excl_all,:),1) + [-1;1]*std(recons_aligned(~excl_all,:),[],1)/sqrt(sum(~excl_all)),'k--');
xlim([-180 180]);
set(gca,'XTick',-180:90:180,'TickDir','out','FontSize',14);
xlabel('Angle (relative to target, \circ)');
ylabel('Reconstructed channel response (a.u.)');
title('Average over all trials');


%% how many trials are excluded from each bin? did this hurt us?

% NOTE: the below is only relevant for EEG analyses - we'll skip this for
% the brain camp tutorial, but I'll leave the code here for you to use

if strcmpi(which_data,'EEG')
    
    
    % now that we're pretty happy we have nice reconstructions during the delay
    % period, we can/should check whether this could possibly be due to having
    % more/less trials in a given bin on average than other bins
    
    figure;
    histogram(pbins(c_all(excl_all==0,2)),8);
    
    % to me, this looks basically ok - but let's implement a version of the
    % analysis that, within each position bin (on each CV fold) uses the same
    % number of trials per bin
    
    
    %% try IEM with fixed # of trials per bin
    %
    % this is almost identical to the leave-one-run-out CV above, except we're
    % going to limit the training set a bit (so lots of this is copy/paste)
    
    chan_resp_balanced = nan(size(X_all)); % fill this in with estimated channel responses
    
    % NOTE: for simplicity, I'm going to just use the first n trials of each
    % bin, where n is the number of trials in the bin with the fewest trials on
    % that CV fold. Ideally, you'd repeat this procedure 10ish times, drawing a
    % different set of trials each time. See Foster et al, 2016 paper & code
    % for details on this procedure
    %
    % we're using the same cross-validation index defined above
    
    % to approximate the Foster et al, 2016 analysis, you'd use 3 folds, and
    % then average all trials within a position bin before training. Advanced
    % students may wish to try and implement that here.
    
    for rr = 1:length(cvu)
        trnidx = cv_all~=cvu(rr) & ~excl_all; % train using all 'included' trials except testing run
        tstidx = cv_all==cvu(rr); % for now, reconstruct with all trials (can drop excluded trials later)
        
        % look in c_all(:,2) at trnidx for each pos bin, see which has min
        % count
        incl_trials_bin = nan(length(pu),1); % how many trials we include for each bin
        binned_trials = cell(length(pu),1); binned_X = cell(length(pu),1);
        for pp = 1:length(pu)
            incl_trials_bin(pp) = sum(c_all(trnidx,2)==pu(pp));
            binned_trials{pp} = delay_data(trnidx & c_all(:,2)==pu(pp),:);
            binned_X{pp} = X_all(trnidx & c_all(:,2)==pu(pp),:);
        end
        
        % within each bin, we'll only include this many
        incl_trials = min(incl_trials_bin); % note: can also do this w/ cellfun...
        binned_trials = cellfun(@(xx) xx(1:incl_trials,:), binned_trials,'UniformOutput',false);
        binned_X      = cellfun(@(xx) xx(1:incl_trials,:), binned_X,'UniformOutput',false);
        trndata = vertcat(binned_trials{:});
        trnX    = vertcat(binned_X{:});
        
        tstdata = delay_data(tstidx,:);
        
        % estimate channel weights
        W_hat = trnX \ trndata;
        
        % invert model to estimate channel responses
        chan_resp_balanced(tstidx,:) = tstdata/W_hat;
        
        clear w_hat trndata trnX tstdata trnidx tstidx;
    end
    
    % here's all the bookkeeping, etc, associated with re-aligning all trials
    % (again)
    % recenter, etc:
    % 1) circularly shift weighted channel profile to align all trials
    tmp_raw = chan_resp_balanced * myb_orig.';
    recons_aligned_balanced = nan(size(tmp_raw));
    
    % we want to adjust so that each position is set to 0
    for tt = 1:size(c_all,1)
        shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
        recons_aligned_balanced(tt,:) = circshift(tmp_raw(tt,:),-1*shift_by);
        clear shift_by;
    end
    clear tmp_raw;
    
    % now plot: we'll re-plot the original chan_resp's, and add the new
    % balanced ones as dashed lines
    
    figure; subplot(3,1,[1 2]);
    hold on;  pbin_plot = pbins; pbin_plot(pbin_plot==0)=360;
    for pp = 1:length(pu)
        plot(chan_centers,mean(chan_resp(c_all(:,2)==pu(pp)&~excl_all,:),1),'-','LineWidth',2,'Color',pos_colors(pp,:));
        plot(chan_centers,mean(chan_resp_balanced(c_all(:,2)==pu(pp)&~excl_all,:),1),':','LineWidth',2,'Color',pos_colors(pp,:));
        plot([1 1]*pbin_plot(pu(pp)),[0 1],'--','LineWidth',1.5,'Color',pos_colors(pp,:));
    end
    xlim([22.5 382.5]);
    set(gca,'XTick',0:90:360,'TickDir','out','FontSize',14);
    xlabel('Channel center (\circ)');
    ylabel('Reconstructed channel response (a.u.)');
    title(sprintf('Delay-period reconstructions: %i to %i ms, %i CV folds',IEM_delay(1),IEM_delay(2),cv_folds));
    
    % and the average across all trials, aligned, smoothed
    subplot(3,1,3); hold on;
    plot(angs,mean(recons_aligned(~excl_all,:),1),'k-','LineWidth',2);
    plot(angs,mean(recons_aligned_balanced(~excl_all,:),1),'k:','LineWidth',2);
    xlim([-180 180]);
    set(gca,'XTick',-180:90:180,'TickDir','out','FontSize',14);
    xlabel('Angle (relative to target, \circ)');
    ylabel('Reconstructed channel response (a.u.)');
    title('Average over all trials');
    legend({'All trials','Balanced'});
    
end


%% Quantifying reconstructions: 'fidelity' (Sprague et al, 2016)
%
% Now that we have reconstructions, what do we do with them? Visualizing
% them as channel response functions and/or spatial reconstructions (as
% plotted above) is certainly useful - but it can be even more useful to
% distill the graded response profile into one or more parameters
% describing it. Here, we'll go over a few different ways to quantify these
% reconstructions (focusing in this tutorial on 1-d reconstructions, but
% similar principles apply to 2-d reconstructions, see Sprague & Serences,
% 2013).
%
% We'll start simple: let's compute a metric that is big when the
% reconstruction is, on average, peaked in the 'correct' position, and
% hovers around zero when there's no information. There are two ways to do
% this, both with substantial precedent in the literature.
%
% First, let's compute what I call "FIDELITY" - this amounts to the
% 'strength' of the reconstruction in the 'direction' of the stimulus on
% that trial. Because, in this tutorial (and with color/orientation/polar
% angle stimuli), we're dealing with stimuli spanning a circular space, we
% can essentially just take the circular mean of the reconstruction. That
% is, if we plot the reconstruction in polar coordinates, the mean value of
% that plot can be considered a point-wise estimate of the reconstruction.
% The direction in which that mean vector points can be used to infer the
% represented feature value on that trial, and its length the 'strength'.
% But, of course the direction can be 'wrong' - so instead of reporting the
% 'strength' (vector length), we report the 'fidelity' - the vector length
% projected onto the 'correct' feature value.

% Let's start with some examples. We'll plot reconstructions as above, averaged across
% each position bin, in polar coordinates. (Note that, for simplicity, I'm
% going to use the first set of reconstructions we computed - not the ones
% computed with 'balanced' training sets. Feel free to swap out data
% variables below to see how things change!)

% (note: we're using binned averages, not perfectly-aligned reconstructions
% here, for simplicity & visualization. you'd actually want to use the
% perfectly-aligned data, as we see below)

figure; 
if ~old_matlab_ver,ax1 = subplot(1,2,1,polaraxes);else,ax1 = subplot(1,2,1);end
hold on;
for pp = 1:length(pu)
    thisidx = ~excl_all & c_all(:,2)==pu(pp);
    
    % plot the reconstruction in polar coords
    if ~old_matlab_ver
      polarplot(deg2rad(angs),mean(recons_raw(thisidx,:),1),'-','Color',pos_colors(pp,:),'LineWidth',1.5);
    else
      hPolar = polar(ax1,deg2rad(angs),mean(recons_raw(thisidx,:),1),'-');
      set(hPolar,'Color',pos_colors(pp,:),'LineWidth',1.5);
    end
    clear thisidx;
end
title('Aligned reconstructions (binned)');

% now let's focus on one bin (which_bin)
if ~old_matlab_ver,ax2 = subplot(1,2,2,polaraxes);else,ax2 = subplot(1,2,2);end
hold on; 
which_bin = 6;
title(sprintf('Example bin (%i)',which_bin));


% first, draw the same as above for 'background'

% get the trials we're using:
thisidx = ~excl_all & c_all(:,2)==which_bin;

if ~old_matlab_ver
  polarplot(deg2rad(angs),mean(recons_raw(thisidx,:),1),'-','Color',pos_colors(which_bin,:),'LineWidth',1.5);
else
  hPolar = polar(deg2rad(angs),mean(recons_raw(thisidx,:),1),'-');
  set(hPolar,'Color',pos_colors(which_bin,:),'LineWidth',1.5);
end
  
% and draw the mean of this position bin:
if ~old_matlab_ver
  polarplot([1 1]*deg2rad(pbins(which_bin)),[0 1],':','Color',pos_colors(which_bin,:),'LineWidth',1);
else
  hPolar = polar([1 1]*deg2rad(pbins(which_bin)),[0 1],':');
  set(hPolar,'Color',pos_colors(which_bin,:),'LineWidth',1);
end
  

% we want to do the vector average; let's add some example vectors
which_angs = (-135:45:180)-22; % (offsetting them from bin center)
% (note: these are JUST FOR VISUALIZATION) - we'll compute the real vector
% average below with all angles of the reconstruction. EXERCISE: is that
% really necessary? how few angles can we get away with? what are the
% constraints on which angles must be included?

for aa = 1:length(which_angs)
    % figure out which bin of recons_raw to average for each plotted angle
    angidx = find(angs==which_angs(aa)); % note: could replace w/ round, etc, for non-matching angles....
    if ~old_matlab_ver
      polarplot([1 1]*deg2rad(which_angs(aa)),[0 1]*mean(recons_raw(thisidx,angidx),1),'-','LineWidth',2,'Color',[0.3 0.3 0.3]);
    else
      hPolar = polar([1 1]*deg2rad(which_angs(aa)),[0 1]*mean(recons_raw(thisidx,angidx),1),'-');
      set(hPolar,'LineWidth',2,'Color',[0.3 0.3 0.3]);
    end
end

% now, we compute a vector sum of ALL points on the reconstruction
%
% but remember, this is a CIRCULAR variable - we can't just take the mean
% of recons_raw of these trials, so what do we do?

% first, convert average reconstruction to X,Y
[tmpx,tmpy] = pol2cart(deg2rad(angs),mean(recons_raw(thisidx,:),1));

% now, take the mean of X,Y
mux = mean(tmpx); muy = mean(tmpy); clear tmpx tmpy;

% and convert back to polar coordinates:
[muth,mur] = cart2pol(mux,muy); % NOTE: muth is in RADIANS here

% and draw:
if ~old_matlab_ver
  polarplot(muth*[1 1],mur*[0 1],'-','LineWidth',4,'Color',pos_colors(which_bin,:));
else
  hPolar = polar(muth*[1 1],mur*[0 1],'-');
  set(hPolar,'LineWidth',4,'Color',pos_colors(which_bin,:));
end

% Hopefully you can see that the average (thick colored line) aligns very
% nicely with the 'true' polar angle (thin dashed line) of this bin. So
% let's quantify this by projecting the vector mean (thick line) onto a
% unit vector in the direction of the bin center. Projection is just the
% dot product of the vector mean (mux,muy)<dot>(TrueX,TrueY), and TrueXY
% are just cosd(th) and sind(th), where theta is the bin center:
this_fidelity = dot([mux muy],  [cosd(pbins(which_bin)) sind(pbins(which_bin))]);
text(-pi/2,2,sprintf('Fidelity = %.03f',this_fidelity),'HorizontalAlignment','center');

% draw fidelity:
if ~old_matlab_ver
  polarplot([1 1]*deg2rad(pbins(which_bin)),this_fidelity*[0 1],'-','LineWidth',4,'Color','k');
else
  hPolar = polar([1 1]*deg2rad(pbins(which_bin)),this_fidelity*[0 1],'-');
  set(hPolar,'LineWidth',4,'Color','k');
end



% EXERCISE: does the order of operations here matter? can you compute
% fidelity on each trial then average, or do you need to compute fidelity
% on the average?

%% simplifying the computation of 'fidelity'
% So - we know how to project the vector average of the reconstruction onto
% the unit vector in the 'correct' direction, which lets us quantify
% fidelity. But, we've also computed 'aligned' reconstructions (indeed, we
% almost always do this in most analyses). In that case, the projection
% onto the 'correct' position amounts to projecting onto the x axis - so we
% can simplify our calculation of fidelity when using aligned data:
% F = mean( r(theta) * cos(theta) )  where r(theta) is 0-centered
% reconstruction value

% let's make a function that does this:
%compute_fidelity = @(rec) mean( rec .* cosd(angs) ,2);
compute_fidelity = @(rec) mean( bsxfun(@times,rec,cosd(angs)) ,2);

% compute fidelity for each trial
all_fidelity_delay = compute_fidelity(recons_aligned);


% plot distribution of fidelity across all trials
figure;
subplot(1,2,1);
histogram(all_fidelity_delay(~excl_all));
xlabel('Fidelity');
title('Fidelity across all trials');
set(gca,'TickDir','out');

% plot average fidelity (+std dev) for each position bin
subplot(1,2,2); hold on;
for pp = 1:length(pu)
    thisidx = ~excl_all & c_all(:,2)==pu(pp);
    thismu  = mean(all_fidelity_delay(thisidx));
    thiserr = std(all_fidelity_delay(thisidx))/sqrt(sum(thisidx));
    
    plot(pbins(pp)*[1 1],thismu + thiserr*[-1 1],'-','LineWidth',1.5,'Color',pos_colors(pp,:));
    plot(pbins(pp), thismu, 'wo','MarkerSize',8,'MarkerFaceColor',pos_colors(pp,:))
    clear thismu thiserr;
end
xlabel('Position bin center (\circ)');
ylabel('Mean fidelity');
title('Fidelity for each position bin');
tmpylim = get(gca,'YLim');
ylim([0 tmpylim(2)]);
xlim([-22.5 315+22.5]);
set(gca,'TickDir','out','XTick',0:90:360);

% what do you notice about the fidelity metric? Which positions are most
% poorly represented? Try using different position bins for the fidelity
% example computation above (change "which_bin")

%% Quantifying reconstructions: 'slope' (Foster et al, 2016)

% It's also common in the literature to just flip reconstructions (or, more
% commonly, channel response functions) over on themselves about the
% midpoint, average, and compute the slope of the channel resp vs channel
% center best-fit line.
%
% In Foster et al, 2016 (and many subsequent publications), they flip the
% channel response values (chan_resp here; n_chans long). But, let's see if
% we can get fancy and do this with the aligned reconstructions
% (recons_aligned) instead.
%
% We'll use only points in angs that have values on both sides of center
% (so -179:-1 and 1:179), and we'll flip the + side of angs (where the
% slope is negative) so that we overlay two up-going curves:

ang_stay = angs>=-179 & angs <= -1;
ang_flip = angs>=1 & angs <= 179;

% make a variable that's n_trials x length(ang_stay) x 2 and insert the
% aligned reconstruction into it
tmp_slope = nan(size(recons_aligned,1),sum(ang_stay),2);

tmp_slope(:,:,1) = recons_aligned(:,ang_stay);
tmp_slope(:,:,2) = fliplr(recons_aligned(:,ang_flip)); % flip these left/right

% look at the slope for all trials
figure;
hold on;
plot(angs(ang_stay),mean(tmp_slope(~excl_all,:,1),1),'-','LineWidth',1.5);
plot(angs(ang_stay),mean(tmp_slope(~excl_all,:,2),1),'-','LineWidth',1.5);

mu_slope = mean(tmp_slope,3); % average both flanks
plot(angs(ang_stay),mean(mu_slope(~excl_all,:),1),'k-','LineWidth',2);

% and compute the slope of best-fit line to the average
F = polyfit(angs(ang_stay),mean(mu_slope(~excl_all,:),1),1); % linear fit
plot(angs(ang_stay),polyval(F,angs(ang_stay)),'k--','LineWidth',2);

legend({'Left flank','Right flank (flipped)','Average','Best-fit line'},'Location','best');

set(gca,'FontSize',14,'Box','off','TickDir','out','XTick',-180:45:0);
xlabel('Relative position (\circ)');
ylabel('Reconstruction activation (a.u.)');
title('Computing slope (mean across trials)');
text(-160,0.8,sprintf('Slope = %0.03f/\\circ',F(1)));


% if you want slopes for each trial individually, need to loop over trials




%% Reconstructions through time

% How do neural representations of stimulus features evolve over the course
% of a trial? In these studies, we likely expect representations to be
% quite stable, since subjects are just maintaining a fixed point in space
% over an extended delay interval. But, some studies may involve asking
% participants to manipulate information in mind, or shift their attention,
% etc. See King & Dehaene 2014 for some neat examples for how to assay the
% temporal evolution of neural representations
%
% There are a bunch of ways we can compare representations over time. The
% most common, simplest, and most likely to work is to apply the IEM
% procedure at each timepoint during the trial in-turn. That is, to compute
% weights using a subset of trials at timepoint t, then reconstruct using
% activity at that same timepoint t. Let's give it a shot:

% which time range do we want to try reconstructing?
if strcmpi(which_data,'EEG')
    recon_tpt_range = [-250 1750]; % from just before delay to end
else
    recon_tpt_range = [-750 12000]; % from just before delay to end
end


% make a variable indexing into 3rd dim of data corresponding to
% interesting delay timepoints:
tptidx = find(tpts>=recon_tpt_range(1) & tpts<=recon_tpt_range(2));
% if we like, we can subsample:
ss_amt = 1; % how much we subsample - we take every n'th sample (1 = all samples, 2 = 1/2, 4 = 1/4, etc)
tptidx = tptidx(1:ss_amt:end);

% let's add an option to average neighboring timepoints/smooth (if we set
% this to the samplikng period, we won't be averaging at all)
avgwindow = 1000/sampling_rate; % average samples within this/2 on each side of tpt (ms)
%avgwindow = 50;

% make a variable to fill with reconstructions at each timepoint
chan_resp_t = nan(size(c_all,1),n_chan,length(tptidx)); % trial x channel x tpts

tpts_recon = tpts(tptidx); % when plotting reconstruction timecourses, use these

for tt = 1:length(tptidx)
    
    % select which tpts to average for this iteration:
    % all timepoints within avgwindow/2 of tpts(tptidx(tt))
    this_tpt_idx = abs(tpts-tpts(tptidx(tt))) <= avgwindow/2;
    
    % we can use this to index into our original data variable: d_all
    this_delay_data = mean(d_all(:,:,this_tpt_idx),3);
    
    for rr = 1:length(cvu)
        trnidx = r_all~=cvu(rr) & ~excl_all; % train using all 'included' trials except testing run
        tstidx = r_all==cvu(rr); % for now, reconstruct with all trials (can drop excluded trials later)
        
        trndata = this_delay_data(trnidx,:);
        tstdata = this_delay_data(tstidx,:);
        trnX    = X_all(trnidx,:);
        
        W_hat = trnX \ trndata; % estimate channel weights
        
        chan_resp_t(tstidx,:,tt) = tstdata/W_hat;
        
        clear w_hat trndata trnX tstdata trnidx tstidx;
    end
    clear this_tpt_idx this_delay_data;
end

% recenter, etc:
tmp_raw = nan(size(chan_resp_t,1),length(angs),size(chan_resp_t,3));
for tt = 1:size(chan_resp_t,3)
    tmp_raw(:,:,tt) = chan_resp_t(:,:,tt) * myb_orig.';
end
recons_aligned_t = nan(size(tmp_raw));

% we want to adjust so that each position is set to 0
for tt = 1:size(c_all,1)
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned_t(tt,:,:) = circshift(tmp_raw(tt,:,:),-1*shift_by);
    clear shift_by;
end
clear tmp_raw;

% plot the reconstructions as an image
figure;
subplot(2,1,1); hold on;
imagesc(tpts_recon,angs,squeeze(mean(recons_aligned_t(~excl_all,:,:),1)));
%axis tight;
ylabel('Position (\circ)');
title(sprintf('Matched training; subsample: %ix, window: %i ms',ss_amt,avgwindow));
xlim(recon_tpt_range);ylim([angs(1) angs(end)]);
set(gca,'TickDir','out');

% and quantify: compute fidelity over time
fidelity_t = squeeze(compute_fidelity(recons_aligned_t)); % trials x tpts

subplot(2,1,2);
plot(tpts_recon,mean(fidelity_t(~excl_all,:),1),'k-','LineWidth',2);
xlabel('Time (ms)');
ylabel('Fidelity');
title('Fidelity through time (matched training)');
xlim(recon_tpt_range);
set(gca,'TickDir','out');

%% 3d surface plot
% and just because it's fun, a surface version
% (I adapted this directly from Josh's scripts)
figure;
surf(tpts_recon,angs,squeeze(mean(recons_aligned_t(~excl_all,:,:),1)),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp
h=findobj('type','patch');
set(h,'linewidth',2)
hold on
set(gca, 'box','off')
set(gca,'color','none')
set(gca,'LineWidth',2,'TickDir','out');
set(gca,'FontSize',14)
view(3)
%axis([x(1) x(end) tpts_recon(1) em.time(end) lim]);
set(gca,'YTick',[-180:90:180])
%set(gca,'XTick',[recon_tpt_range(1):500:recon_tpt_range(2)])
title(which_data,'FontSize',16)
ylabel('Channel Offset (\circ)');
xlabel('Time (ms)');
set(get(gca,'xlabel'),'FontSize',14,'FontWeight','bold')
set(get(gca,'ylabel'),'FontSize',14,'FontWeight','bold')
zlabel({'Channel'; 'response (a.u.)'}); set(get(gca,'zlabel'),'FontSize',14,'FontWeight','bold')
%set(get(gca,'ylabel'),'rotation',90); %where angle is in degrees
grid off



%% Fixed encoding model across time
%
% Above, we recomputed the encoding model at each point in time, and this
% gave us nice, stable reconstructions through the entire delay. But how do
% we interpret that? If reconstructions are 'stronger' at one point in time
% compared to another, what does that mean?
%
% It's hard to tell in this case. The reason is, the encoding model itself
% changes at each point in time, and so the space to which the IEM projects
% the electrode-space data is different! This tells us that there is
% information at each point in time, but does not say anything about the
% stability of this information - we cannot directly compare
% reconstructions at the beginning and end, except to say that there's
% more/less accessible information, not necessarily that the representation
% itself is waxing/waning. This is a subtle, but important point, and it
% applies to ALL comparisons of IEM results across conditions. If the model
% used to reconstruct across a comparison differs, we don't conclusively
% know whether the change in the reconstruction is due to a change in model
% fit, a change in the neural representation, or a combination of the two.
% This is a similar problem to seeing someone report a change in a ratio
% metric, without any report about the numerator or denominator
% individually - it is imposisble to know whether the numerator,
% denominator, or both changed across the manipulation of interest.
%
% One way around this is to use a 'fixed' encoding model - one estimated at
% a particular point in time, or using a separate set of data entirely
% (perhaps you have a localizer task, or a single-item set of trials, etc).
% Here, because there is only one 'condition', we'll estimate a fixed
% encoding model using the delay timepoints we were analyzing above. Then,
% we'll reconstruct using activity at each timepoint in turn. This will let
% us directly compare early- vs late-delay reconstructions, since they were
% computed using the same encoding model.


% we'll use the same timepoints as above for reconstruction, but only a
% single epoch for training:
%trn_window = delay_window;
trn_window = IEM_delay;
trn_tpts = tpts >= trn_window(1) & tpts <= trn_window(2); delay_tpts; % from above - but experiment w/ different values!

chan_resp_t_fix = nan(size(c_all,1),n_chan,length(tptidx));


for rr = 1:length(cvu)
    
    % on every cross-validation fold, train with a single set of exemplars
    % across trials averaged over trn_tpts
    
    trnidx = r_all~=cvu(rr) & ~excl_all; % train using all 'included' trials except testing run
    tstidx = r_all==cvu(rr); % for now, reconstruct with all trials (can drop excluded trials later)
    
    trndata = mean(d_all(trnidx,:,trn_tpts),3);
    trnX    = X_all(trnidx,:);
    
    % estimate IEM once per CV fold
    W_hat = trnX \ trndata; % estimate channel weights
    
    for tt = 1:length(tptidx)
        
        % select which tpts to average for this iteration:
        % all timepoints within avgwindow/2 of tpts(tptidx(tt))
        this_tpt_idx = abs(tpts-tpts(tptidx(tt))) <= avgwindow/2;
        
        tstdata = mean(d_all(tstidx,:,this_tpt_idx),3);
        
        % use that estimated IEM (same across time!) to reconstruct each
        % timepoint
        chan_resp_t_fix(tstidx,:,tt) = tstdata/W_hat;
        
        clear this_tpt_idx tstdata;
    end
    clear w_hat trndata trnX trnidx tstidx;
end

% recenter, etc:
tmp_raw = nan(size(chan_resp_t_fix,1),length(angs),size(chan_resp_t_fix,3));
for tt = 1:size(chan_resp_t_fix,3)
    tmp_raw(:,:,tt) = chan_resp_t_fix(:,:,tt) * myb_orig.';
end
recons_aligned_t_fix = nan(size(tmp_raw));

% we want to adjust so that each position is set to 0
for tt = 1:size(c_all,1)
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned_t_fix(tt,:,:) = circshift(tmp_raw(tt,:,:),-1*shift_by);
    clear shift_by;
end
clear tmp_raw;


% plot the reconstructions from fixed IEM as an image
figure;
subplot(2,1,1); hold on;
imagesc(tpts_recon,angs,squeeze(mean(recons_aligned_t_fix(~excl_all,:,:),1)));
%axis tight;
ylabel('Position (\circ)');
title(sprintf('Fixed training (%i to %i ms)',trn_window(1),trn_window(2)));
xlim(recon_tpt_range);ylim([angs(1) angs(end)]);
set(gca,'TickDir','out');

% and quantify: compute fidelity over time
fidelity_t_fix = squeeze(compute_fidelity(recons_aligned_t_fix)); % trials x tpts

subplot(2,1,2); hold on;
plot(tpts_recon,mean(fidelity_t_fix(~excl_all,:),1),'k-','LineWidth',2);
plot(trn_window,[1 1]*-0.05,'k-','LineWidth',5);
xlabel('Time (ms)');
ylabel('Fidelity');
title('Fidelity through time (fixed training)');
xlim(recon_tpt_range);
set(gca,'TickDir','out');

%% Comparing reconstructions from 'dynamic' vs 'fixed' encoding models
figure;
hold on;
plot(tpts_recon,mean(fidelity_t(~excl_all,:),1),'-','LineWidth',2);
plot(tpts_recon,mean(fidelity_t_fix(~excl_all,:),1),'-','LineWidth',2);
plot(trn_window,[1 1]*-0.05,'k-','LineWidth',5); % training window
legend({'Matched training','Fixed training','Training window'});

ylabel('Fidelity');
xlabel('Time (ms)');
title('Comparing fidelity across IEM training regimes');


%% are these results by chance? shuffle trial labels at training
%
% to compute 'quality' of reconstructions on each iteration, we'll just use
% fidelity metric above (though slope, etc, also usable).
%
% How do we know the results we've shown above are 'real'? Especially for a
% single subject! A common way of doing this is removing structure from the
% 'training' dataset (typically by shuffling trial labels), computing an
% encoding model, and using that 'null' encoding model to reconstruct each
% trial of the test data (with their labels intact). Wash, rinse, repeat
% ~100-1000x. You can use this surrogate dataset to compute any relevant
% stats, like fidelity/slope, and compare the 'real' value you see (at each
% timepoint even!) to this surrogate distribution to get a sense for the
% likelihood you'd see a real result by chance.
%
% For simplicity, we'll just demonstrate this on the 'delay period' version
% of the analysis. But, it's quite easy to extend this approach to any of
% the temporal analyses we discussed (just computationally-intensive).

niter = 100; % how many times do we shuffle?

% seed the random number generator
rng(2783829); % so we get the same result each time!

% initialize a few variables for 'null' reconstruction distributions
chan_resp_shuf = nan(size(c_all,1),n_chan,niter); % trials x channels x shuffling iterations
tmp_raw = nan(size(chan_resp_shuf,1),length(angs),niter);
for ii = 1:niter
    
    for rr = 1:length(cvu)
        trnidx = r_all~=cvu(rr) & ~excl_all; % train using all 'included' trials except testing run
        tstidx = r_all==cvu(rr); % for now, reconstruct with all trials (can drop excluded trials later)
        
        trndata = delay_data(trnidx,:);
        tstdata = delay_data(tstidx,:);
        trnX    = X_all(trnidx,:);
        
        % shuffle rows of trnX before training
        trnX = trnX(randperm(size(trnX,1)),:);
        
        % now, there's no reliable correspondence between the predicted
        % response of each channel on a trial and the observed data on that
        % trial - any accidental structure that was there should be
        % removed!
        
        W_hat = trnX \ trndata; % estimate channel weights
        
        chan_resp_shuf(tstidx,:,ii) = tstdata/W_hat;
        
        clear w_hat trndata trnX tstdata trnidx tstidx;
    end
    
    tmp_raw(:,:,ii) = chan_resp_shuf(:,:,ii) * myb_orig.';
    
end

% now compute reconstructions, align
recons_aligned_shuf = nan(size(tmp_raw));

% we want to adjust so that each position is set to 0
for tt = 1:size(c_all,1)
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned_shuf(tt,:,:) = circshift(tmp_raw(tt,:,:),-1*shift_by);
end
clear tmp_raw;

% and compute their fidelity
fidelity_shuf = squeeze(compute_fidelity(recons_aligned_shuf));

null_dist = mean(fidelity_shuf(~excl_all,:),1); % 100 'null' mean fidelity values
p = mean( null_dist >= mean(all_fidelity_delay(~excl_all)) );
figure;histogram(null_dist,10); hold on;
plot([1 1] * mean(all_fidelity_delay(~excl_all)), [0 15],'r--','LineWidth',2);
xlabel('Fidelity');
title(sprintf('Comparing true fidelity to shuffled (p = %0.03f)',p));




%% Using channel-based models to 'decode' stimulus features
%
% One particularly useful application of the IEM technique is visualizing
% model-based reconstructed channel response profiles that show how
% stimulus features might be represented, given a particular set of model
% assumptions. (ASIDE: see Gardner & Liu, 2019 for an example of how
% changing model assumptions will necessarily change results. But, note
% that knowing properties of the model allows one to interpret channel
% response profiles, even if they've been changed). 
%
% Above, we talked about quantifying the 'strength' of stimulus
% representations within channel response profiles using a 'fidelity' and a
% 'slope' metric. In addition to knowing how strong the representation is,
% we may also be interested in learning now 'accurate' it is - what feature
% value is represented by the channel response profile? And how sure are we
% that's the accurate feature value?
%
% We'll approach this problem in two stages. First, decoding feature
% values.
%
% There are *lots* of ways to do this. For simplicity, let's apply a
% population vector approach. A high-level assumption of the IEM method is
% that we're projecting activation patterns from signal space
% (electrodes/voxels) into a model-based 'information space', with the
% slightly-accurate interpretation that each modeled channel is roughly
% analogous to a neural population (again, see Gardner & Liu, 2019). Thus,
% we'll pretend we can treat the reconstructed channel responses as
% pseudo-populations, and so we'll apply some relatively straightforward
% population decoding metrics (see:
% http://www.math.pitt.edu/~bdoiron/assets/dayan_abbott.pdf)
%
% Let's start with the linear version - we're going to let each 'channel'
% vote for its preferred feature value based on how strongly it responds
% during each test trial. 
%
% So - the decoded feature value is the angle of the vector sum of a unit vector
% pointing in the direction of each channel's preferred value weighted by
% its response


% chan_resp is our variable corresponding to the 'raw', unrotated/aligned
% channel response profiles - one number for each trial, for each modeled
% channel. chan_centers is our variable listing the tuned position of each
% channel, in degrees. 

% let's decode an example trial first:

decode_trial = 18; % 18 is good...
figure;subplot(1,2,1);hold on;

% plot the channel response profile for this trial
plot(chan_centers,chan_resp(decode_trial,:),'-','LineWidth',2);
plot([1 1]*mod(c_all(decode_trial,1),360),[0 1],'r--','LineWidth',1.5);
xlabel('Modeled channel center (\circ)');
ylabel('Channel response profile');
title(sprintf('Trial %i',decode_trial));
set(gca,'XTick',chan_centers(2:2:end),'TickDir','out','FontSize',16);

% plot the vector average visually (in polar coords)

% first, the channel responses in polar coords
subplot(1,2,2);
polarplot(deg2rad(chan_centers),chan_resp(decode_trial,:),'k:','LineWidth',1.0);
hold on;

polarplot([1 1]*deg2rad(c_all(decode_trial,1)),[0 1],'r--','LineWidth',.75);

for cc = 1:length(chan_centers)
    polarplot([1 1]*deg2rad(chan_centers(cc)),[0 chan_resp(decode_trial,cc)],'-','LineWidth',2);
end

% to compute the poulation vector, we just compute the sum of each
% population's preferred feature value (unit vector in each channel's tuned
% direction) weighted by its reconstructed response:

% this computes each channel's unit vector (given by [cosd(th) sind(th)]
% and multiplies it by the corresponding channel's response, and sums
% across channels
vec_sum_ex = mean(chan_resp(decode_trial,:).*[cosd(chan_centers);sind(chan_centers)],2);

% NOTE: above I'm computing/plotting the mean, so it's easier to see on the
% polar plot - the angle is of course the same for mean and sum

% vec_sum_ex is the (x) and (y) coords of the vector sum - we'll convert
% them to polar and add to the plot
[vec_sum_ex_pol(1),vec_sum_ex_pol(2)] = cart2pol(vec_sum_ex(1),vec_sum_ex(2));
polarplot([1 1]*vec_sum_ex_pol(1),[0 vec_sum_ex_pol(2)],'k-','LineWidth',2.5);

%title(sprintf('Actual: %i \circ ; Decoded: %i \circ',mod(c_all(decode_trial,1),360),round(mod(rad2deg(vec_sum_ex_pol(1)),360))));

fprintf('Trial %i\nActual position:\t%f deg\n',decode_trial,mod(c_all(decode_trial,1),360));
fprintf('Decoded position:\t%f deg\n',mod(rad2deg(vec_sum_ex_pol(1)),360));

% add the decoded value to the original channel response profile plot
subplot(1,2,1); hold on;
plot(mod(rad2deg(vec_sum_ex_pol(1)),360)*[1 1],[0 1],'k--','LineWidth',1.5);
legend({'Data','Actual','Decoded'},'Location','best');


% EXERCISE: look at other trials - not every one is perfect, but often they
% look quite good! (remember, this is *single-trial* decoding of a
% *continuous feature value*



%% Now, let's decode ALL trials!

% Turns out this is pretty easy overall - just a for loop (exercise for the
% reader), or some matrix math

vec_sum_all = chan_resp * [cosd(chan_centers.') sind(chan_centers.')];

decoded_angs = atan2d(vec_sum_all(:,2),vec_sum_all(:,1));

% and plot them - compare each trial's decoded feature value to its actual value 
figure; subplot(1,2,1); hold on;
try
    scatter(mod(c_all(:,1),360),mod(decoded_angs,360),15,'k','filled','MarkerFaceAlpha',0.5);
catch
    scatter(mod(c_all(:,1),360),mod(decoded_angs,360),15,'k','filled'); % in case matlab version can't handle transparent points
end
%plot unity line
plot([0 360],[0 360],'k--');

xlabel('Actual position (\circ)');
ylabel('Decoded position (\circ)');
title('Channel-based vector average');

xlim([0 360]);ylim([0 360]);

axis equal square;

set(gca,'TickDir','out','FontSize',14,'XTick',0:90:360,'YTick',0:90:360);

% and let's look at a histogram of decoding errors
tmperr = decoded_angs - c_all(:,1);
decoding_err = mod((tmperr+180),360)-180;
clear tmperr;

subplot(1,2,2);hold on;
histogram(decoding_err); xlim([-180 180]);
xlabel('Decoding error (\circ)'); 
title('Single-trial decoding error');
set(gca,'TickDir','out','FontSize',14,'XTick',-180:90:180,'YTick',[]);

% Overall, it looks like the decoder is doing a really good job! Try
% changing aspects of the analysis above - such as cross-validation scheme,
% shape of basis functions, number of basis functions, etc. In general,
% especially for 'large-scale' features like spatial position, single-trial
% fMRI decoding is very effective and accurate. (tip: to run things more
% quickly, set the niter variable above to "0", or comment out that section
% entirely - then you can change params and run the whole script at once)
%
% QUESTION: what does this method achieve that a categorical method, like
% SVM/LDA, cannot? 
%
% ANSWER: on each trial, we decode a continuous feature value (here, a
% point estimate) - not which category a trial is a member of. As an
% exercise, try also using a categorical classifier, like the old
% classify.m in matlab, etc. 


%% Decoding "uncertainty" - how sure are we about what we're decoding?
%
% The linear IEM technique described throughout the tutorial is best
% thought
% of as a model-based dimensionality reduction method. We estimate a
% mapping between measured signal space, with one dimension for each
% measurement (voxel, electrode, etc), and a model-based 'information
% space', with one dimension for each modeled feature-selective channel.
% Because, so far, we're doing simple linear operations, there are
% straightforward ways to compute this mapping, and the act of
% reconstructing modeled channel responses for novel measured neural
% activity patterns is really just projecting from one space to another.
% Importantly, we're projecting a *point estimate* to recover the
% activation of each modeled channel.
%
% See Sprague, Boynton & Serences, 2019, bioRxiv for a discussion of how
% changes in model assumptions (e.g., channel shapes) can predictably
% change properties of these model-based channel resposne profiles.
%
% However, in some cases, we're interested in not just which feature value
% is represented as read out from this point estimate using a linear
% decoder (the population vector method described above), but the *range* of
% feature values with which our measurement is consistent. 
%
% That is - what is the likelihood over stimulus features, given the
% measured activation pattern and previously-estimated encoding model (W)?
%
% To estimate a likelihood function, we need to know something about the
% *noise* associated with our measurements. For the tutorial today, we're
% only going to consider independent noise for each measurement dimension
% (electrode/voxel). In real applications, you also need to consider
% correlated noise across measurement dimensions - that is, their
% covariance. 
%
% (NOTE: this is all adapted from a very simplified version of the method
% presented by van Bergen et al, 2015; van Bergen & Jehee, 2017; Liu et al,
% 2018)
%
% Let's start looking into noise for each dimension - we'll implement
% model-fitting stage of the IEM as described above, but let's focus on the
% *residuals* - the parts of the signal the best-fit encoding model does
% not explain. 

% to demonstrate, we'll estimate a model with all trials except run 1,
% compute residuals, use those residuals to estimate a noise model, and
% compute stimulus likelihoods using the trials in run 1.

decode_run = 1;

trnidx = r_all~=decode_run;
tstidx = r_all==decode_run;

W_hat = X_all(trnidx,:) \ delay_data(trnidx,:);

% W_hat is n_chans x n_vox

resid = delay_data(trnidx,:) - X_all(trnidx,:)*W_hat;

resid_dim = 10; % examine residuals for one signal dimension

figure;hold on;
histogram(resid(:,resid_dim),'Normalization','pdf');
xlabel('Residual');
title(sprintf('Signal dimension %i',resid_dim));
set(gca,'TickDir','out','FontSize',14);

% put a normal distribution on top...
myx = linspace(min(resid(:,resid_dim)),max(resid(:,resid_dim)),101);
plot(myx,normpdf(myx,mean(resid(:,resid_dim)),std(resid(:,resid_dim))),'r--');


% It looks like our residuals are normally distributed within a single
% voxel - cool! (EXERCISE: look at scatterplots of residuals across
% several voxels - e.g., using plotmatrix; you'll see *lots* of correlated
% variability, which we won't address here - see van Bergen et al, 2015)


% ok - so now we want to compute a likelihood function P(r|s) over s, where
% s is the position in polar angle and r is the vector of observed
% responses across signal dimensions (voxels). In computing the likelihood,
% we'll use our previously-estimated encoding model, W_hat, and a summary
% statistic of the residuals (their variance). 

% because we're assuming independent normally-distributed noise for each voxel, for a single stimulus value, s_i, and voxel, r_j, the likelihood
% P(r_j|s_i) is the normal distribution evaluated at the residual computed
% for that voxel (r_j) for that stimulus value (s_i), using the
% previously-estimated W_hat and noise model (std(s_i))


this_vox = resid_dim; this_stim_value = 45; % deg
% compute the voxel's predicted response given W_hat, basis
tmp_C = build_basis_polar_mat(this_stim_value,chan_centers,chan_width); % use same values as before...

% tmp_C is the predicted response of each channel for this possible
% stimulus value. we use that, along with W_hat and our encoding model B =
% C * W_hat, to compute how this voxel should respond to stimulus s_i:
% tmp_C*W_hat(:,this_vox)
%
% Because we know there's some noise (our independent voxel noise model),
% we can compute the likelihood of the observed state on a given trial
% (decode_trial) supposing stimulus s_i was presented - because we're
% assuming normally-distributed noise, this is just evaluating the normal
% distribution at the difference between the observed response and the
% predicted response:
tmp_lik = normpdf(delay_data(decode_trial,this_vox)-tmp_C*W_hat(:,this_vox),0,std(resid(:,this_vox)));

% equivalently written as:
%myp = delay_data(decode_trial,this_vox)-tmp_C*W_hat(:,this_vox); % predicted response
%myn = std(resid(:,this_vox));  % noise (std = sigma)
%tmp_lik2 = (1/(sqrt(2*pi)*myn)) * exp(-((myp).^2)/(2*myn.^2));

% now - we want to compute the likelihood of the response for each voxel,
% and the overall likelihood of response pattern r (vector) is the product
% of the likelihoods of s_i for each voxel r_i. Because we're computing the
% product of a probability, it can sometimes be easier to work in logs -
% instead of computing the likelihood, we compute the log likelihood. By
% taking the log on each side, the product becomes a sum, which is more
% numerically tractable.

tmp_LL_allvox = sum(log(normpdf(delay_data(decode_trial,:)-tmp_C*W_hat,0,std(resid))));
% or, equivalently written as a simplified form based on log of Gaussian:
% log(Gaussian) = log(1/(sqrt(2*pi)*myn)) + -1*(myp.^2)./(2*myn.^2);
%allp = delay_data(decode_trial,:)-tmp_C*W_hat; % predicted response
%alln = std(resid);
%tmp_LL_allvox2 = sum(  log(1./(sqrt(2*pi).*alln)) + -1*(allp.^2)./(2*alln.^2)  );

% if you inspect the values of tmp_lik_allvox, you'll see that the log
% likelihood is very negative - of course! we're summing logs of
% probabilities, which are negative. so, often, we use the negative log
% likelihood (I'll use that below). When using negative log likelihood
% (NLL), we could adopt minimization strategies to minimize NLL (rather
% than maximiizing LL) - this is useful in some contexts, but will not be
% discussed further here. just keep the -1's straight!


% finally, we want to compare the likelihood for *all* stimulus values -
% we'll just do 2:2:360 to keep things simple (must be 'round')

stim_vals = 2:2:360;
tmp_NLL_s = nan(length(stim_vals),1);

% for each stimulus value, compute its negative log likelihood
for ss = 1:length(stim_vals)
    tmp_C = build_basis_polar_mat(stim_vals(ss),chan_centers,chan_width); 
    tmp_NLL_s(ss) = sum(-1*log(normpdf(delay_data(decode_trial,:)-tmp_C*W_hat,0,std(resid))));
end

tmp_lik_s = exp(-tmp_NLL_s + min(tmp_NLL_s));

% and normalize to a probability distribution (so integral == 1)
tmp_lik_s = tmp_lik_s/(trapz(deg2rad(stim_vals),tmp_lik_s));

% look at the decoded likelihood function
figure; hold on;
%plot(stim_vals,exp(-tmp_NLL_s + min(tmp_NLL_s)),'k-'); % normalize to max 1
plot(stim_vals,tmp_lik_s,'k-'); % normalize to max 1
% above - we computed the negative log likelihood - probabilities are < 1,
% so negative log gives positive numbers - here, by removing the minimum of
% the NLL, and plotting the likelihood (exp(-1*NLL)), we normalize so that
% the 'maximum' value of the likelihood is 1 (exp(0) = 1).

plot(mod(c_all(decode_trial,1),360)*[1 1],[0 max(tmp_lik_s)],'k--');

xlabel('Position (\circ)');
ylabel('Likelihood (p(r|s))');
title(sprintf('Trial %i',decode_trial));
set(gca,'TickDir','out','FontSize',14,'XTick',0:90:360);
xlim([0 360]);

% one way to think about what's happening is we're 'trying out' different
% stimulus values, and comparing our actual observed response pattern to
% the noisy distribution of response patterns that would be consistent with
% each stimulus value we're trying out. Given a noise model and estimated
% encoding model, we can predict the *distribution* in voxel space of
% likely (and unlikely) response patterns for each stimulus - by trying
% different stimuli and comparing observed to predicted response patterns,
% we can find the set of stimuli that are most consistent with our observed
% response patterns. Because we can do this across all stimulus values, we
% can get a measure not just of the "most likely" stimulus value (that is -
% the decoded stimulus position) - but of the *uncertainty* with which the
% brain represents that stimulus position. Operationally, we'll define that
% as the circular standard deviation of the likelihood function. Previous
% work - especially that led by van Bergen et al (2015) - has shown a few
% ways in which this 'neural uncertainty' is related to behavior. 


% so - let's decode the stimulus feature value (approximate circular integral), and its
% uncertainty (circular std dev)

tmp_vec = sum(deg2rad(stim_vals(2)-stim_vals(1))*tmp_lik_s.*[cosd(stim_vals);sind(stim_vals)].');
tmp_pos = mod(atan2d(tmp_vec(2),tmp_vec(1)),360);

% for uncertainty, we want the circular standard deviation - or at least
% something directly proportional to it. I took this computation from van
% Bergen et al, 2015's code
tmp_unc = 180*sqrt(-2*log(sqrt(sum(tmp_vec.^2))))/pi;

% add the decoded feature value and uncertainty
hold on;
plot(tmp_pos,max(tmp_lik_s)+1,'kv','MarkerFaceColor','k','MarkerSize',5);
if (tmp_pos - tmp_unc) >= 0 && (tmp_pos + tmp_unc) <= 360
    plot(tmp_pos + [-1 1]*tmp_unc,[1 1] * max(tmp_lik_s)+1,'k-','LineWidth',1.5);
else
    % to be able to 'wrap around' the unc bar
    if tmp_pos - tmp_unc < 0
        % draw the left side
        plot([0 tmp_pos+tmp_unc], [1 1] * max(tmp_lik_s)+1,'k-','LineWidth',1.5);
        % and the right
        plot([mod(tmp_pos-tmp_unc,360) 360],[1 1] * max(tmp_lik_s)+1,'k-','LineWidth',1.5);
    elseif tmp_pos + tmp_unc > 360
        % draw the left side
        plot([0 mod(tmp_pos+tmp_unc,360)],[1 1] * max(tmp_lik_s)+1,'k-','LineWidth',1.5);
        % and the right
        plot([tmp_pos-tmp_unc 360],[1 1] * max(tmp_lik_s)+1,'k-','LineWidth',1.5);
    end
end


%% Decoding all trials

% Great! We've decoded a single trial's likelihood function, and from
% there, figured out the represented feature value (circular mean) and its
% uncertainty (circular standard deviation). Let's do this for every trial
% in the experiment, cross-validating by runs as above

% make a variable to save into
lf_all = nan(size(delay_data,1),length(stim_vals));


cv_all = r_all;
cvu = unique(cv_all);

% to make things more efficient, we'll compute the predicted channel
% responses for all feature values in stim_vals here - so we can do some
% matrix math below
c_pred_all = build_basis_polar_mat(stim_vals,chan_centers,chan_width); % n_trials x n_chans

% loop over cross-validation folds
for cv_idx = 1:length(cvu)
    
    trnidx = cv_all~=cvu(cv_idx);
    tstidx = cv_all==cvu(cv_idx);
    
    W_hat = X_all(trnidx,:) \ delay_data(trnidx,:);
    
    % W_hat is n_chans x n_vox
    resid = delay_data(trnidx,:) - X_all(trnidx,:)*W_hat;
    
    this_noise_model = std(resid); % our noise model is just the standard deviations of the residuals for each voxel
    
    this_delay_data = delay_data(tstidx,:);
    this_trial_idx = find(tstidx==1); % for inserting into lf_all, pos_all, unc_all
    
    % decode likelihood for each trial
    for tt = 1:size(this_delay_data,1)
        
        % for this trial (tt of this_delay_data), get predicted mean
        % response for each possible stim_val
    
        all_pred = this_delay_data(tt,:)-c_pred_all*W_hat; % predicted response
       
        % sum likelihoods over voxels
        tmp_NLL_s = -1 * sum(  log(1./(sqrt(2*pi).*this_noise_model)) + -1*(all_pred.^2)./(2*this_noise_model.^2) , 2 );
        
        % convert to a likelihood
        tmp_lik_s = exp(-tmp_NLL_s + min(tmp_NLL_s));
        
        % now store it as a pdf
        
        % need to compute integral - which means we need to 'wrap it
        % around' since we're doing approxiamte trapezoidal integration
        lik_for_int = [tmp_lik_s(end); tmp_lik_s];
        
        lf_all(this_trial_idx(tt),:) = tmp_lik_s/(trapz(deg2rad([0 stim_vals]),lik_for_int));
        
        clear tmp_NLL_s tmp_lik_s all_pred lik_for_int;
    end
     
end

% decode position, uncertainty
tmp_vec = deg2rad(stim_vals(2)-stim_vals(1))*lf_all*[cosd(stim_vals);sind(stim_vals)].';

pos_all = mod(atan2d(tmp_vec(:,2),tmp_vec(:,1)),360);

unc_all = 180*sqrt(-2*log(sqrt(sum(tmp_vec.^2,2))))/pi;


%% plot them!
%
% we'll plot likelihoods from 5 trials based on their uncertainty
% percentile (lowest, 25th, 50th, 75th, and highest percentile)

figure;
subplot(2,2,1); hold on;
prct_to_plot = [0 25 50 75 100];
unc_prctile = prctile(unc_all,prct_to_plot);

% look up trials at each percentile
trials_to_plot = nan(size(prct_to_plot));
for pp = 1:length(prct_to_plot)
    [~,trials_to_plot(pp)] = min(abs(unc_all-unc_prctile(pp)));
end


mycolors = lines(2*length(trials_to_plot));mycolors = mycolors(1:length(trials_to_plot),:);
for tt = 1:length(trials_to_plot)
    plot(stim_vals,lf_all(trials_to_plot(tt),:),'-','LineWidth',1.5,'Color',mycolors(tt,:));
    plot([1 1]*mod(c_all(trials_to_plot(tt),1),360),[0 max(lf_all(trials_to_plot(tt),:))],'--','Color',mycolors(tt,:));
    plot(pos_all(trials_to_plot(tt)),1.5+max(lf_all(trials_to_plot(tt),:)),'v','Color','k','MarkerFaceColor',mycolors(tt,:),'MarkerSize',5);
end
set(gca,'TickDir','out','FontSize',12,'XTick',0:90:360);
xlim([0 360]);
xlabel('Position (\circ)');
ylabel('Likelihood (p(r|s))');


% and a bar graph comparing their relative uncertainties
subplot(2,2,3);
hold on;
for tt = 1:length(trials_to_plot)
    bar(tt,unc_all(trials_to_plot(tt)),'FaceColor',mycolors(tt,:));
end

set(gca,'TickDir','out','XTick',1:length(trials_to_plot),'XTickLabel',prct_to_plot','FontSize',12);
xlim([0 length(trials_to_plot)+1]);
ylabel('Uncertainty (\circ)');
xlabel('Percentile among all trials');


% next - a scatterplot of decoded versus actual position, color-coded by
% uncertainty
subplot(2,2,2); hold on;
try % account for matlab's that can't do translucent datapoints
    scatter(mod(c_all(:,1),360),pos_all,25,unc_all,'filled','MarkerFaceAlpha',0.5);
catch
    scatter(mod(c_all(:,1),360),pos_all,25,unc_all,'filled');
end
plot([0 360],[0 360],'k--');
colormap parula;
xlim([0 360]);ylim([0 360]);
axis equal square;
set(gca,'TickDir','out','FontSize',12,'XTick',0:180:360,'YTick',0:180:360);
xlabel('Actual position (\circ)');
ylabel('Decoded position (\circ)');
title('Likelihood decoding');


% finally - histogram of decoding errors, like for the linear decoder
subplot(2,2,4); hold on;
tmperr = pos_all - c_all(:,1);
decoding_err = mod((tmperr+180),360)-180;
clear tmperr;

histogram(decoding_err);
set(gca,'TickDir','out','FontSize',12,'XTick',-180:90:180);
title('Accuracy');
xlabel('Decoding error (\circ)');