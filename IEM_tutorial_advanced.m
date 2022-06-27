% IEM_tutorial_advanced.m
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
% - implementation of IEM
% - Decoding stimulus feature values from reconstructed channel response
%   profiles - we'll do this with a simple population vector model, and with
%   a Bayesian extension which implements a simplified noise model to recover
%   a likelihood function over stimulus feature space
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