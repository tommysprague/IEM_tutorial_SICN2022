function bb = build_basis_polar_mat(eval_at,chan_th,chan_size,chan_pow )
% build_basis_polar_mat creates a basis set of selective 1-d information
% channels in a 360-deg space, computed along the points in EVAL_AT, with
% channels centered at CHAN_TH, with size given by CHAN_SIZE, and raised to
% the power indicated by CHAN_POW
%
% bb = build_basis_polar_mat(eval_at,chan_th) computes a basis set with the
% default size (180 deg) and the default power (defined by # of basis
% functions per Freeman & Adelson, 1992).
%
% bb = build_basis_polar_mat(eval_at,chan_th,chan_size) uses the size
% specified by chan_size (one number; positive; deg) for all channels
%
% bb = build_basis_polar_mat(eval_at,chan_th,[],chan_pow) uses chan_pow
% instead of the Freeman & Adelson recommendation. This will alter the
% shape of the channel profile slightly, and can decrease its FWHM.
%
% CHAN_WIDTH adjusts a 'size' constant - the distance from the channel
% center (in polar angle deg) at which the filter reaches zero. This will
% scale linearly with FWHM, but the ratio will depend on the power. 
%
% the basis set is returned as a matrix, with each column containing one
% basis function/modeled information channel. The resulting dimensions are
% length(eval_at) x length(chan_th)
%
% Example:
% my_angs = linspace(-179,180,360);
% my_chan_centers = linspace(360/8,360,8);
% my_chan_width = 180;
%
% my_basis = build_basis_polar_mat(my_angs,my_chan_centers,my_chan_width);
% figure; plot(my_angs,my_basis,'-','LineWidth',1.5);
% xlabel('Position (\circ)'); ylabel('Channel sensitivity');
%
% Tommy Sprague (tommy.sprague@gmail.com)


% adapted from build_basis_polar - this does everything at once, useful for
% generating a set of transformed/rotated basis functions

% will return bb, length(evalAt) x length(rfTh)
n_basis = length(chan_th);
[chan_th,eval_at] = meshgrid(squeeze(chan_th),squeeze(eval_at));

% TODO: make sure only one chan_size...


if nargin < 3 || isempty(chan_size)
    chan_size = 180; % default value
end

if nargin < 4 || isempty(chan_pow)
    chan_pow = n_basis-mod(n_basis,2);
end

% utility function to compute distance between two angles
ang_dist = @(a,b) min(mod(a-b,360),mod(b-a,360));


bb = (cosd( 180 * ang_dist(eval_at,chan_th) ./ (2*chan_size) ).^chan_pow) .* (ang_dist(eval_at,chan_th)<=chan_size) ;


return