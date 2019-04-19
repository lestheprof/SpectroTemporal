function weightarray = spectrotemporal(SD,filelistfile, varargin)
%spectrotemporal finds the weights that arise from this set of sound
%inputsto a set of LIF neurons
%   creates a 3D array of weights, weightarray. The files named in
%   filelist in directory sd are used as input. N is the number of bandpassed channels,
%   M the number of neurons, and K the numnber of timesteps used. The
%   returned array is M by N by K. All the other parameters are set in
%   varargin, having had a default value associated first.
%
% started LSS 16 April 2019
% last updated 17 April 2019
debug = true ;
% defaults
N = 100 ; % number of bandpass channels
M = 50 ; % number of LIF neurons ;
K = 30 ; % number of timesetps used
% Filterbank parameters
Fs = 44100 ; % sampling rate
minCochFreq = 200 ; % minimum gammatone frequency
maxCochFreq = 5000 ; %  maximum gammatone freuency
N_erbs = 1 ; % default bandwidth
smoothlength = 0.01 ; % Bartlett filter parameter length of triangular (bartlett) window used to smooth
%   rectified signal
% parameters for onset and offset signal generation (if used)
sigma1 = 0.02 ; %half difference of Gaussians std
sigmaratio = 1.2 ; % ratio of the two HDOG's
nsamples = 4000 ; % number of samples used in onset/offset convolving function
% LIF parameters
dissipation = 100 ; % dissipation, 1/tau
rp = 0.002 ; % refractory period
% weight initialisation
weightnorm = 1 ; % normalised value of weight
weightssupplied = false ; % default is to randomly initialise weights
% weight adaptation
k_fired = 0.1 ; % adaptation learning rate for fired neuron
k_notfired = 0.01 ; % adaptation learnking rate for non-firing neurons

% varargin parameter setting
i = 1 ;
while(i<=size(varargin,2))  
    switch lower(varargin{i})
        case 'fs' % sampling rate
            Fs = varargin{i+1};
            i=i+1 ;
       case 'mincochfreq' % minimum gammatone frequency
            minCochFreq = varargin{i+1};
            i=i+1 ;
        case 'maxcochfreq' %  maximum gammatone freuency
            maxCochFreq = varargin{i+1};
            i=i+1 ;
        case 'n_erbs' % default bandwidth
            N_erbs = varargin{i+1};
            i=i+1 ;
        case 'n'
            N = varargin{i+1}; % number of bandpass channels
            i=i+1 ;
        case 'smoothlength'
            smoothlength = varargin{i+1}; % smoothing for signal
            i=i+1 ;
        case 'sigma1'
            sigma1 = varargin{i+1}; % HDoG std (in seconds)
            i=i+1 ;
        case 'sigmaratio'
            sigmaratio  = varargin{i+1}; % ratio of std's
            i=i+1 ;
        case 'nsamples' % number of samples in colvolving function
            nsamples = varargin{i+1};
            i=i+1 ;
        case 'm'
            M = varargin{i+1}; % number of LIF neurons
            i=i+1 ;
        case 'dissipation'
            dissipation = varargin{i+1}; % dissipation for all LIF neurons
            i=i+1 ;
        case 'rp'
            rp = varargin{i+1}; % refractory period (in seconds)
            i=i+1 ;
        case 'weightssupplied'
            weightssupplied = true ;
            weightarray  = varargin{i+1}; % supplied weights
            i=i+1 ;
        case 'weightnorm'
            weightnorm = varargin{i+1}; % normalisation for weight
            i=i+1 ;
        case 'k_fired'
            k_fired = varargin{i+1}; % learning rate for fired neuron
            i=i+1 ;
        case 'k_notfired'
            k_notfired = varargin{i+1}; % learning rate for neurons not firing
            i=i+1 ;

         otherwise
            error('spectrotemporal: Unknown argument %s given',varargin{i});
    end
    i=i+1 ;
end

% outline of processing
% precomputation of data used more than once

% Wholegaussian and smoothing should probably be precomputed if there are many signals to be
% processed
dtperelement = 1.0/Fs ; % sampling interval
wholegaussian = diffofgaussians(sigma1, sigma1 * sigmaratio, nsamples * 2 + 1,dtperelement) ;
hdog = wholegaussian(nsamples: end) ;
if (debug)
    disp(['spectrotemporal: sum of half difference of Guassians = ' num2str(sum(hdog))]) ;
end
% Calculate smoothing for signal
bartlettlength = floor(smoothlength/dtperelement) ;
bartlettwindow = bartlett(bartlettlength) ;
bartlettwindow = bartlettwindow/sum(bartlettwindow) ; % normalise to sum of 1

% initialise the weights (or were they supplied in varargin?)
if ~(weightssupplied) % default
    weightarray = zeros([M, N, K]) ;
else % weights were supplied
    % check dimensions of weightarray
    if ~(isequal(size(weightarray), [M, N, K]))
        disp(['spectrotemporal: supplied weight array and size of bandpass filter, neural array, time lapse do not agree', ...
            num2str(size(weightarray,1)), ' ' , num2str(size(weightarray,2)), ' ', num2str(size(weightarray,3))]) ;
        return ;
    end
end

% read input_filelist to get the list of files to be processed
inputfid = fopen([SD '/' filelistfile]) ;
% inputfid = fopen(filelist) ;
fline = fgetl(inputfid) ;
nooffiles = 1 ;
while (ischar(fline) && (~isempty(fline)))
    filelist{nooffiles} = fline ;
    fline = fgetl(inputfid) ;
    nooffiles = nooffiles + 1 ;
end
nooffiles = nooffiles - 1 ;
fclose(inputfid) ;

% for each file...% process each file, one by one
for i = 1:nooffiles
% read the sound and bandpass the signal
% calculate the smoothed absolute signal for each band
% possibly calculate oosignal for each band
% possibly calculate the onset signal for each band
% possibly calculate the offset signal for each band
% (use one of the above)
% initialise the network at the start of the sound
% run the network, updating the weights

end

% return the final weight matrix.

end

