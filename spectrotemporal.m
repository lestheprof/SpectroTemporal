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
MAXDURATION = 100 ; % maximal duration of a single sound file
smoothlength = 0.01 ; % Bartlett filter parameter length of triangular (bartlett) window used to smooth
%   rectified signal
useonset = false ; % will we use onset signals in processing?
useoffset = false ; % will we use offset signals in processing?
% parameters for onset and offset signal generation (if used)
sigma1 = 0.01 ; %half difference of Gaussians std
sigmaratio = 1.2 ; % ratio of the two HDOG's
nsamples = 4000 ; % number of samples used in onset/offset convolving function
logabs = false ; % use log of absolute value
logonset = false ; % use log of onset/offset values
% LIF parameters
LIFtimestep = 0.001 ; % timestep for use with LIF network
LIFdissipation = 100 ; % dissipation, 1/tau
LIFrp = 0.002 ; % refractory period
rptimesteps = ceil(LIFrp/LIFtimestep) ; % rp in timesteps
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
        case 'maxduration'
            MAXDURATION = varargin{i+1}; % maximal duration of a single sound file
            i=i+1 ;
        case 'smoothlength'
            smoothlength = varargin{i+1}; % smoothing for signal
            i=i+1 ;
        case 'useonset'
            useonset = varargin{i+1}; % will we use the onset signal?
            i=i+1 ;
        case 'useoffset'
            useoffset = varargin{i+1}; % will we use the offset signal?
            i=i+1 ;
        case 'logabs'
            logabs = varargin{i+1}; % do we take the log of the absolute value?
            i=i+1 ;
        case 'logonset'
            logonset = varargin{i+1}; % do we use the log of the onset/offset values
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
        case 'liftimestep'
            LIFtimestep = varargin{i+1}; % timestep for use with LIF neurons
            i=i+1 ;
        case 'lifdissipation'
            LIFdissipation = varargin{i+1}; % dissipation for all LIF neurons
            i=i+1 ;
        case 'lifrp'
            LIFrp = varargin{i+1}; % refractory period (in seconds)
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
    % initialise weights randomly within some range
    weightarray = random('uniform', -1, +1, [M, N, K]) ;
    % normalise weight array so that norm(i, :, :) = 1
    for neuron = 1: M
        weightarray(neuron, :, :) = weightarray(neuron, :, :)/norm(squeeze(weightarray(neuron, :, :))) ;
    end
    % weightarray = zeros([M, N, K]) ; % currently assumes only abs input, otherwise needs another K for onset, and another for offset
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
    [bmSig, sig, fs, datalength, cochCFs, delayVector] = ...
        bmsigmono([SD '/' filelist{i}], N, minCochFreq, maxCochFreq, MAXDURATION, 'gamma', N_erbs) ;
    if (Fs ~= fs)
        disp(['spectrotemporal: provided sample rate = ' num2str(Fs) ' differs from sound sample rate = ' num2str(fs)]) ;
    end
    % calculate the smoothed absolute signal for each band
    absSig = zeros(size(bmSig)) ; % initialise
    if (useonset || useoffset) % if either onset or offset signal is to be used we'll also need the onset/offset signal
        ooSig = zeros(size(bmSig)) ; % initialise
    end
    for band = 1:N % calculate the abs, and possibly onset and offset signals for each band
        absSig(band,:) = conv(abs(bmSig(band,:)),bartlettwindow, 'same') ;
        if (useonset || useoffset) % only calculate if required
            ooSig(band,:) = conv(absSig(band,:),hdog, 'same') ;
        end
    end
    % possibly calculate the onset signal for each band
    % possibly calculate the offset signal for each band
    if useonset
        onset_signal = abs(max(0,ooSig)) ; % onset signal, positive or 0, 1 per band
    end
    if useoffset
        offset_signal = abs(max(0, -ooSig)) ; % offset signal, positive or 0, 1 per band
    end
    if (useonset || useoffset)
        if logonset % logarithmic adjustment?
            onset_signal = log(1 + onset_signal) ;
            offset_signal = log(1 + offset_signal) ;
        end
    end
    if logabs
        absSig = log(absSig) ;
    end
    % repackage the inputs into buckets of length LIFtimestep
    % always do absSig
    r_absSig = resample(absSig', 1, floor(LIFtimestep * fs) )' ;
    % normalise: adjust all to mean of 1
    r_absSig = r_absSig/mean(r_absSig, 'all') ;
    if useonset % only do onset if it's to be used
        r_onset_signal = resample(onset_signal', 1, floor(LIFtimestep * fs) )' ;
        r_onset_signal = r_onset_signal/mean(r_onset_signal, 'all') ;
    end
    if useoffset % only do offset if it's to be used
        r_offset_signal = resample(offset_signal', 1, floor(LIFtimestep * fs) )' ;
        r_offset_signal = r_offset_signal/mean(r_offset_signal, 'all') ;
    end
    
    % initialise the network at the start of each sound
    % initialise LIF neurons: each is represented by an array representing the
    % activity of each neuron, an array representing whether the neuron is
    % firing, (and if we use a refractory period, an array representing the
    % amount of time left in thar RP)
    LIFactivity = zeros([1 M]) ; % currently initialised to 0
    LIFfiring = zeros([1 M])  ; % nothing firing
    LIFrefpdleft = zeros([1 M])  ; % no RP right now
    
    % run the network, updating the weights for firing neurons
    
    % number of timesteps to run is length of sound , in LIFtimeteps, less K 
    % temporary
    if debug
        LIFactivities = zeros([M size(r_absSig, 2) - K]) ; % record neuron activities
    end
    for ts = 1:(size(r_absSig, 2) - K)
        % calculate the incoming activity for each of the neurons
        % W matrix is M (neurons) by B (channels) by K (timesteps), input
        % is timesteps by B; output is M by 1
        newLIFactivity = zeros([M 1]) ;
        for tt = 1:K
            newLIFactivity = newLIFactivity + squeeze(weightarray(:,:,tt)) * r_absSig(:, ts + tt -1) ;
        end
       
        % calculate current activity: should the newLIFactivity be *ed by
        % LIFtimestep?
        LIFactivity = LIFactivity + LIFtimestep * (newLIFactivity' -  LIFdissipation * LIFactivity) ;
        % keep it nonnegative
        LIFactivity = max(0, LIFactivity) ;

        % temporary: record neuron activities
        if debug
             LIFactivities(:, ts) = LIFactivity ;
        end
        % What's firing?
        % update weights
    end
    
end

% return the final weight matrix.

end

