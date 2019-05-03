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
% last updated 26 April 2019
% Issues with the weight initialisation and updating ensuring a fixed norm,
% updating arithmetically results in -ve weights for non-firing neurons,
% and +ve weights for firing neurons.
% So let's try keeping the mean for each neuron 0 as well as normalising.
%
% There are other possibilities too: 
%   make the mean of each input 0 (perhaps at each timestep)
%   Not yet tried out. 
%
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
useabs = true ; % wilkl use the absolute value signal
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
LIFthreshold = 1.0 ; % LIF neuron threshold

% weight initialisation
weightnorm = 1 ; % normalised value of weight
weightssupplied = false ; % default is to randomly initialise weights
% weight adaptation
k_fired = 0.01 ; % adaptation learning rate for fired neuron
k_notfired = 0.001 ; % adaptation learnking rate for non-firing neurons
singleupdate = true ; % update only first firing neuron at any time step

% varargin parameter setting
i = 1 ;
while(i<=size(varargin,2))
    switch lower(varargin{i})
        case 'debug' 
            debug = varargin{i+1};
            i=i+1 ;
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
        case 'useabs'
            useabs = varargin{i+1}; % will we use the absolute value  signal?
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
        case 'lifthreshold'
            LIFthreshold = varargin{i+1}; % LIF neuron threshold
            i=i+1 ;
        case 'weightssupplied'
            weightssupplied = true ;
            weightarray  = varargin{i+1}; % supplied weights
            i=i+1 ;
        case 'weightnorm'
            weightnorm = varargin{i+1}; % -weightmaxval to weightmaxval is weight range
            i=i+1 ;
        case 'k_fired'
            k_fired = varargin{i+1}; % learning rate for fired neuron
            i=i+1 ;
        case 'k_notfired'
            k_notfired = varargin{i+1}; % learning rate for neurons not firing
            i=i+1 ;
        case 'singleupdate'
            singleupdate = varargin{i+1}; % if true update only first neuron to fire at any time step
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
if ~(useabs || useonset || useoffset)
    % no inputs are set to be used
    disp('spectrotemporal: no input sets selected') ;
    return  ;
end
% calculate size of weight array, and what goes where
ranges = 0 ;
if useabs
    absrange = 1:N;
    ranges = ranges + 1 ;
end
if useonset
    onsetrange = (ranges*N + 1): (ranges+1)*N  ;
    ranges = ranges + 1 ;
end
if useoffset
    offsetrange = (ranges*N + 1): (ranges+1)*N  ;
    ranges = ranges + 1 ;
end

% initialise the weights (or were they supplied in varargin?)
if ~(weightssupplied) % default
    % initialise weights randomly within some range
    weightarray = random('uniform', -1, +1, [M, N*ranges, K]) ;
    %  updated 26 4 19: normalise weight array so that  mean (i, :, :) = 0 and then norm(i, :, :) = weightnorm
    for neuron = 1: M
        weightarray(neuron, :, :) = weightarray(neuron, :, :) - mean(weightarray(neuron, :, :), 'all') ;
        weightarray(neuron, :, :) = weightnorm * (weightarray(neuron, :, :)/norm(squeeze(weightarray(neuron, :, :)))) ;
    end
else % weights were supplied
    % check dimensions of weightarray
    if ~(isequal(size(weightarray), [M, N*ranges, K]))
        disp(['spectrotemporal: supplied weight array and size of bandpass filter, neural array, time lapse do not agree', ...
            num2str(size(weightarray,1)), ' ' , num2str(size(weightarray,2)), ' ', num2str(size(weightarray,3))]) ;
        return ;
    end
end
%put ref period into timesteps
rptimesteps = ceil(LIFrp/LIFtimestep) ; % rp in timesteps

% read input_filelist to get the list of files to be processed
inputfid = fopen([SD '/' filelistfile]) ;
% read 1st line
fline = fgetl(inputfid) ;
nooffiles = 1 ;
while (ischar(fline) && (~isempty(fline)))
    filelist{nooffiles} = fline ; % allow this warning about sizing. 
    fline = fgetl(inputfid) ;
    nooffiles = nooffiles + 1 ;
end
nooffiles = nooffiles - 1 ;
fclose(inputfid) ;

% for each file...
% process each file, one by one

for i = 1:nooffiles
    % read the sound and bandpass the signal (gammatone filterbank)
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
    % could add other transforms as well
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
    if useabs
        r_absSig = resample(absSig', 1, floor(LIFtimestep * fs) )' ;
        % normalise: adjust all to mean of 1
        r_absSig = r_absSig/mean(r_absSig, 'all') ;
    end
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
    LIFrefpdleft = zeros([1 M])  ; % no RP right now
    
    % run the network, updating the weights for firing neurons
    
    % number of timesteps to run is length of sound , in LIFtimeteps, less K 
    % temporary
    if debug
        LIFactivities = zeros([M size(r_absSig, 2) - K]) ; % record neuron activities
        LIFfirings = zeros([M size(r_absSig, 2) - K]) ; % record neuron firings
    end
    for ts = 1:(size(r_absSig, 2) - K) % for each possible time step - without overlapping the end of the signal
        % calculate the incoming activity for each of the neurons
        % W matrix is M (neurons) by B (channels) by K (timesteps), input
        % is timesteps by B; output is M by 1
        newLIFactivity = zeros([M 1]) ;
        LIFfiring = zeros([1 M])  ; % nothing firing
        for tt = 1:K
            if useabs
                newLIFactivity = newLIFactivity + squeeze(weightarray(:,absrange,tt)) * r_absSig(:, ts + tt -1) ;
            end
            if useonset
                newLIFactivity = newLIFactivity + squeeze(weightarray(:,onsetrange,tt)) * r_onset_signal(:, ts + tt -1) ;
            end
            if useoffset
                newLIFactivity = newLIFactivity + squeeze(weightarray(:,offsetrange,tt)) * r_offset_signal(:, ts + tt -1) ;
            end
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
        % decrement refractory period 
        inperiod = find(LIFrefpdleft > 0) ;
        LIFrefpdleft(inperiod) = LIFrefpdleft(inperiod) - 1 ;
        % What's firing?        
        % find where LIFactivity exceeds threshold
        firingset = find(LIFactivity > LIFthreshold) ;
        % reset activity
        LIFactivity(firingset) = 0 ; % should be done throughout - or at end of - refractory period. 

        for firingno = 1:length(firingset) % process firing neurons
            if (LIFrefpdleft(firingset(firingno)) == 0) % omit refractory period neurons
                LIFfiring(firingset(firingno)) = 1 ; % mark firing neurons: exclude thore in refractory period
                LIFrefpdleft(firingset(firingno)) = rptimesteps ; % set refractory period timer
            end
        end

        % temporary: record firing neurons
        if debug
            LIFfirings(:,ts) = LIFfiring ;
        end
        % update weights: LIFfiring contains the neurons that fired this ts
        if (any(LIFfiring)) % do nothing if nothing fires
            if debug
                   disp(['spectrotemporal: at ' num2str(ts) ' firing ' num2str(find(LIFfiring))]) ;
            end
            for nno = 1:M
                if LIFfiring(nno)
                    % neuron fired
                    if useabs
                        weightarray(nno, absrange, :) = weightarray(nno, absrange, :) + reshape(k_fired *  r_absSig(:, (ts + 1):(ts+K)), [1 N K]) ;
                    end
                    if useonset
                        weightarray(nno, onsetrange, :) = weightarray(nno, onsetrange, :) + reshape(k_fired * r_onset_signal(:, (ts + 1):(ts+K)), [1 N K]) ;
                    end
                    if useoffset
                        weightarray(nno, offsetrange, :) = weightarray(nno, offsetrange, :) + reshape(k_fired * r_offset_signal(:, (ts + 1):(ts+K)), [1 N K]) ;
                    end
                    % renormalise
                    weightarray(nno, :, :) = weightarray(nno, :, :) - mean(weightarray(nno, :, :), 'all') ; % mean 0
                    weightarray(nno, :, :) = weightnorm * (weightarray(nno, :, :)/norm(squeeze(weightarray(nno, :, :)))) ; % norm = weightnorm
                else
                    % neuron did not fire
                    if useabs
                        weightarray(nno, absrange, :) = weightarray(nno, absrange, :) - reshape(k_notfired *  r_absSig(:, (ts + 1):(ts+K)), [1 N K]) ;
                    end
                    if useonset
                        weightarray(nno, onsetrange, :) = weightarray(nno, onsetrange, :) - reshape(k_notfired * r_onset_signal(:, (ts + 1):(ts+K)), [1 N K]) ;
                    end
                    if useoffset
                        weightarray(nno, offsetrange, :) = weightarray(nno, offsetrange, :) - reshape(k_notfired * r_offset_signal(:, (ts + 1):(ts+K)), [1 N K]) ;
                    end
                    % renormalise
                    weightarray(nno, :, :) = weightarray(nno, :, :) - mean(weightarray(nno, :, :), 'all') ; % mean 0
                    weightarray(nno, :, :) = weightnorm * (weightarray(nno, :, :)/norm(squeeze(weightarray(nno, :, :)))) ;
                end
                if (singleupdate && LIFfiring(nno)) 
                    break ; % only upodate the firing from a single neuron
                end
            end
        end
        
    end
    
end

% return the final weight matrix.

end

