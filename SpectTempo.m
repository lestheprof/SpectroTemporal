function [weights, params] = SpectTempo(filedir, filelist, ncomponents,  varargin)
% SpectTempo finds the spectrotemporal components one by one, subtracting
% the back-projection  of those already found.
% aim is an industrial-strength component finder.
% also would be good to put segmentation in as an option.
% also would be good to return a structure with all the paramters for use in
% re-creating the results.
%
% parameters: filedir: directory where the files are to be found
% filelist: name of file with list of files to be run
% ncomponents: the number of compoinents tio be found.
%
% Many other values can be set using varargin

% find the location to be used: redundant now that filedir and filelist are
% paramaters.,
% [x, hname] = system('hostname');
% hname = deblank(hname);
% switch hname
%     case 'peseta.cs.stir.ac.uk'
%         f101dir = '/Volumes/Extra/AllenCorpus/AllenCorpus/f101' ;
%         m102dir = '/Volumes/Extra/AllenCorpus/AllenCorpus/m102' ;
%     case 'Leslie-Smiths-MacBook-Pro.local'
%         f101dir = '../SoundStimuli_2019/f101' ;
%         m102dir = '../SoundStimuli_2019/m102' ;
%     otherwise
%         disp('Computer name not known: add paths to case statement') ;
% end

% default values
debug = true ;
Fs = 16000 ;
N = 100 ; % number of bandpass channels
% M = 1 ; % number of LIF neurons: 1 at a time here: forced in
% spectrotempporal2
K = 30 ; % number of timesetps used
% Filterbank parameters
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
weightnorm1 = 4 ; % normalised value of weight
weightnormsubseq = 8 ;
weightnormtype = 1 ;
weightarray = [] ; % if weightarray stays null, weights will be initialised randomly
% this could allow all the neurons to be initialised to the same weight
% values
% weight adaptation
k_fired = 0.01 ; % adaptation learning rate for fired neuron
posonly = true ; % use only positive values in signals 


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
        case 'k'
            K = varargin{i+1}; % number of timesteps
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
        case 'nsamples' % number of samples in convolving function
            nsamples = varargin{i+1};
            i=i+1 ;
            %         case 'm' % not permitted to change M in this version
            %             M = varargin{i+1}; % number of LIF neurons
            %             i=i+1 ;
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
            weightarray  = varargin{i+1}; % supplied weights
            i=i+1 ;
        case 'weightnorm1'
            weightnorm1 = varargin{i+1}; % -weightmaxval to weightmaxval is weight range
            i=i+1 ;
        case 'weightnormsubseq'
            weightnormsubseq = varargin{i+1}; % -weightmaxval to weightmaxval is weight range
            i=i+1 ;
        case 'weightnormtype'
            weightnormtype = varargin{i+1}; % -weightmaxval to weightmaxval is weight range
            i=i+1 ;
        case 'k_fired'
            k_fired = varargin{i+1}; % learning rate for fired neuron
            i=i+1 ;
        case 'posonly'
            posonly = varargin{i+1}; % if true (default) use only the positive-going parts of signals after feedback
            i=i+1 ;
            
        otherwise
            error('SpectTempo: Unknown argument %s given',varargin{i});
    end
    i=i+1 ;
end

% initialise weights array
weights = zeros([ncomponents N K]) ;
% find 1st component: needs all the varargin values specified
weights(1,:,:) = spectrotemporal2(filedir, filelist, 'fs', Fs, 'mincochfreq',minCochFreq, ...
    'maxCochFreq', maxCochFreq, 'N_erbs', N_erbs, 'N', N, 'K', K, 'maxduration', MAXDURATION, ...
    'smoothlength', smoothlength, 'useabs', useabs, 'logabs', logabs, 'logonset', logonset, ...
    'useonset', useonset, 'useoffset', useoffset, 'weightssupplied', weightarray, 'weightnorm', weightnorm1, 'LIFrp', LIFrp, ...
    'k_fired', k_fired, 'liftimestep', LIFtimestep,'LIFthreshold', LIFthreshold,  ...
    'lifdissipation', LIFdissipation, 'posonly', posonly, 'weightnormtype', weightnormtype, 'debug', debug) ;
for compno = 2:ncomponents
    weights(compno,:,:) = spectrotemporal2(filedir, filelist, 'fs', Fs, 'mincochfreq',minCochFreq, ...
    'maxCochFreq', maxCochFreq, 'N_erbs', N_erbs, 'N', N, 'K', K, 'maxduration', MAXDURATION, ...
    'smoothlength', smoothlength, 'useabs', useabs, 'logabs', logabs, 'logonset', logonset, ...
    'useonset', useonset, 'useoffset', useoffset, 'weightssupplied', weightarray, 'weightnorm', weightnormsubseq, 'LIFrp', LIFrp, ...
    'k_fired', k_fired, 'liftimestep', LIFtimestep,'LIFthreshold', LIFthreshold,  ...
    'lifdissipation', LIFdissipation, 'posonly', posonly, 'debug', debug, 'weightnormtype', weightnormtype, 'existingweights', weights(1:compno-1,:,:)) ;
    
end
% create a structure with all the parameters for returning
params.filedir = filedir ;
params.filelist = filelist ;
params.Fs = Fs ;
params.minCochFreq = minCochFreq ;
params.maxCochFreq = maxCochFreq ;
params.N_erbs = N_erbs ;
params.N = N ;
params.K = K ;
params.maxduration = MAXDURATION ;
params.smoothlength = smoothlength ;
params.useabs = useabs ;
params.logabs = logabs ;
params.logonset = logonset ;
params.useonset = useonset ;
params.useoffset = useoffset ;
params.weightnorm1 = weightnorm1 ;
params.LIFrp = LIFrp ;
params.k_fired = k_fired ;
params.LIFtimestep = LIFtimestep ;
params.LIFthreshold = LIFthreshold ;
params.LIFdissipation = LIFdissipation ;
params.posonly = posonly ;
params.weightnormtype = weightnormtype ;


end