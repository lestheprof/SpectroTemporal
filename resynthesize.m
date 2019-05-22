function outsignal = resynthesize(params,inputarray, varargin)
%resynthesize resynthesises a sound from the weight array using the
% parameters supplied
% params is the parameters used for creating the weight array. Not all
% of them are necessary
% inputarray is a no of bands by no of timesteps array
%
% started LSS 20 May 2019
%
outFs = 44100 ; % sample rate for the output sound
basetype = 1 ; % 1 is sine wave base type; 2 is noise base type
% varargin parameter setting
i = 1 ;
while(i<=size(varargin,2))
    switch lower(varargin{i})
        case 'outfs'
            outFs = varargin{i+1};
            i=i+1 ;
        case 'basetype'
            basetype = varargin{i+1};
            i=i+1 ;
        otherwise
            error('SpectTempo: Unknown argument %s given',varargin{i});
    end
    i=i+1 ;
end
% create the baseline signal which will be modulated to produce the final
% signal
[nbands, nsegments] = size(inputarray) ;
% sanity check
if ~(params.N == nbands)
    error(['Number of bands in array = ' num2str(nbands) ' differs from value in paramaters = ' num2str(params.N)]) ;
end
sigduration = nsegments * params.LIFtimestep ; % signal duration
sigsamples = ceil(sigduration * outFs) ;
cochCFs=MakeErbCFs(params.minCochFreq,params.maxCochFreq,params.N); % get centre frequencies
if (basetype == 1) % sin wave base signal
    rawsignal = sin((2 * pi * cochCFs)' * ((0:sigsamples-1) / outFs)) ; % N (bands) times samples in each band
else
    if (basetype == 2) % noise base type
        % generatwe white nboise of appropriate length
        whitenoise = wgn(sigsamples, 1, 0) ;
        % bandpass filter the noise
        rawsignal = zeros([params.N sigsamples]) ;
        for chno=1:1:params.N  % which channel
            rawsignal(chno, :) = gammatone1(whitenoise,outFs,cochCFs(chno),params.N_erbs);
            % normalise to max value of 1
            rawsignal(chno, :) = rawsignal(chno, :) / max(abs(rawsignal(chno, :))) ;
        end
    end
end
        
outsignal = zeros(size(rawsignal)) ;

% modulate the signal, band by band (can I do it all at once?)
for bno = 1:nbands
    modulator = spline( (1:size(inputarray,2)) * params.LIFtimestep, inputarray(bno,:), (1:size(rawsignal, 2))/outFs) ;
    outsignal(bno,:) = modulator .* rawsignal(bno,:) ;
end
outsignal = sum(outsignal) ; % add up the different bands
end

