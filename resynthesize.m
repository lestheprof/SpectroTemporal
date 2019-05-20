function outsignal = resynthesize(params,inputarray)
%resynthesize resynthesises a sound from the weight array using the
% parameters supplied
% params is the parameters used for creating the weight array. Not all
% of them are necessary
% inputarray is a no of bands by no of timesteps array
%
% started LSS 20 May 2019
%
outFs = 44100 ; % sample rate for the output sound

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
rawsignal = sin((2 * pi * cochCFs)' * ((0:sigsamples-1) / outFs)) ; % N (bands) times samples in each band
outsignal = zeros(size(rawsignal)) ;

% modulate the signal, band by band (can I do it all at once?)
for bno = 1:nbands
    modulator = spline( (1:size(inputarray,2)) * params.LIFtimestep, inputarray(bno,:), (1:size(rawsignal, 2))/outFs) ;
    outsignal(bno,:) = modulator .* rawsignal(bno,:) ;
end
outsignal = sum(outsignal) ; % add up the different bands
end

