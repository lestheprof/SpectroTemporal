function soundout = resynthesize(params,inputarray)
%resynthesize resynthesises a sound from the weight array using the
% parameters supplied
%   Detailed explanation goes here
% params is the parameters used for creating the weight array. Not all
% of them are necessary
% inputarray is a no of bands by no of timesteps array
outFs = 44100 ; % sample rate for the output sound

% create the basline singal which will be modulated to produce the final
% signal
[nbands, nsegments] = size(inputarray) ;
% sanity check
if ~(params.N == nbands)
    error(['Number of bands in array = ' num2str(nbands) ' differs from value in paramaters = ' num2str(params.N)]) ;
end
sigduration = nsegments * params.LIFtimestep ; % signal duration
sigsamples = ceil(sigduration * outFs) ;
cochCFs=MakeErbCFs(params.minCochFreq,params.maxCochFreq,params.N); % get centre frequencies
rawsignal = sin(2 * pi * cochCFs * ([1:sigsamples] / outFs)) ;

% 

end

