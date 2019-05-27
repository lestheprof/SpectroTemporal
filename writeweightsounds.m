function  writeweightsounds(filebasename, params, wtmatrix , numtowrite, varargin)
%writeweightsounds write a nuymber of wav filkes from the weight matrix
%   generate sounds from 1D arrays contained in the 2Darray wtmatrix. 
%  basetype 1 uses sin waves at fc from gammatone filter
%  basetype 2 uses bandpassed limited noise
% risetime and falltime are linear for now.
%
% LSS May 22 2019
%
risetime = 0.005 ;
falltime = 0.005 ;
fs = 44100 ;
n_erbs = 0.25 ; % smaller is smaller band size for noise
basetype = 2 ; % bandlimited noise is default
% varargin parameter setting
i = 1 ;
while(i<=size(varargin,2))
    switch lower(varargin{i})
        case 'risetime'
           risetime = varargin{i+1}; % output sound rise time
            i=i+1 ; 
        case 'falltime'
            falltime = varargin{i+1}; % output sound fall time
            i=i+1 ; 
        case 'n_erbs'
            n_erbs = varargin{i+1}; % noise band size adjustment
            i=i+1 ; 
        case 'basetype'
            basetype = varargin{i+1}; % base type for resynthesis
            i=i+1 ;  
        otherwise
            error('writeweightsounds: Unknown argument %s given',varargin{i});
    end
    i=i+1 ;
end
for sno = 1: numtowrite
    filedata = resynthesize(params, squeeze(wtmatrix(sno,:,:)),'basetype', basetype, 'n_erbs', n_erbs) ;
    % normalise filedata to max value of 0.9
    filedata = (0.9/max(abs(filedata))) * filedata ;
    % modulate envelope
    risefallmodulator = ones([1 length(filedata)]) ;
    for sampno = 1: floor(risetime * fs)
        risefallmodulator(sampno) = (sampno/floor(risetime * fs)) ; % ramp up
    end
    fallsamp = floor(falltime * fs) ;
    for sampno = (length(filedata) - floor(falltime * fs)) : length(filedata)
        risefallmodulator(sampno) = fallsamp/floor(falltime * fs) ;
        fallsamp = fallsamp - 1 ;
    end
    filedata = risefallmodulator .* filedata ;    
    audiowrite([filebasename '_' num2str(sno) '.wav'], filedata, fs) ;
end

end

