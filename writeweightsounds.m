function  writeweightsounds(filebasename, params, wtmatrix , numtowrite)
%writeweightsounds write a nuymber of wav filkes from the weight matrix
%   Detailed explanation goes here
risetime = 0.005 ;
falltime = 0.005 ;
fs = 44100 ;
for sno = 1: numtowrite
    filedata = resynthesize(params, squeeze(wtmatrix(sno,:,:)),'basetype', 1) ;
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

