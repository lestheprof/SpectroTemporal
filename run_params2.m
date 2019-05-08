
[x, hname] = system('hostname');
hname = deblank(hname);
switch hname
    case 'peseta.cs.stir.ac.uk'
        f101dir = '/Volumes/Extra/AllenCorpus/AllenCorpus/f101' ;
        m102dir = '/Volumes/Extra/AllenCorpus/AllenCorpus/m102' ;
    case 'Leslie-Smiths-MacBook-Pro.local'
        f101dir = '../SoundStimuli_2019/f101' ;
        m102dir = '../SoundStimuli_2019/m102' ;
    otherwise
        disp('Computer name not known: add paths to case statement') ;
end

% from using params1.m, <2, 1> seemed good: that is,
% kf = 0.00015, knf = 0.00001.
% this one looks at the time overall. 

LIFrpMultiplier = 10 ;
k_fired = 0.00015 ;
k_notfired = 0.00001 ;
ktimestep = 0 ;
kdissipation = 0 ;
weightnorm = 4 ;
debug = true ;
for LIFtimestep = 0.005:0.01:0.025
    ktimestep = ktimestep + 1 ;
    kdissipation = 0 ;
    for LIFdissipation = 20:40:100
        kdissipation = kdissipation + 1 ;
        fname = [num2str(ktimestep) '_' num2str(kdissipation) '.mat'] ;
        if debug
            disp(['LIF_timestep = ' num2str(LIFtimestep) 'LIF_dissipation = ' num2str(LIFdissipation)]) ;
        end
        LIFrp = LIFrpMultiplier * LIFtimestep ;
        f101_2a = spectrotemporal(f101dir, 'filelist_10.txt', 'fs', 16000, 'useabs', true, ...
            'useonset', false, 'useoffset', false, 'weightnorm', weightnorm, 'LIFrp', 10 * LIFtimestep, ...
            'k_fired', k_fired, 'k_notfired', k_notfired, 'M', 5, 'liftimestep', LIFtimestep, ...
            'lifdissipation', LIFdissipation, 'debug', debug) ;
        m102_2a = spectrotemporal(m102dir, 'filelist_10.txt', 'fs', 16000, 'useabs', true, ...
            'useonset', false, 'useoffset', false, 'weightnorm', weightnorm, 'LIFrp', 10 * LIFtimestep, ...
            'k_fired', k_fired, 'k_notfired', k_notfired, 'M', 5, 'liftimestep', LIFtimestep, ...
            'lifdissipation', LIFdissipation, 'debug', debug) ;
        % save files
        save([f101dir '/' 'run2_' fname], 'f101_2', 'k_fired', 'k_notfired', 'LIFtimestep', 'LIFdissipation','weightnorm', 'LIFrp')  ;
        save([m102dir '/' 'run2_' fname], 'm101_2', 'k_fired', 'k_notfired', 'LIFtimestep', 'LIFdissipation','weightnorm', 'LIFrp')  ;
    end
end
