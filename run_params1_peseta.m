f101dir = '/Volumes/Extra/AllenCorpus/AllenCorpus/f101' ;
m102dir = '/Volumes/Extra/AllenCorpus/AllenCorpus/m102' ;
kf = 1 ;
for k_fired = 0.0001:0.00005:0.0003
    knf = 1 ;
    for k_notfired = 0.00001:0.00001:0.00006
        fname = [num2str(kf) '_' num2str(knf) '.mat'] ;
        f101_2 = spectrotemporal(f101dir, 'filelist_all.txt', 'fs', 16000, 'useabs', true, ...
            'useonset', false, 'useoffset', false, 'weightnorm', 2, 'LIFrp', 0.025, ...
            'k_fired', k_fired, 'k_notfired', k_notfired, 'M', 5, 'liftimestep', 0.0025, ...
            'lifdissipation', 20, 'debug', false) ;
        m102_2 = spectrotemporal(m102dir, 'filelist_all.txt', 'fs', 16000, 'useabs', true, ...
            'useonset', false, 'useoffset', false, 'weightnorm', 2, 'LIFrp', 0.025, ...
            'k_fired', k_fired, 'k_notfired', k_notfired, 'M', 5, 'liftimestep', 0.0025, ...
            'lifdissipation', 20, 'debug', false) ;
        % save files
        save([f101dir '/' 'aaarun_' fname], 'f101_2', 'k_fired', 'k_notfired')  ;
        save([m102dir '/' 'aaarun_' fname], 'm101_2', 'k_fired', 'k_notfired')  ;
       
        knf = knf + 1 ;
    end
    kf = kf + 1 ;
end
