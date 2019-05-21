f101dir = '/Volumes/Extra/AllenCorpus/AllenCorpus/f101' ;
m102dir = '/Volumes/Extra/AllenCorpus/AllenCorpus/m102' ;
kf = 1 ;
for k_fired = [0.0001, 0.0005, 0.001 , 0.005, 0.01]
        fname = [num2str(kf) '.mat'] ;
        [wnew5allf101 paramsf101] = SpectTempo(FD, 'filelist_10.txt', 7, 'Fs', 16000, 'weightnormtype', 2, 'k' ,50,'weightnormsubseq', 16, 'debug', 0, 'liftimestep', 0.0025, 'k_fired', k_fired, 'debug', 0) ;
        [wnew5allm102 paramsm102] = SpectTempo(FD, 'filelist_10.txt', 7, 'Fs', 16000, 'weightnormtype', 2, 'k' ,50,'weightnormsubseq', 16, 'debug', 0, 'liftimestep', 0.0025, 'k_fired', k_fired, 'debug', 0) ;

        % save files
        save([f101dir '/' 'may212019_' fname], 'wnew5allf101', 'paramsf101')  ;
        save([m102dir '/' 'may212019_' fname], 'wnew5allm102', 'paramsm102')  ;
       
    kf = kf + 1 ;
    disp(['k_fired = ' num2str(k_fired)]) ;
end
