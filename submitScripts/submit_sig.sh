#perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/DYJets_mumuFilter_test.txt isMC=True --prodSpace /data/cmszfs1/user/revering/ZmmSim --batch 10 --jobname ZmmSim_testBarrelCheck
perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/0p2sig_gen_data.txt isMC=True isSig=True --prodSpace /local/cms/user/auerb029/Run3Analyzed7_20_23 --batch 10 --jobname 0p2signal_analysis

perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/0p4sig_gen_data.txt isMC=True isSig=True --prodSpace /local/cms/user/auerb029/Run3Analyzed7_20_23 --batch 10 --jobname 0p4signal_analysis

perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/0p6sig_gen_data.txt isMC=True isSig=True --prodSpace /local/cms/user/auerb029/Run3Analyzed7_20_23 --batch 10 --jobname 0p6signal_analysis

perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/0p8sig_gen_data.txt isMC=True isSig=True --prodSpace /local/cms/user/auerb029/Run3Analyzed7_20_23 --batch 10 --jobname 0p8signal_analysis

perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/1p0sig_gen_data.txt isMC=True isSig=True --prodSpace /local/cms/user/auerb029/Run3Analyzed7_20_23 --batch 10 --jobname 1p0signal_analysis


