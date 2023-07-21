#perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/DYJets_mumuFilter_test.txt isMC=True --prodSpace /data/cmszfs1/user/revering/ZmmSim --batch 10 --jobname ZmmSim_testBarrelCheck
perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/bkg_gen_data.txt isMC=True --prodSpace /local/cms/user/auerb029/Run3Analyzed7_20_23 --batch 10 --jobname background_analysis
