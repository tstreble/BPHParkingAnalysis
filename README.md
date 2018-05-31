# BPHParkingAnalysis
Analysis package for analysis of BPH parked data

### Instructions for 10_1_5
```
cmsrel CMSSW_10_1_5
cd CMSSW_10_1_5/src/
cmsenv

git clone https://github.com/ICBPHCMS/BPHParkingAnalysis.git
cd NtupleProducer
git checkout master
```

### Producing ntuples
```
sh Make.sh BToKpipiNtupleProducer.cc
./BToKpipiNtupleProducer.exe mc --input test94X_NANO.root --output test.root
```
