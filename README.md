### Description

This package provides an implementation of the Matrix Element Method for the ttH H->tautau analysis.


### Instructions for 8_0_X

```
cmsrel CMSSW_8_0_12
cd CMSSW_8_0_12/src
cmsenv
git clone https://github.com/tstreble/ttH_Htautau_MEM_Analysis
scram b -j 4
```


### Run test code

```
cd CMSSW_8_0_12/bin/slc6_amd64_gcc530
cmsenv
./test ttH_Htautau_MEM_Analysis/MEM/small_nomin_122016.py
```
