ALL_TOOLS      += self
self_EX_INCLUDE := /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/combinetool/CMSSW_10_2_1/src /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/combinetool/CMSSW_10_2_1/include/slc6_amd64_gcc700/src /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/combinetool/CMSSW_10_2_1/include/LCG /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_1/src /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_1/include/slc6_amd64_gcc700/src /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_1/include/LCG
self_EX_LIBDIR := /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/combinetool/CMSSW_10_2_1/lib/slc6_amd64_gcc700 /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/combinetool/CMSSW_10_2_1/external/slc6_amd64_gcc700/lib /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_1/lib/slc6_amd64_gcc700 /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_1/external/slc6_amd64_gcc700/lib
self_EX_FLAGS_SYMLINK_DEPTH_CMSSW_SEARCH_PATH  := 2
self_EX_FLAGS_LLVM_ANALYZER  := llvm-analyzer
self_EX_FLAGS_SKIP_TOOLS_SYMLINK  := cxxcompiler ccompiler f77compiler gcc-cxxcompiler gcc-ccompiler gcc-f77compiler llvm-cxxcompiler llvm-ccompiler llvm-f77compiler llvm-analyzer-cxxcompiler llvm-analyzer-ccompiler icc-cxxcompiler icc-ccompiler icc-f77compiler x11 dpm
self_EX_FLAGS_DEFAULT_COMPILER  := gcc
self_EX_FLAGS_EXTERNAL_SYMLINK  := PATH LIBDIR CMSSW_SEARCH_PATH
self_EX_FLAGS_NO_EXTERNAL_RUNTIME  := LD_LIBRARY_PATH PATH CMSSW_SEARCH_PATH
TOOLS_OVERRIDABLE_FLAGS  +=CPPDEFINES CXXFLAGS FFLAGS CFLAGS CPPFLAGS LDFLAGS CUDA_FLAGS CUDA_LDFLAGS
self_ORDER := 20000
IS_PATCH:=

