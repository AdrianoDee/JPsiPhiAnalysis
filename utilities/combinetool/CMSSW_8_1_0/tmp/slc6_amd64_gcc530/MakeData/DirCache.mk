ALL_PACKAGES += HiggsAnalysis/CombinedLimit
subdirs_src_HiggsAnalysis_CombinedLimit := src_HiggsAnalysis_CombinedLimit_bin src_HiggsAnalysis_CombinedLimit_python src_HiggsAnalysis_CombinedLimit_test src_HiggsAnalysis_CombinedLimit_macros src_HiggsAnalysis_CombinedLimit_data src_HiggsAnalysis_CombinedLimit_scripts src_HiggsAnalysis_CombinedLimit_src
ifeq ($(strip $(PyHiggsAnalysisCombinedLimit)),)
PyHiggsAnalysisCombinedLimit := self/src/HiggsAnalysis/CombinedLimit/python
src_HiggsAnalysis_CombinedLimit_python_parent := 
ALL_PYTHON_DIRS += $(patsubst src/%,%,src/HiggsAnalysis/CombinedLimit/python)
PyHiggsAnalysisCombinedLimit_files := $(patsubst src/HiggsAnalysis/CombinedLimit/python/%,%,$(wildcard $(foreach dir,src/HiggsAnalysis/CombinedLimit/python ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
PyHiggsAnalysisCombinedLimit_LOC_USE := self  
PyHiggsAnalysisCombinedLimit_PACKAGE := self/src/HiggsAnalysis/CombinedLimit/python
ALL_PRODS += PyHiggsAnalysisCombinedLimit
PyHiggsAnalysisCombinedLimit_INIT_FUNC        += $$(eval $$(call PythonProduct,PyHiggsAnalysisCombinedLimit,src/HiggsAnalysis/CombinedLimit/python,src_HiggsAnalysis_CombinedLimit_python,1,1,$(SCRAMSTORENAME_PYTHON),$(SCRAMSTORENAME_LIB),,))
else
$(eval $(call MultipleWarningMsg,PyHiggsAnalysisCombinedLimit,src/HiggsAnalysis/CombinedLimit/python))
endif
ALL_COMMONRULES += src_HiggsAnalysis_CombinedLimit_python
src_HiggsAnalysis_CombinedLimit_python_INIT_FUNC += $$(eval $$(call CommonProductRules,src_HiggsAnalysis_CombinedLimit_python,src/HiggsAnalysis/CombinedLimit/python,PYTHON))
ifeq ($(strip $(combine)),)
combine := self/src/HiggsAnalysis/CombinedLimit/bin
combine_files := $(patsubst src/HiggsAnalysis/CombinedLimit/bin/%,%,$(foreach file,combine.cpp,$(eval xfile:=$(wildcard src/HiggsAnalysis/CombinedLimit/bin/$(file)))$(if $(xfile),$(xfile),$(warning No such file exists: src/HiggsAnalysis/CombinedLimit/bin/$(file). Please fix src/HiggsAnalysis/CombinedLimit/bin/BuildFile.))))
combine_BuildFile    := $(WORKINGDIR)/cache/bf/src/HiggsAnalysis/CombinedLimit/bin/BuildFile
combine_LOC_USE := self  HiggsAnalysis/CombinedLimit boost_program_options
combine_PACKAGE := self/src/HiggsAnalysis/CombinedLimit/bin
ALL_PRODS += combine
combine_INIT_FUNC        += $$(eval $$(call Binary,combine,src/HiggsAnalysis/CombinedLimit/bin,src_HiggsAnalysis_CombinedLimit_bin,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_BIN),bin,$(SCRAMSTORENAME_LOGS)))
combine_CLASS := BINARY
else
$(eval $(call MultipleWarningMsg,combine,src/HiggsAnalysis/CombinedLimit/bin))
endif
ALL_COMMONRULES += src_HiggsAnalysis_CombinedLimit_bin
src_HiggsAnalysis_CombinedLimit_bin_parent := HiggsAnalysis/CombinedLimit
src_HiggsAnalysis_CombinedLimit_bin_INIT_FUNC += $$(eval $$(call CommonProductRules,src_HiggsAnalysis_CombinedLimit_bin,src/HiggsAnalysis/CombinedLimit/bin,BINARY))
ALL_COMMONRULES += src_HiggsAnalysis_CombinedLimit_test
src_HiggsAnalysis_CombinedLimit_test_parent := HiggsAnalysis/CombinedLimit
src_HiggsAnalysis_CombinedLimit_test_INIT_FUNC += $$(eval $$(call CommonProductRules,src_HiggsAnalysis_CombinedLimit_test,src/HiggsAnalysis/CombinedLimit/test,TEST))
src_HiggsAnalysis_CombinedLimit_scripts_files := $(filter-out \#% %\#,$(notdir $(wildcard $(foreach dir,$(LOCALTOP)/src/HiggsAnalysis/CombinedLimit/scripts,$(dir)/*))))
$(eval $(call Src2StoreCopy,src_HiggsAnalysis_CombinedLimit_scripts,src/HiggsAnalysis/CombinedLimit/scripts,$(SCRAMSTORENAME_BIN),*))
ALL_SUBSYSTEMS+=HiggsAnalysis
subdirs_src_HiggsAnalysis = src_HiggsAnalysis_CombinedLimit
