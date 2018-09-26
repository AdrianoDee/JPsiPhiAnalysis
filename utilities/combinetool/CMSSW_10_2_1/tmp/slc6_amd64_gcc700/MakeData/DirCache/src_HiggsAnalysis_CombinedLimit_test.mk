ALL_COMMONRULES += src_HiggsAnalysis_CombinedLimit_test
src_HiggsAnalysis_CombinedLimit_test_parent := HiggsAnalysis/CombinedLimit
src_HiggsAnalysis_CombinedLimit_test_INIT_FUNC += $$(eval $$(call CommonProductRules,src_HiggsAnalysis_CombinedLimit_test,src/HiggsAnalysis/CombinedLimit/test,TEST))
