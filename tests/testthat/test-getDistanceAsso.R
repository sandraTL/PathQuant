
context("test-getDistanceAsso")



emptyDF = data.frame();

oneColDF <- data.frame("gene" = c("aa","sdd","saaa"))
twoColWrongDF <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                            "metabo" = c("a","C19615","a"))
twoColWrongDF1 <- data.frame("gene" = c("hsa:3711","hsa:1579","hsa:34"),
                             "metabo" = c("a","C19615","a"))
twoColWrongDF2 <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                             "metabo" = c("C00001","C19615","C05271"))
testWrigthInputDF = data.frame("genes" = as.vector(c("hsa:1579","hsa:34")),
                               "metabolites"= as.vector(c("C19615","C05271")));


testWrigthOutputDF = data.frame(
                        "geneCommonName" = c("CYP4A11","ACADM"),
                        "geneKEGGId" = c("hsa:1579","hsa:34"),
                        "isGeneInMap" = c(TRUE,TRUE),
                        "metaboliteCommonName" = c("Hexadecanedioate",
                                                   "trans-Hex-2-enoyl-CoA"),
                        "metaboliteKEGGId" = c("C19615","C05271"),
                        "isMetaboliteInMap" = c(TRUE,TRUE),
                        "distance" = c(Inf,0));


test_that("getDistanceAsso", {


    #input an non existant pathway
    expect_error(getDistanceAsso("hsa0110",testWrigthInputDF, F))

    #input an empty associatedGeneMetaboDF
    expect_error(getDistanceAsso("hsa01100",emptyDF, F))

    #test expected good output
    expect_equivalent(getDistanceAsso("hsa01100",testWrigthInputDF, F),
                      testWrigthOutputDF)

    expect_equivalent(getDistanceAsso("hsa01100",testWrigthInputDF),
                      testWrigthOutputDF)

    #test associatedGeneMetaDF wrong values as Input
    expect_error(getDistanceAsso("hsa01100",twoColWrongDF, F))

    expect_error(getDistanceAsso("hsa01100",twoColWrongDF1, F))

    expect_error(getDistanceAsso("hsa01100",twoColWrongDF2, F))




})




