context("test-getDistanceAll")

######### ADD good outcome test

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



test_that("getDistanceAll", {


    #input an non existant pathway
    expect_error(getDistanceAll("hsa0110",testWrigthInputDF, completeMetaboDF))

    #input an empty data.frame
    expect_error(getDistanceAll("hsa01100",emptyDF, completeMetaboliteDF))

    #input an empty data.frame
    expect_error(getDistanceAll("hsa01100",testWrigthInputDF, emptyDF))

    #test completeMetaboliteDF
    expect_error(getDistanceAll("hsa01100",testWrigthInputDF, oneColDF))
    #
    #test associatedGeneMetabo values of inputs
    expect_error(getDistanceAll("hsa01100",oneColDF,completeMetaboDF))

    #test associatedGeneMetabo values of inputs
    expect_error(getDistanceAll("hsa01100",twoColWrongDF,completeMetaboDF))

    #test associatedGeneMetabo values of inputs
    expect_error(getDistanceAll("hsa01100",twoColWrongDF1,completeMetaboDF))

    #test associatedGeneMetabo values of inputs
    expect_error(getDistanceAll("hsa01100",twoColWrongDF2,completeMetaboDF))



})
