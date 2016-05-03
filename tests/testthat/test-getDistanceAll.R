context("test-getDistanceAll")

######### ADD good outcome test

emptyDF = data.frame();

oneColDF <- data.frame("gene" = c("aa","sdd","saaa"))
twoColWrongDFGene <- data.frame("gene" = c("aa","hsa:1579","saaa"))
twoColWrongDFMetabo <- data.frame("metabo" = c("a","C19615","a"))
testWrigthInputDFGene <- data.frame("genes" = as.vector(c("hsa:1579","hsa:34")))
testWrigthInputDFMetabo <- data.frame("metabolites"= as.vector(c("C19615","C05271")))

test_that("getDistanceAll", {

    #input an non existant pathway
    expect_error(getDistanceAll("hsa0110",testWrigthInputDFGene, completeMetaboDF))

    #input an empty data.frame
    expect_error(getDistanceAll("hsa01100",emptyDF, completeMetaboDF))

    #input an empty data.frame
    expect_error(getDistanceAll("hsa01100",testWrigthInputDFGene, emptyDF))

    #test completeMetaboliteDF
    expect_error(getDistanceAll("hsa01100",testWrigthInputDFGene, oneColDF))
    #
    #test associatedGeneMetabo values of inputs
    expect_error(getDistanceAll("hsa01100",oneColDF,completeMetaboDF))

    #test associatedGeneMetabo values of inputs
    expect_error(getDistanceAll("hsa01100",twoColWrongDFGene,completeMetaboDF))


})
