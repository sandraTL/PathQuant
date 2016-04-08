
context("test-distanceGeneToAllMetabolite")

#used variables for tests
notExistingKeggPathway <- "hsa001110"
existingKeggPathway <- "hsa01100"
notExistingKeggGene <- "hsa:2716545"
existingKeggGene <- "hsa:1579"

emptyDF <- data.frame();
oneColDF <- data.frame("gene" = c("aa","sdd","saaa"))

twoColWrongDF <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                            "metabo" = c("a","C19615","a"))
twoColWrongDF1 <- data.frame("gene" = c("hsa:3711","hsa:1579","hsa:34"),
                             "metabo" = c("a","C19615","a"))
twoColWrongDF2 <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                             "metabo" = c("C00001","C19615","C05271"))

assoDataGoodDF <- data.frame("gene" = as.vector(c("hsa:1579","hsa:34")),
                            "metabolite"= as.vector(c("C19615","C05271")));
metaboGoodDF <- data.frame("metbaolite" = c("C19615","C05271"))


test_that("distanceGeneToAllMetabolite", {

    #test pathwayID
    expect_error(distanceGeneToAllMetabolite(notExistingKeggPathway,
                 assoDataGoodDF, metaboGoodDF,existingKeggGene))

    #test associatedGeneMetaDF
    expect_error(distanceGeneToAllMetabolite(existingKeggPathway,
                 twoColWrongDF, metaboGoodDF,existingKeggGene))

    expect_error(distanceGeneToAllMetabolite(existingKeggPathway,
                 twoColWrongDF1, metaboGoodDF,existingKeggGene))
    expect_error(distanceGeneToAllMetabolite(existingKeggPathway,
                 twoColWrongDF2, metaboGoodDF,existingKeggGene))

    #test completeMetaboDF
    expect_error(distanceGeneToAllMetabolite(existingKeggPathway,
                 assoDataGoodDF, oneColDF,existingKeggGene))

    # test if gene is in associatedGeneMetaDF
    expect_error(distanceGeneToAllMetabolite(existingKeggPathway,
                assoDataGoodDF, metaboGoodDF,notExistingKeggGene))


})

test_that("barplotFunctionGeneToAllMetabo",{

})
test_that("getFrequenceAssociationsDF",{})
test_that("getAssociationsDF",{})
test_that("getAssociatedMetaboByGene",{})
test_that("barplot_adjustMaximalDistance",{})
