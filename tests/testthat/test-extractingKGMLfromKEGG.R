
context("test-extractingKGMLfromKEGG")

#used variables for tests
notExistingKeggPathway <- "hsa001110"
existingKeggPathway <- "hsa01100"
notExistingKeggGene <- "hsa:2716545"
existingKeggGene <- "hsa:27165"

emptyDF <- data.frame();
oneColDF <- data.frame("gene" = c("aa","sdd","saaa"))
twoColWrongDF <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                            "metabo" = c("a","C19615","a"))
AssoDataGoodDF <- data.frame("gene" = as.vector(c("hsa:1579","hsa:34")),
                            "metabolite"= as.vector(c("C19615","C05271")));
MetaboGoodDF <- data.frame("metbaolite" = c("C19615","C05271"))


test_that("getPathwayKGML", {

#    arg1 <-  getPathwayKGML(notExistingKeggPathway)
#      expect_error(arg1)

#   testArg1<-distanceGeneToAllMetabolite(notExistingKeggPathway,
#                 AssoDataGoodDF, MetaboGoodDF,existingKeggGene)
#   testArg2 <- distanceGeneToAllMetabolite(notExistingKeggPathway,
#                            emptyDF,MetaboGoodDF,existingKeggGene)
#   testArg2_1 <- distanceGeneToAllMetabolite(notExistingKeggPathway,
#                           oneColDF,MetaboGoodDF,existingKeggGene)
#   testArg2_2 <- distanceGeneToAllMetabolite(notExistingKeggPathway,
#                          twoColWrongDF,MetaboGoodDF,existingKeggGene)
#   testArg3 <- distanceGeneToAllMetabolite(notExistingKeggPathway,
#                         AssoDataGoodDF, MetaboGoodDF,existingKeggGene)
#
#   print(distanceGeneToAllMetabolite(notExistingKeggPathway,
#                                     AssoDataGoodDF, MetaboGoodDF,existingKeggGene))
#   expect_error(distanceGeneToAllMetabolite(notExistingKeggPathway,
#                         AssoDataGoodDF, MetaboGoodDF,existingKeggGene))
 # expect_equal(testArg1,1)
  ## to come testArg4
})

test_that("barplotFunctionGeneToAllMetabo",{



})
test_that("getFrequenceAssociationsDF",{})
test_that("getAssociationsDF",{})
test_that("getAssociatedMetaboByGene",{})
test_that("barplot_adjustMaximalDistance",{})
