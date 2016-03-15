context("test-getDistanceAsso")
notExistingKeggPathway <- "hsa:001110"
existingKeggPathway <- "hsa:01100"
emptyDF <- data.frame();
oneColDF <- data.frame("gene" = c("aa","sdd","saaa"))
twoColWrongDF <- data.frame("gene" = c("hsa:1579","dd","saaa"),
                            "metabo" = c("a","a","a"))
twoColWrongDF_1 <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                            "metabo" = c("a","a","a"))
twoColWrongDF_2 <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                            "metabo" = c("a","C19615","a"))
AssoDataGoodDF <- data.frame("gene" = as.vector(c("hsa:1579","hsa:34")),
                             "metabolite"= as.vector(c("C19615","C05271")));
MetaboGoodDF <- data.frame("metbaolite" = c("C19615","C05271"))


test_that("heatmapFunction", {


      #  arg1 <-  heatmapFunction(existingKeggPathway, emptyDF)
        expect_error(heatmapFunction(existingKeggPathway, emptyDF))
        expect_error(heatmapFunction(existingKeggPathway, oneColDF))
        expect_error(heatmapFunction(existingKeggPathway, twoColWrongDF))
        expect_error(heatmapFunction(existingKeggPathway, twoColWrongDF_1))

})
