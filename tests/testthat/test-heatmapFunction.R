context("test-heatmapFunction")

notExistingKeggPathway <- "hsa:001110"
existingKeggPathway <- "hsa:01100"

emptyDF <- data.frame();
oneColDF <- data.frame("gene" = c("aa","sdd","saaa"))

twoColWrongDF <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                            "metabo" = c("a","C19615","a"))
twoColWrongDF1 <- data.frame("gene" = c("hsa:3711","hsa:1579","hsa:34"),
                             "metabo" = c("a","C19615","a"))
twoColWrongDF2 <- data.frame("gene" = c("aa","hsa:1579","saaa"),
                             "metabo" = c("C00001","C19615","C05271"))

AssoDataGoodDF <- data.frame("gene" = as.vector(c("hsa:1579","hsa:34")),
                             "metabolite"= as.vector(c("C19615","C05271")));
MetaboGoodDF <- data.frame("metbaolite" = c("C19615","C05271"))


test_that("heatmap", {

        # test associatedGeneMetaDF values of input
        expect_error(heatmap(existingKeggPathway, emptyDF))
        expect_error(heatmap(existingKeggPathway, oneColDF))
        expect_error(heatmap(existingKeggPathway, twoColWrongDF))
        expect_error(heatmap(existingKeggPathway, twoColWrongDF1))
        expect_error(heatmap(existingKeggPathway, twoColWrongDF2))

        #test pathway
        expect_error(heatmap(notExistingKeggPathway, AssoDataGoodDF))
})
