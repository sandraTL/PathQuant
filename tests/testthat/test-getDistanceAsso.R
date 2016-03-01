library(PADIGM)
context("test-getDistanceAsso")

#Small exemple of associated data.frame entry
# testWrigthInputDF = data.frame("genes" = as.vector(c("hsa:5105","hsa:47")),
#                                "metabolites"= as.vector(c("C00042","C00036")));

testWrigthInputDF = data.frame("genes" = as.vector(c("hsa:1579","hsa:34")),
                               "metabolites"= as.vector(c("C19615","C05271")));

#Expected output of getDistanceAsso("hsa01100", testWrightInputDF, F,
#                                                     output= "data.frame")
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

    # Normal Output
    output <- getDistanceAsso("hsa01100",testWrigthInputDF, F,"data.frame")


    #test if output is a data.frame
    expect_that( output, is_a("data.frame"))

    # test expected output values of the function
    expect_equivalent(output, testWrigthOutputDF)
    #test length of each columns and type of data in each column
    expect_equal(length(output[,1]), nrow(testWrigthInputDF))
    expect_that( output[,1], is_a("factor"))

    expect_equal(length(output[,2]), nrow(testWrigthInputDF))
    expect_that( output[,2], is_a("factor"))

    expect_equal(length(output[,3]), nrow(testWrigthInputDF))
    expect_that( output[,3], is_a("logical"))

    expect_equal(length(output[,4]), nrow(testWrigthInputDF))
    expect_that( output[,4], is_a("factor"))

    expect_equal(length(output[,5]), nrow(testWrigthInputDF))
    expect_that( output[,5], is_a("factor"))

    expect_equal(length(output[,6]), nrow(testWrigthInputDF))
    expect_that( output[,6], is_a("logical"))

    expect_equal(length(output[,7]), nrow(testWrigthInputDF))
    expect_that( output[,7], is_a("numeric"))
    #test if all distance is < 0 including Inf values
    expect_that( all(output[,7] >= 0), is_true())

    # Pathway not in database

    expect_error(getDistanceAsso("hsa0110",testWrigthInputDF, F,"data.frame")
                 ,"pathway doesn't exist in KEGG database", fixed=T)

})




