
get.gene.product <- function(associations_annotated){

   # print("getGeneClassification")

  briteDF <- NULL;

  for(i in 1:nrow(associations_annotated)){

      r <- strsplit(as.character(associations_annotated[i,2]), "\\s+")

      re <- NULL;
      if(is.na(r[[1]][1])){ re <- c(re, "no brite classification")

      }
      for(j in 2:length(r[[1]])){
          if(!is.na(r[[1]][j])){
              if(r[[1]][j] == "br:hsa01000"){re <-c(re, "enzyme")}
              if(r[[1]][j] == "br:hsa02000"){re <-c(re, "transporter")}
              if(r[[1]][j] == "br:hsa04131"){re <-c(re, "protein")}
          }

      }
      if(is.null(re)) re <- NA

      re <-  paste(as.vector(unlist(re)), collapse=" ")
      briteDF <- rbind(briteDF, "Classification" = data.frame(as.vector(re)))
  }

    briteDF <- data.frame("Classification" = as.vector(briteDF))
    return <- briteDF
}
