
assoDataXlsx <- function(distanceAssoDataDF){

# exporting data.frame to excel is easy with xlsx package

sheetname <- "Distance"

xlsx::write.xlsx(distanceAssoDataDF, "assoDistance.xlsx",
                 row.names=FALSE,sheetName=sheetname,showNA=TRUE)
file <- "assoDistance.xlsx"
# but we want to highlight cells if value greater than or equal to 5
wb <-xlsx::loadWorkbook(file)# load workbook
#wb <- xlsx::createWorkbook()# load workbook

# print(wb)
#
fo1 <- xlsx::Fill(foregroundColor="#225EA8")
 print("1")
 print(class(fo1))
 print("2")
 print(class(wb))# create fill object
cs1 <- xlsx::CellStyle(wb, fo1)           # create cell style DARK BLUE
fo2 <- xlsx::Fill(foregroundColor="#7FCDBB")   # LIGHT BLUE
cs2 <- xlsx::CellStyle(wb, fill=fo2)
fo3 <- xlsx::Fill(foregroundColor="#CB181D")    # DARK RED
cs3 <- xlsx::CellStyle(wb, fill=fo3)
fo4 <- xlsx::Fill(foregroundColor="#FC9272")    # LIGHT RED
cs4 <- xlsx::CellStyle(wb, fill=fo4)
fo5 <- xlsx::Fill(foregroundColor="#238443")    # DARK GREEN
cs5 <- xlsx::CellStyle(wb, fill=fo5)
fo6 <- xlsx::Fill(foregroundColor="#ADDD8E")    # LIGHT GREEN
cs6 <- xlsx::CellStyle(wb, fill=fo6)
#
# #
 sheets <- xlsx::getSheets(wb)
 print(sheets)# get all sheets
 sheet <- xlsx::sheets[[sheetname]]          # get specific sheet
 rows <- xlsx::getRows(sheet, rowIndex=1:(nrow(distanceAssoDataDF)+1))# get rows
# headerRows <- xlsx::getRows(sheet, rowIndex=1)
#
#
#Color selected cells
getCellsToHighligh(1,rows,cs2)
getCellsToHighligh(2,rows,cs2)
getCellsToHighligh(3,rows,cs2)
getCellsToHighligh(1,headerRows,cs1)
getCellsToHighligh(2,headerRows,cs1)
getCellsToHighligh(3,headerRows,cs1)
getCellsToHighligh(4,rows,cs4)
getCellsToHighligh(5,rows,cs4)
getCellsToHighligh(6,rows,cs4)
getCellsToHighligh(4,headerRows,cs3)
getCellsToHighligh(5,headerRows,cs3)
getCellsToHighligh(6,headerRows,cs3)
getCellsToHighligh(7,rows,cs6)
getCellsToHighligh(7,headerRows,cs5)


cells <- xlsx::getCells(rows, colIndex=1:8)
values <- lapply(cells, xlsx::getCellValue) # extract the cell values

xlsx::saveWorkbook(wb, file)
}


getCellsToHighligh <- function(col, rows, colour){

    cells <- xlsx::getCells(rows, colIndex=col)         # get cells
    valuesGene <- lapply(cells, xlsx::getCellValue)
    highlightGene <- as.numeric(as.vector(names(valuesGene)))

    lapply(names(cells[highlightGene]),
           function(ii)xlsx::setCellStyle(cells[[ii]],colour))


}
