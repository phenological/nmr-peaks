#' Methods for importing/exporting to/from JSON

#' Method for creating a JSON file out of an object, that contains S4 objects of the types contained in
#' this file.
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid exporting some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @export
#'
setGeneric("toJSONFile", function(obj, control=NA, con="ANY") standardGeneric("toJSONFile"))

#' Method for creating a JSON file out of an NMRPeak1D
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#' @importFrom methods is
#' @export
#'
setMethod("toJSONFile", signature(obj="NMRPeak1D", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            write("{", con, append = TRUE, sep="")
            sep <- ""
            for(slotName in names(getSlots(is(obj)))) {
              value <- slot(obj, slotName)
              if (length(value) > 0 && !all(is.na(value))) {
                write(paste0(sep, '"',slotName, '":'), con, append = TRUE, sep="")
                toJSONFile(value, control, con)
                sep <- ","
              }
            }
            write("}", con, append = TRUE, sep="")
          }
)

#' Method for creating a JSON file out of an NMRSignal1D
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @importFrom methods getSlots slot
#' @export
#'
setMethod("toJSONFile", signature(obj="NMRSignal1D", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            write("{", con, append = TRUE, sep="")
            sep <- ""
            for(slotName in names(getSlots(is(obj)))) {
              value <- slot(obj, slotName)
              if (length(value) > 0 && !all(is.na(value))) {
                write(paste0(sep, '"',slotName, '":'), con, append = TRUE, sep="")
                toJSONFile(value, control, con)
                sep <- ","
              }
            }
            write("}", con, append = TRUE, sep="")
          }
)

#' Method for creating a JSON file out of an NMRSignalModel
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @export
#'
setMethod("toJSONFile", signature(obj="NMRSignalModel", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            write("{", con, append = TRUE, sep="")
            sep <- ""
            slotNames = names(getSlots(is(obj)))
            
            # A hack to avoid the xy being exported
            if ("no_xy" %in% names(control)) {
              if (control["no_xy"])
                slotNames <- slotNames[!(slotNames %in% c('experimental', 'ppm', "fitted"))]
            }
            
            for(slotName in slotNames) {
              value <- slot(obj, slotName)
              if (length(value) > 0 && !all(is.na(value))) {
                write(paste0(sep, '"',slotName, '":'), con, append = TRUE, sep="")
                toJSONFile(value, control, con)
                sep <- ","
              }
            }
            write("}", con, append = TRUE, sep="")
          }
)

#' Method for creating a JSON file out of a list
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @export
#'
setMethod("toJSONFile", signature(obj="list", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            lnames <- names(obj)
            if (is.null(lnames)) {
              write("[", con, append = TRUE, sep="")
              sep <- ""
              for(element in obj) {
                #if(!all(is.na(element))) {
                write(sep, con, append = TRUE, sep="")
                toJSONFile(element, control, con)
                sep <- ","
                #}
              }
              write("]", con, append = TRUE, sep="")
            } else {
              write("{", con, append = TRUE, sep="")
              sep <- ""
              i <- 0
              for (slotName in lnames) {
                if (slotName == "")
                  slotName = i
                #if(!all(is.na(obj[[slotName]]))) {
                write(paste0(sep, '"',slotName, '":'), con, append = TRUE, sep="")
                toJSONFile(obj[[slotName]], control, con)
                sep <- ","
                #}
                i <- i + 1
              }
              write("}", con, append = TRUE, sep="")
            }
          }
)

#' Method for creating a JSON file out of a vector
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @export
#'
setMethod("toJSONFile", signature(obj="vector", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            write(jsonlite::toJSON(obj, control), con, append = TRUE, sep="")
          }
)

#' Method for creating a JSON file out of a numeric
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @export
#'
setMethod("toJSONFile", signature(obj="numeric", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            if (length(obj) > 1) {
              write("[", con, append = TRUE, sep="")
              sep <- ""
              for (element in obj) {
                #if(!all(is.na(element))) {
                write(sep, con, append = TRUE, sep="")
                toJSONFile(element, control, con)
                sep <- ","
                #}
              }
              write("]", con, append = TRUE, sep="")
            } else {
              if(is.na(obj)) {
                write("null", con, append = TRUE, sep="")
              } else if (is.infinite(obj)) {
                if (obj < 0) {
                  write("-2e52", con, append = TRUE, sep="")
                } else {
                  write("2e52", con, append = TRUE, sep="")
                }
              } else {
                write(as.character(obj), con, append = TRUE, sep="")
              }
            }
          }
)

#' Method for creating a JSON file out of a logical
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @export
#'
setMethod("toJSONFile", signature(obj="logical", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            if (length(obj) > 1) {
              write("[", con, append = TRUE, sep="")
              sep <- ""
              for (element in obj) {
                #if(!all(is.na(element))) {
                write(sep, con, append = TRUE, sep="")
                toJSONFile(element, control, con)
                sep <- ","
                #}
              }
              write("]", con, append = TRUE, sep="")
            } else {
              if( is.na(obj)) {
                write("null", con, append = TRUE, sep="")
              } else {
                if (obj == TRUE) {
                  write("true", con, append = TRUE, sep="")
                } else {
                  write("false", con, append = TRUE, sep="")
                }
              }
            }
          }
)

#' Method for creating a JSON file out of a character
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @export
#'
setMethod("toJSONFile", signature(obj="character", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            if (length(obj) > 1) {
              sep <- ""
              write("[", con, append = TRUE, sep="")
              for (element in obj) {
                #if(!all(is.na(element))) {
                write(sep, con, append = TRUE, sep="")
                toJSONFile(element, control, con)
                sep <- ","
                #}
              }
              write("]", con, append = TRUE, sep="")
            } else {
              write(paste0('"', obj, '"'), con, append = TRUE, sep="")
            }
          }
)         

#' Method for creating a JSON file out of a matrix
#'
#' @param obj A data object to be parsed (list, array or S4)
#' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
#' @param con A connection to the output file
#' @return void
#'
#' @export
#'
setMethod("toJSONFile", signature(obj="matrix", control="ANY", con="ANY"),
          function(obj, control=NA, con) {
            if (length(obj) > 1) {
              sep <- ""
              write("[", con, append = TRUE, sep="")
              for (i in 1:dim(obj)[[1]]) {
                #if(!all(is.na(obj[[i]]))) {
                write(sep, con, append = TRUE, sep="")
                toJSONFile(obj[i,], control, con)
                sep <- ","
                #}
              }
              write("]", con, append = TRUE, sep="")
            } else {
              write(paste0('"', obj, '"'), con, append = TRUE, sep="")
            }
          }
)

#' Create and simplify a JSON file out of a data object.
#'
#' @param data A data object to be parsed (list, array or S4)
#' @param fileName The name of the file for storing the result
#' @return void
#'
#' @export
#'
writeToJSON <- function(data, fileName) {
  file.create(fileName)
  fileConn<-file(fileName, "wb")
  toJSONFile(data, control=c(no_xy=TRUE), con=fileConn)
  close(fileConn)
  
  fileLines <-readLines(fileName, encoding="UTF-8")
  fileConn<-file(fileName,"wb")
  write(paste(fileLines, collapse = ""), fileConn, sep = "")
  close(fileConn)
}

#' Introspect a data object and transform any matching structure into the
#' corresponding S4 object
#'
#' @param input A data object to be parsed (list, array or S4)
#' @return an object
#'
#' @export
#'
setGeneric("fromVector", function(input) standardGeneric("fromVector"))

#' Introspect a data object and transform any matching structure into the
#' corresponding S4 object
#'
#' @param input A data object to be parsed (list, array or S4)
#' @return an object
#'
#' @importFrom methods getSlots slot<- new
#' @export
#'
setMethod("fromVector", signature(input="ANY"),
          function(input) {
            if (is.null(input)) return(NA)
            listNames <- names(input)
            if (is.null(listNames)) {
              if (length(input)==1 && any(c("boolean", "character", "logical", "numeric") %in%  is(input))) {
                if (input == "NA") return(NA)
                return(input)
              } else {
                tmp <- lapply(input, function(row) {fromVector(row)})
                if (length(tmp) > 0 && length(tmp[[1]]) == 1) {
                  return(unlist(tmp))
                } else {
                  return (tmp)
                }
              }
            }else if ("type" %in% listNames) {
              output <- tryCatch(new(input[["type"]])
                                 ,error=function(e){
                                   lapply(input, function(row) fromVector(row))
                                 }
              )
              if(is.object(output)){
                slotNames <- names(getSlots(is(output)))
                #if (all(!is.na(slotNames))) {
                for (slotName in slotNames) {
                  if (slotName %in% listNames) {
                    slot(output, slotName) <- fromVector(input[[slotName]])
                  }
                }
                #}  
              }
              return(output)
            } else {
              return(lapply(input, function(row) {fromVector(row)}))
            }
          }
)

# #' Method for creating a JSON file out of an Analyte
# #'
# #' @param obj A data object to be parsed (list, array or S4)
# #' @param control A set of control parameters for the transformation. We use it to avoid export some S4 slots
# #' @param con A connection to the output file
# #' @return void
# #'
# #' @export
# #'
# setMethod("toJSONFile", signature(obj="Analyte", control="ANY", con="ANY"),
#           function(obj, control=NA, con) {
#             write("{", con, append = TRUE, sep="")
#             sep <- ""
#             for(slotName in names(getSlots(is(obj)))) {
#               value <- slot(obj, slotName)
#               if (length(value) > 0 && !all(is.na(value))) {
#                 write(paste0(sep, '"',slotName, '":'), con, append = TRUE, sep="")
#                 toJSONFile(value, control, con)
#                 sep <- ","
#               }
#             }
#             write("}", con, append = TRUE, sep="")
#           }
# )