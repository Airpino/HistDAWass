## ---- Show overridding for Time series of TdistributionH and TMatH ----
setMethod("show",
          signature(object="TdistributionH"),
          definition=function(object){
            newdist=distributionH(object@x, object@p)
            show(newdist)
          }
)
setMethod("show",
          signature(object="TMatH"),
          definition=function(object){
            cat("to be implemented")
          }
)
setMethod("show",
          signature(object="HTS"),
          definition=function(object){
            cat("to be implemented")
          }
)

## --- Plot overloading for Timed TdistributionH  TMatH and HTS----
#' plot for a TdistributionH object
#' @name plot-TdistributionH
#' @docType methods
#' @aliases plot,TdistributionH-method
#' @description A plot function for a \code{TdistributionH} object. The function returns a representation 
#' of the histogram.
#' @param x  a \code{TdistributionH} object
#' @param type (optional) a string describing the type of plot, default="HISTO".\cr Other allowed types are 
#' \cr"CDF"=Cumulative distribution function, \cr"QF"= quantile function, \cr"DENS"=a density approximation, 
#' \cr"HBOXPLOT"=horizontal boxplot, \cr"VBOXPLOT"= vertical boxplot,
#'  @param col (optional) a string the color of the plot, default="green".
#'  @param border (optional) a string the color of the border of the plot, default="black".
#'  @export
setMethod("plot",
          signature(x = "TdistributionH"),
          function (x,  type="HISTO",col="green",border="black") 
          {
            newdist=distributionH(x@x, x@p)
            plot(newdist,  type=type,col=col,border=border)
          }
)
# setMethod("plot",
#           signature(x = "TMatH"),
#           function (x, y="missing", type="HISTO",col="green",border="black") 
#           {
#             stop("to be implemented")            
#           }
# )
# setMethod("plot",
#           signature(x = "HTS"),
#           function (x, y="missing", type="HISTO",col="green",border="black") 
#           {
#             stop("to be implemented")
#           }
# )