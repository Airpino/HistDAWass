## ---- Show overridding for Time series of TdistributionH and TMatH ----
setMethod("show",
          signature(object="TdistributionH"),
          definition=function(object){
            cat("to be implemented")
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

## --- Plot overridding for Timed TdistributionH  TMatH and HTS----
setMethod("plot",
          signature(x = "TdistributionH"),
          function (x, y="missing", type="HISTO",col="green",border="black") 
          {
            stop("to be implemented")
          }
)
setMethod("plot",
          signature(x = "TMatH"),
          function (x, y="missing", type="HISTO",col="green",border="black") 
          {
            stop("to be implemented")            
          }
)
setMethod("plot",
          signature(x = "HTS"),
          function (x, y="missing", type="HISTO",col="green",border="black") 
          {
            stop("to be implemented")
          }
)