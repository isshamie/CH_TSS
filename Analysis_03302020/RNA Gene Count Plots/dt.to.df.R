#' convert a data.table to a data.frame with unique rownames
#'
#' @param dt 
#' @param rn.col name of the column in the data table to be used as row nanes. If unspecified, the first column will be used
#'
#' @return
#' @export
#'
#' @examples
dt.to.df <- function(dt, rn.col = NULL, remove.rn = T) {
  if(is.null(rn.col)){
    rn.col <- names(dt)[1]
  }
  if (sum(names(dt) == rn.col)!=1){
    message('row name column matched to multiple columns in data.table')
    stop()
  }
  if (remove.rn){
    df <- as.data.frame(dt[, -which(names(dt)==rn.col),  with = F])
  }else{
    df <- as.data.frame(dt)
  }
  rownames(df) <- dt[, get(rn.col)]
  return(df)
}