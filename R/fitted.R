#' @export
fitted.ordinalbayes<-function(object, neww = NULL, newdata, newx = NULL, model.select = "average", ...) {
  predict.ordinalbayes(object, neww, newdata, newx, model.select, ...)
}
