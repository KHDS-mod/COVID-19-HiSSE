panel.lollipop <- function(x, y,
  groups = NULL, subscripts,
  group.number = NULL, col_handle = NULL,
  alpha = 1, lty = 1, lwd = 1,
  origin = "bottom",
  scales = list(y = list(rot = 0)),
  custom_threshold = NULL, ...)
{
    panel.segments(y0 = y, x0 = 0, x1 = x, y1 = y, col = col_handle,
      alpha = alpha, lty = lty, lwd = lwd)
    #panel.segments(y0 = 0, x0 = 1.1, x1 = 1.1, y1 = 1000, col = col_handle,
    #  alpha = alpha, lty = 2, lwd = lwd)
    panel.xyplot(x, y, cex=0.4, ...)
  
}
