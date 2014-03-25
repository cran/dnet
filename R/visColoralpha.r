#' Function to add transparent (alpha) into colors
#'
#' \code{visColoralpha} is supposed to add transparent (alpha) into colors.
#'
#' @param col input colors. It can be vector of R color specifications, such as a color name (as listed by 'colors()), a hexadecimal string of the form "#rrggbb" or "#rrggbbaa"
#' @param alpha numeric vector of values in the range [0, 1] for alpha transparency channel (0 means transparent and 1 means opaque)
#' @return 
#' a vector of colors (after transparent being added)
#' @note none
#' @export
#' @seealso \code{\link{visColormap}}
#' @include visColoralpha.r
#' @examples
#' # 1) define "blue-white-red" colormap
#' palette.name <- visColormap(colormap="bwr")
#'
#' # 2) use the return function "palette.name" to generate 10 colors spanning "bwr"
#' col <- palette.name(10)
#'
#' # 3) add transparent (alpha=0.5)
#' cols <- visColoralpha(col, alpha=0.5)

visColoralpha <- function(col, alpha)
{
    
    a <- rgb2hsv(col2rgb(col))
    cols <- hsv(h=a[1,], s=a[2,], v=a[3,], alpha=alpha)

    return(cols)
}