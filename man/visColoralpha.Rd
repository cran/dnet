\name{visColoralpha}
\alias{visColoralpha}
\title{Function to add transparent (alpha) into colors}
\usage{
visColoralpha(col, alpha)
}
\arguments{
\item{col}{input colors. It can be vector of R color
specifications, such as a color name (as listed by
'colors()), a hexadecimal string of the form "#rrggbb" or
"#rrggbbaa"}

\item{alpha}{numeric vector of values in the range [0, 1]
for alpha transparency channel (0 means transparent and 1
means opaque)}
}
\value{
a vector of colors (after transparent being added)
}
\description{
\code{visColoralpha} is supposed to add transparent (alpha)
into colors.
}
\note{
none
}
\examples{
# 1) define "blue-white-red" colormap
palette.name <- visColormap(colormap="bwr")

# 2) use the return function "palette.name" to generate 10 colors spanning "bwr"
col <- palette.name(10)

# 3) add transparent (alpha=0.5)
cols <- visColoralpha(col, alpha=0.5)
}
\seealso{
\code{\link{visColormap}}
}

