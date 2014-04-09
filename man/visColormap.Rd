\name{visColormap}
\alias{visColormap}
\title{Function to define a colormap}
\usage{
visColormap(colormap = c("bwr", "jet", "gbr", "wyr", "br", "yr",
"rainbow",
"wb"))
}
\arguments{
\item{colormap}{short name for the colormap}
}
\value{
\itemize{ \item{\code{palette.name}: a function that takes
an integer argument for generating that number of colors
interpolating the given sequence} }
}
\description{
\code{visColormap} is supposed to define a colormap. It
returns a function, which will take an integer argument
specifying how many colors interpolate the given colormap.
}
\note{
The input colormap includes: \itemize{ \item{"jet": jet
colormap} \item{"bwr": blue-white-red} \item{"gbr":
green-black-red} \item{"wyr": white-yellow-red} \item{"br":
black-red} \item{"yr": yellow-red} \item{"wb": white-black}
\item{"rainbow": rainbow colormap, that is,
red-yellow-green-cyan-blue-magenta} \item{Alternatively,
any hyphen-separated HTML color names, e.g.
"blue-black-yellow", "royalblue-white-sandybrown",
"darkblue-lightblue-lightyellow-darkorange",
"darkgreen-white-darkviolet",
"darkgreen-lightgreen-lightpink-darkred". A list of
standard color names can be found in
\url{http://html-color-codes.info/color-names}} }
}
\examples{
# 1) define "blue-white-red" colormap
palette.name <- visColormap(colormap="bwr")

# 2) use the return function "palette.name" to generate 10 colors spanning "bwr"
palette.name(10)
}
\seealso{
\code{\link{visColoralpha}}
}

