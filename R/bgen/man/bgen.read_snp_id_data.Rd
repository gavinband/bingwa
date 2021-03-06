\name{bgen.read_snp_id_data}
\alias{bgen.read_snp_id_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bgen.read_snp_id_data(header)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{
%%     ~~Describe \code{header} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (header) 
{
    require(bitops)
    con = header$connection
    endian = header$endian
    data = list()
    data$number_of_samples = readBin(con, integer(), size = 4, 
        n = 1, endian = endian)
    if (length(data$number_of_samples) == 0) {
        return(NULL)
    }
    id_storage_size = readBin(con, integer(), size = 1, n = 1, 
        endian = endian, signed = FALSE)
    SNPID_size = readBin(con, integer(), size = 1, n = 1, endian = endian, 
        signed = FALSE)
    data$SNPID = readChar(con, nchars = SNPID_size)
    if (id_storage_size > SNPID_size) {
        readChar(con, nchar = id_storage_size - SNPID_size)
    }
    rsid_size = readBin(con, integer(), size = 1, n = 1, endian = endian, 
        signed = FALSE)
    data$rsid = readChar(con, nchars = rsid_size)
    if (id_storage_size > rsid_size) {
        readChar(con, nchar = id_storage_size - rsid_size)
    }
    chromosome = readBin(con, integer(), size = 1, n = 1, endian = endian, 
        signed = FALSE)
    if (chromosome < 23) {
        data$chromosome = sprintf("\%02d", chromosome)
    }
    else if (chromosome == 23) {
        data$chromosome = "0X"
    }
    else if (chromosome == 24) {
        data$chromosome = "0Y"
    }
    else if (chromosome == 253) {
        data$chromosome = "XY"
    }
    else if (chromosome == 254) {
        data$chromosome = "MT"
    }
    else if (chromosome == 255) {
        data$chromosome = NA
    }
    else {
        stop(sprintf("Unknown chromosome \%d", chromosome))
    }
    data$position = readBin(con, integer(), size = 4, n = 1, 
        endian = endian)
    if (bitAnd(header$flags, 2)) {
        alleleA_size = readBin(con, integer(), size = 1, n = 1, 
            endian = endian, signed = FALSE)
        data$alleleA = readChar(con, nchars = alleleA_size)
        alleleB_size = readBin(con, integer(), size = 1, n = 1, 
            endian = endian, signed = FALSE)
        data$alleleB = readChar(con, nchars = alleleB_size)
    }
    else {
        data$alleleA = readChar(con, nchars = 1)
        data$alleleB = readChar(con, nchars = 1)
    }
    return(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
