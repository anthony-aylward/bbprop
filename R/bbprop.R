#===============================================================================
# bbprop.R
#===============================================================================

#' @importFrom VGAM dbetabinom


# Functions ====================================================================

#' @title Region of integration for pbbprop or bbprop_test
#' 
#' @param x double. value of difference of proportions
#' @param size1,size2 integer. number of trials for each sample
#' @param func characeter. function to determine region for. should be "cdf" or
#'     "test".
#' @param exhaustive logical. if true check in exaustive mode. false by default.
#' @return matrix, the coordinates comprising the region.
region <- function(
    x,
    size1,
    size2,
    func = "cdf",
    exhaustive = FALSE
) {
    region_indicators <- matrix(FALSE, nrow = size1, ncol = size2)
    if (!exhaustive) {
        if (func == "cdf") {
            for (r1 in seq_len(size1)) {
                if (r1 < x * size1) {
                    r2_range <- seq_len(floor(size2/size1 * (r1 + x * size1)))
                } else if (r1 > (1 - x) * size1) {
                    r2_range <- seq(ceiling(size2/size1 * (r1 - x * size1)), size2)
                } else {
                    r2_range <- seq(
                        ceiling(size2/size1 * (r1 - x * size1)),
                        floor(size2/size1 * (r1 + x * size1))
                    )
                }
                for (r2 in r2_range) region_indicators[r1, r2] <- TRUE
            }
        } else if (func == "test") {
            for (r1 in seq_len(size1)) {
                if (r1 < x * size1) {
                    r2_range <- seq(ceiling(size2/size1 * (r1 + x * size1)), size2)
                } else if (r1 > (1 - x) * size1) {
                    r2_range <- seq_len(floor(size2/size1 * (r1 - x * size1)))
                } else {
                    r2_range <- c(
                        seq_len(floor(size2/size1 * (r1 - x * size1))),
                        seq(ceiling(size2/size1 * (r1 + x * size1)), size2)
                    )
                }
                for (r2 in r2_range) region_indicators[r1, r2] <- TRUE
            }
        } else {
            stop('invalid "func" option')
        }
    } else {
        if (func == "cdf") {
            for (r1 in seq_len(size1)) {
                for (r2 in seq_len(size2)) {
                     if (abs(r1 * size2 - r2 * size1) <= size1 * size2 * x) {
                         region_indicators[r1, r2] <- TRUE
                     }
                }
            }
        } else if (func == "test") {
            for (r1 in seq_len(size1)) {
                for (r2 in seq_len(size2)) {
                     if (abs(r1 * size2 - r2 * size1) >= size1 * size2 * x) {
                         region_indicators[r1, r2] <- TRUE
                     }
                }
            }
        } else {
            stop('invalid "func" option')
        }
    }
    region_indicators
}

count_probability <- function(
    indicator,
    r1,
    r2,
    size,
    prob = rep(0.5, 2),
    rho = rep(0, 2),
    shape1 = NULL,
    shape2 = NULL,
    ...
) {
    if (!indicator) {
        return(0)
    } else {
        r <- c(r1, r2)
        prod(
            vapply(
                seq_len(2),
                function(index){
                    if (is.null(shape1) && is.null(shape2)) {
                        dbetabinom(
                            r[[index]],
                            size[[index]],
                            prob = prob[[index]],
                            rho = rho[[index]],
                            ...
                        )
                    } else if (is.numeric(shape1) && is.numeric(shape2)) {
                        dbetabinom.ab(
                            x[[index]],
                            size[[index]],
                            shape1 = shape1[[index]],
                            shape2 = shape2[[index]],
                            ...
                        )
                    }
                },
                double(length=1)
            )
        )
    }
}

#' @title CDF for difference of beta-binomial proportions
#' 
#' @param x,q vector of quantiles
#' @param size vector giving the number of trials
#' @param prob vector giving probability of success
#' @param rho vector giving correlation parameter
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'     beta distribution. See the documentation for \code{Betabinom} in the
#'     \code{VGAM} package.
#' @param lower_tail logical. If TRUE (default), probabilities are
#'     \eqn{P[X <= x]} otherwise, \eqn{P[X > x]}.
#' @param exhaustive logical. if true check in exaustive mode. false by default.
#' @return numeric, the value of the cumulative distribution function
#' @export
pbbprop <- function(
    q,
    size,
    shape1 = 1,
    shape2 = 1,
    lower_tail = TRUE,
    exhaustive = FALSE
) {
    if (length(size) == 1) size <- rep(size, 2)
    if (length(shape1) == 1) shape1 <- rep(shape1, 2)
    if (length(shape2) == 1) shape1 <- rep(shape2, 2)

    region_indicators <- region(
        q, size[[1]], size[[2]], exhaustive = exhaustive
    )

    X <- cbind(
        as.integer(region_indicators),
        matrix(
            unlist(
                lapply(
                    seq_len(size[[1]]),
                    function(x) lapply(seq_len(size[[2]]), function(y) c(x,y))
                )
            ), 
            ncol = 2,
            byrow = T
        )
    )
    count_probability_partial <- function(indicator, r1, r2) {
        count_probability(
            indicator, r1, r2, size = size, shape1 = shape2, shape2 = shape2
        )
    }
    apply(X, 1, count_probability_partial)
}