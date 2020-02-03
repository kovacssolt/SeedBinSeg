##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# Below you can find some code modifying the "wbs" package (version 1.3) by Baranowski and Fryzlewicz 
# which can be used to perform seeded binary segmentation (with greedy selection).

packageurl <- "http://cran.r-project.org/src/contrib/Archive/wbs/wbs_1.3.tar.gz"
install.packages(packageurl, repos = NULL, type = "source")
library(wbs)

# or
# require(devtools)
# install_version("wbs", version = "1.3", repos = "http://cran.r-project.org")
# library(wbs)

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################






##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#' seeded (deterministic) intervals
#' 
#' The function below generates seeded intervals with specified decay with endpoints in 1,2,...,n. This routine 
#' is used inside the seedBS function.
#'  
#' @details ... 
#' @param n specifying the length of the series
#' @param decay 1/a (in the notation of the paper) in the seeded interval generation 
#' @param unique.int T if unique intervals should be returned (potentially slow, hence, default is F)
#' @return a 2-column matrix with start (first column) and end (second column) points of an interval in each row
#' @examples
#' seeded.intervals(10)
##########################################################################################################################

seeded.intervals <- function(n, decay = sqrt(2), unique.int = F){
	n	<- as.integer(n)
	depth	<- log(n, base = decay)
	depth	<- ceiling(depth)
	M	<- sum(2^(1:depth)-1)
	
	boundary_mtx           <- matrix(NA, ncol = 2)
	colnames(boundary_mtx) <- c("st", "end")
	boundary_mtx[1, ]      <- c(1, n)

	depth	<- log(n, base = decay)
	depth	<- ceiling(depth)


    for(i in 2:depth){
		int_length	<- n * (1/decay)^(i-1)
		
	    	n_int		<- ceiling(round(n/int_length, 14))*2-1		# sometimes very slight numerical inaccuracies
	    	
	    	boundary_mtx	<- rbind(boundary_mtx,
			    			cbind(floor(seq(1, n-int_length, length.out = (n_int))), 
						    	ceiling(seq(int_length, n, length.out = (n_int)))))
    }

	if(unique.int){return(unique(boundary_mtx))}
	boundary_mtx
}

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#' random intervals (with a bug fixed)
#'
#' The function generates M intervals, whose endpoints are 
#' are drawn uniformly without replacements from 1,2,...,n. This routine can be
#' used inside the wbs function and is typically not called directly by the user.
#' @param n a number of endpoints to choose from
#' @param M a number of intervals to generate
#' @param unique.int T if unique intervals should be returned (potentially slow, hence, default is F)
#' @return a M by 2 matrix with start (first column) and end (second column) points of an interval in each row
#' @examples
#' random.intervals(10,100)
##########################################################################################################################

random.intervals <- function (n, M, unique.int = F){
    n <- as.integer(n)
    M <- as.integer(M)
    intervals <- matrix(0, nrow = M, ncol = 2)
    intervals[, 1] <- floor(runif(M, min = 1, max = n))		# intended starting point
    intervals[, 2] <- ceiling(runif(M, min = 1, max = n))	# intended end point
    for (i in 1:M){
    	tmp <- intervals[i, ] 
    	while(tmp[1]==tmp[2]){	                                # in case the two happen to be equal, take new ones
    		
    		tmp[1] <- floor(runif(1, min = 1, max = n)) 
    		tmp[2] <- ceiling(runif(1, min = 1, max = n))
    		intervals[i, ] <- tmp
    	}
    	if(tmp[1] > tmp[2]){intervals[i, ] <- tmp[2:1]}    	# make sure they are ordered
    }
    if(unique.int){return(unique(intervals))}
    intervals
}

# original code from wbs package (version 1.3) with a bug
# random.intervals <-	function(n,M) {
#	
#	n <- as.integer(n)
#	M <- as.integer(M)
#	intervals <- matrix(0,nrow=M,ncol=2)
#	intervals[,1] <- ceiling(runif(M)*(n-1))
#	intervals[,2] <- intervals[,1]+ ceiling(runif(M)*(n-intervals[,1]))	# !!!bug!!!
#	
#	intervals
#	
# }


##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#' A function performing seeded binary segmentation. 
#' It is slightly modifying the "wbs" function from the "wbs" package (version 1.3) to use seeded intervals.
#'
#' @param decay a decay for seeded intervals (1/a in the notation of the paper, to be precise) used in the Seeded Binary Segmentation algorithm 
#' @param ... not in use
#' @return an object of class "wbs", which contains the following fields
#' \item{x}{the input vector provided}
#' \item{n}{the length of x}
#' \item{decay}{the decay for seeded intervals}
#' \item{rand.intervals}{a logical variable indicating type of intervals}
#' \item{integrated}{a logical variable indicating use of standard binary segmentation}
#' \item{res}{a 6-column matrix with results, where 's' and 'e' denote start-
#' end points of the intervals in which change-points candidates 'cpt' have been found;
#' column 'CUSUM' contains corresponding value of CUSUM statistic; 'min.th' is the smallest 
#' threshold value for which given change-point candidate would be not added to the set of estimated 
#' change-points; the last column is the scale at which the change-point has been found} 
##########################################################################################################################

seedBS <- function (x, decay = sqrt(2), ...){
    results <- list()
    results$x <- as.numeric(x)
    results$n <- length(results$x)
    results$M <- NA
    results$integrated <- as.logical(FALSE)
    results$rand.intervals <- as.logical(FALSE)
    results$res <- matrix(nrow = 0, ncol = 6)
    if (results$n < 2) 
        stop("x should contain at least two elements")
    if (NA %in% results$x) 
        stop("x vector cannot contain NA's")
    if (var(x) == 0) 
        stop("x is a constant vector, change-point detection is not needed")
    if (decay>2) 
        stop("decay should be <= 2")
    if (decay<=1) 
        stop("decay should be > 1")
#    if (is.na(results$M)) 
#        stop("M cannot be NA")
#    if (length(results$M) > 1) 
#        stop("M should be a single integer")
#    if (results$M < 0) 
#        stop("M should be an integer > 0")
#    if (results$rand.intervals) 
#        intervals <- matrix(random.intervals(results$n, results$M), 
#            ncol = 2)
#    else {
        intervals <- matrix(seeded.intervals(results$n, decay = decay), 
            ncol = 2)
        results$M <- nrow(intervals)
#    }
#    if (results$integrated) {
#        results$res <- matrix(.C("wbs_int_rec_wrapper", x = as.double(results$x), 
#            n = as.integer(results$n), res = double(6 * (results$n - 
#                1)), intervals = as.integer(intervals), M = as.integer(results$M))$res, 
#            results$n - 1, 6)
#    }
#    else {
        results$res <- matrix(.C("wbs_rec_wrapper", x = as.double(results$x), 
            n = as.integer(results$n), res = double(6 * (results$n - 
                1)), intervals = as.integer(intervals), M = as.integer(results$M))$res, 
            results$n - 1, 6)
        results$res <- matrix(results$res[as.integer(results$res[, 
            1]) > 0, ], ncol = 6)
#    }
    colnames(results$res) <- c("s", "e", "cpt", "CUSUM", "min.th", 
        "scale")
    class(results) <- "wbs"
    results$cpt <- changepoints(results)
    return(results)
}

##########################################################################################################################

# an example
set.seed(1)
X		<- rep(rep(0:1, each = 100), 5) + rnorm(1000)
results_seedBS	<- seedBS(X)
results_wbs	<- wbs(X)
changepoints(results_seedBS)$cpt.ic[[2]]
changepoints(results_wbs)$cpt.ic[[2]]
# In this example the change point estimates are the same for seeded and wild binary segmentation, but the order
# of detection is apparently somewhat different.

