library("corpcor")
library(Rcpp)
library("lme4")
library(parallel)

power <- function(k, n) {
    return(log(k) / log(n))
}
Givsa <- function(y, X, M = NULL, N = 1e3, bknot = NULL, start = NULL, groups = NULL, min_group_size = 2, target_size = 0, delta = 1e-2, r = NULL, c0 = 1e-2, c = 0.8, al = 2 / 3) {
    # Initializing parameters
    if (is.null(bknot)) {
        bknot <- N %/% 4
    }

    if (is.null(M)) {
        M <- Inf
    }

    if (is.null(start)) {
        start <- numeric()
    }

    p <- ncol(X)
    n <- nrow(X)

    if (is.null(groups)) {
        groups <- rep(1, p)
    }
    var2group <- groups

    # The variables t, r, b, delta, c0 are not changed, only c is changed according to th desired effect
    t <- power(p, n)
    b <- 1 / 3

    if (is.null(r)) {
        r <- 1e-1 + max(b - (t * delta / 2), 1e-1)
    }

    g <- c0 * (p^(2 * (1 + delta)))
    k <- c * (n^r)

    # Pre-Processing
    YTX <- t(y) %*% X
    norm_y <- as.numeric(t(y) %*% y)

    # Clubbing independent variables in a single cluster
    group_sizes <- numeric(max(var2group))
    for (i in 1:length(group_sizes)) {
        group_sizes[i] <- length(which(var2group == i))
    }
    count <- 1
    for (i in 1:length(group_sizes)) {
        if (group_sizes[i] >= min_group_size) {
            count <- count + 1
        }
    }
    count2 <- 1
    for (i in 1:length(group_sizes)) {
        if (group_sizes[i] >= min_group_size) {
            var2group[var2group == i] <- count2
            count2 <- count2 + 1
        } else {
            var2group[var2group == i] <- count
        }
    }

    # Size of each group
    group_sizes <- table(var2group)

    # No of groups
    group_count <- length(group_sizes)

    # Initializing parmeters for each group
    count <- rep(0, group_count)
    Groups <- 1:group_count
    group_freq <- rep(1, group_count)
    group_elements <- mclapply(Groups, function(x) {
        which(var2group == x)
    })
    group_potential <- rep(0, group_count)

    group_mu <- numeric(group_count)
    for (i in Groups) {
        y_g <- X[, group_elements[[i]]] %*% (pseudoinverse(t(X[, group_elements[[i]]]) %*% X[, group_elements[[i]]] + 1e-3 * diag(1, length(group_elements[[i]]))) %*% (YTX[1, group_elements[[i]]]))
        group_mu[i] <- max(1e-20, as.numeric(t(y_g) %*% y_g))
    }
    group_potential <- group_mu / max(group_mu)
    mode <- numeric()
    mode_val <- 0

    # Helper functions

    # Probability of choosing a move given the size of the current model
    pdf <- function(x, move) {
        if (move == 1) {
            if (x == 0) {
                return(1)
            }
            if (x >= n) {
                return(0)
            }
            return(1 / 3)
        }
        if (move == 0) {
            if (x == 0) {
                return(0)
            }
            if (x >= n) {
                return(1 / 2)
            }
            return(1 / 3)
        }
        if (move == -1) {
            if (x == 0) {
                return(0)
            }
            if (x >= n) {
                return(1 / 2)
            }
            return(1 / 3)
        }
    }

    # An old function, currently it works as follows:
    # Action refer to either add move or remove move
    # if action = add, if the variate is already in the model return 0 else return 1
    # if action = remove, if the variate is in the model return 1 else return 0

    w <- function(var, block_size, action, divider = 1) {
        var <- as.numeric(var)
        f <- function(block_size, divider) {
            if ((block_size * divider) > M) {
                return(M / (block_size * divider))
            }
            return(1)
        }
        g <- function(block_size) {
            return(1)
        }
        prob <- 0
        if (action == 1) {
            prob <- (1 - var) * f(block_size, divider / 3)
        }
        if (action == -1) {
            prob <- var * g(block_size)
        }
        prob <- min(1, max(0, prob))
        return(c(prob, 1 - prob))
    }

    ## The posterior function returns density of the model upto a proportionality constant
    posterior <- function(model) {
        cardinality <- length(model)
        R_squared <- 0
        prob <- (1 + g)^(-cardinality / 2)
        prob <- prob * (k^(-abs(cardinality - target_size)))
        if (cardinality != 0) {
            xtx <- t(X[, model]) %*% X[, model]
            tryCatch(
                {
                    Px <- solve(xtx)
                    ytx <- YTX[, model]
                    R_squared <- as.numeric(t(ytx) %*% Px %*% ytx) / norm_y
                },
                error = function(e) {
                    prob <- 0
                }
            )
        }
        prob <- max(0, prob * ((1 - (g / (1 + g)) * R_squared)^(-n / 2)) / choose(p, cardinality))
        if (prob > mode_val) {
            mode <<- model
            mode_val <<- prob
            group_freq <<- rep(0, group_count)
            group_freq[unique(var2group[mode])] <<- 1
        }
        return(prob)
    }

    # Functions for executing the moves
    add <- function(i, block_size, gamma) {
        foo <- NULL
        eta <- sample(c(1, 0), 1, prob = w(0, block_size, 1))
        if (eta == 1) {
            foo <- c(gamma, i)
        }
        return(foo)
    }
    remove <- function(i, block_size, gamma) {
        foo <- NULL
        eta <- sample(c(1, 0), 1, prob = w(1, block_size, -1))
        if (eta == 1) {
            foo <- gamma[gamma != i]
        }
        return(foo)
    }
    swap <- function(j, i, block_size, gamma) {
        foo <- NULL
        card <- length(gamma)
        eta_a <- sample(c(1, 0), 1, prob = w(0, block_size, 1, divider = card))
        if (eta_a == 1) {
            foo <- c(gamma[gamma != i], j)
        }
        return(foo)
    }

    # An old function no longer in use
    # adapt <- function(size) {
    #     size[size <= M] <- min(size)
    #     size <- size / max(size)
    #     return(size)
    # }
    sourceCpp("./Hash.cpp", env = environment())

    MOVES <- c(-1, 0, 1)
    run <- function(gamma, t) {
        # Setting up values for the current run

        # cardinality of the current model
        card <- length(gamma)

        # pick the current move
        move <- sample(MOVES, 1, prob = c(pdf(card, -1), pdf(card, 0), pdf(card, 1)))

        # Code to select a group
        prob <- rep(0, group_count)
        prob <- group_potential
        prob <- prob / sum(prob)
        prob <- al * prob + (1 - al) * (1 / group_count)
        prob <- prob / sum(prob)
        current_block <- sample(Groups, 1, prob = prob)
        block_size <- group_sizes[[current_block]]

        N_f <- list()
        # Code for forward addition neighbourhood
        if (move == 1) {
            N_f <- mclapply(setdiff(group_elements[[current_block]], gamma), add, block_size = block_size, gamma = gamma)
            N_f <- Filter(length, N_f)
        }
        # Code for forward removal neighbourhood
        if (move == -1) {
            N_f <- mclapply(gamma, remove, block_size = block_size, gamma = gamma)
        }
        # Code for forward swap neighbourhood
        if (move == 0) {
            for (i in gamma) {
                eta_r <- sample(c(1, 0), 1, prob = w(1, block_size, -1))
                if (eta_r == 1) {
                    temp <- mclapply(setdiff(group_elements[[current_block]], gamma), swap, i = i, block_size = block_size, gamma = gamma)
                    temp <- Filter(length, temp)
                    N_f <- append(temp, N_f)
                }
            }
        }
        if (length(N_f) == 0) {
            # Taking care of problematic corner cases
            print("Empty class")
        }
        fwd_scores <- unlist(mclapply(N_f, FUN = posterior))
        fwd_scores[is.na(fwd_scores)] <- 0

        # Choosing the proposal based on forward scores
        if (all(fwd_scores == 0)) {
            proposal <- sample(N_f, 1)[[1]]
        } else {
            proposal <- (sample(N_f, 1, prob = fwd_scores))[[1]]
        }
        N_b <- list()

        # Creating backward neighbourhoods
        r <- 0
        removed_block <- 0
        rblock_size <- 0
        # backward addition neighbourhood
        if (move == -1) {
            # r is the removed covariate
            r <- setdiff(gamma, proposal)
            removed_block <- var2group[r]
            rblock_size <- group_sizes[removed_block]
            N_b <- mclapply(setdiff(group_elements[[removed_block]], proposal), add, block_size = rblock_size, gamma = proposal)
            N_b <- Filter(length, N_b)
        }
        # backward removal neighbourhood
        else if (move == 1) {
            N_b <- mclapply(proposal, remove, block_size = block_size, gamma = proposal)
        }
        # backward swap neighbourhood
        else if (move == 0) {
            r <- setdiff(gamma, proposal)
            removed_block <- var2group[r]
            rblock_size <- group_sizes[var2group[r]]
            for (i in proposal) {
                eta_r <- sample(c(1, 0), 1, prob = w(1, rblock_size, -1))
                if (eta_r == 1) {
                    temp <- mclapply(setdiff(group_elements[[removed_block]], proposal), swap, i = i, block_size = rblock_size, gamma = proposal)
                    temp <- Filter(length, temp)
                    N_b <- append(N_b, temp)
                }
            }
        }

        if (!(list(gamma) %in% N_b)) {
            N_b <- append(N_b, list(gamma))
        }
        bwd_scores <- unlist(mclapply(N_b, FUN = posterior))

        # Calculating probability of acceptance
        alpha <- 1
        if (move == 0) {
            alpha <- (prob[removed_block] * w(0, rblock_size, 1, card)[1] * w(1, rblock_size, -1)[1] * sum(fwd_scores)) / (prob[current_block] * w(0, block_size, 1, card)[1] * w(1, block_size, -1)[1] * sum(bwd_scores))
        }
        if (move == -1) {
            alpha <- (pdf(card - 1, 1) * prob[removed_block] * w(0, rblock_size, 1)[1] * sum(fwd_scores)) / (pdf(card, -1) * w(1, block_size, -1)[1] * sum(bwd_scores))
        }
        if (move == 1) {
            alpha <- (pdf(card + 1, -1) * w(1, block_size, -1)[1] * sum(fwd_scores)) / (prob[current_block] * pdf(card, 1) * w(0, block_size, 1)[1] * sum(bwd_scores))
        }
        alpha <- min(1, alpha)

        # Adding the new sample to the chain
        u <- runif(1)
        if (u <= alpha) {
            gamma <- proposal
            # accepted <- TRUE
        }
        return(gamma)
    }
    # The for loop was pushed to C++ for faster computation
    samples <- list()
    samples <- mcmc(start, run, N, bknot)

    # res contains all the models visited by the chain and the model with highest posterior probability
    res <- list()
    res$samples <- samples
    res$mode <- mode
    return(res)
}
