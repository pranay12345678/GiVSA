Models <- res$sample
p <- ncol(X)
chain <- matrix(0, nrow = length(Models), ncol = p)
Model <- list()
N <- nrow(chain)
for (i in 1:length(Models)) {
    chain[i, Models[[i]]] <- 1
    Model[i] <- paste(Models[[i]], collapse = "_")
}
Model <- unlist(Model)
Model <- sort(table(Model), decreasing = TRUE)
visits <- (apply(chain, MARGIN = 2, sum))
prob <- visits / nrow(chain)
coefficients <- rep(0, p)

HPM_param <- sort(as.numeric(unlist(strsplit(names(Model)[1], "_"))))
if (length(HPM_param) == 0)  {
    HPM <- coefficients
} else {
    X_mask <- X[, HPM_param]
    HPM2 <- solve(t(X_mask) %*% X_mask) %*% (t(X_mask) %*% y)
    HPM <- coefficients
    HPM[HPM_param] <- HPM2
}

MPM_param <- which(prob >= 0.5)
if (length(MPM_param) == 0) {
    MPM <- coefficients
} else {
    X_mask <- X[, MPM_param]
    MPM2 <- solve(t(X_mask) %*% X_mask) %*% (t(X_mask) %*% y)
    MPM <- coefficients
    MPM[MPM_param] <- MPM2
}

mode_param <- res$mode
# print(mode)
if (length(mode_param) == 0) {
    mode <- coefficients
} else {
    X_mask <- X[, mode_param]
    mode2 <- solve(t(X_mask) %*% X_mask) %*% (t(X_mask) %*% y)
    mode <- coefficients
    mode[mode_param] <- mode2
}

print("The highest probability model is:")
print(HPM_param)

print("The median probability model is:")
print(MPM_param)

print("The model with highest posterior probability is:")
print(mode_param)