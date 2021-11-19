source("helpers/kdepairs.default.R")
source("helpers/kdepairs.R")

library(tools)

cat("Did you delete all the variables you wanted? If no, run 'rm(list=ls())'")

# ----------------------------------------------------------------------
# name of the data file
#fname <- "../01/DSR01_RUN_1_CONTROL/curgen_db_005.txt"
fname <- "../CAL_KSC_FE_SET_2_RUN_1_IC_EXP/final.txt"
#fres <-"res_tmcmc_FW.tex"

# names of variables
l = c(expression(D[u]),expression(s),expression(chi),expression(D[f]),expression(r), expression(sigma ^2))
#l = c(expression(alpha),expression(phi),expression(mu[sr]),expression(mu[rs]),expression(r[1]),expression(r[2]), expression(sigma ^2))
#l = c(expression(r[s]),expression(r[r]),expression(K),expression(epsilon),expression(gamma),expression(sigma ^2))
#l = c(expression(r[s]),expression(r[r]),expression(K),expression(sigma ^2))

discrepancy <- 0  # is the last column of the data discrepancy of likelihood?
truelik <- 1  # use likelihood or log-likelihood for coloring

# ----------------------------------------------------------------------
print_tab <- function(text, x, nd, w) {
  cat(text, formatC(x, digits=nd, format='e', width=w), "\n", sep="\t")
}

print_nd <- function(x, nd) {
  cat(format(round(x, digits=nd), nsmall=nd))
}

cat("Reading file", fname,"\n")
data <- read.table(fname)
data <- array(data=unlist(data), dim=dim(data))
nd <- dim(data)[2]  # dimension of the samples

# mean
md <- colMeans(data[, 1:nd-1])
print_tab("means:\t", md, 4, 8)

# standard deviation
m2d <- colMeans(data[, 1:nd-1]^2)
sd <- sqrt(m2d-md^2)
print_tab("stds:\t", sd, 4, 8)

# most probable parameters
best_id <- which.max(data[, nd-1])
best <- data[best_id, 1:nd-1]
print_tab("MPVs:\t", best, 4, 8)

# quantiles
q1 <- c(); for (i in seq(1, nd-1)) q1[i] <- quantile(data[, i], probs=c(0.05))
print_tab("q-0.05:  ", q1, 4, 8)
q2 <- c(); for (i in seq(1, nd-1)) q2[i] <- quantile(data[, i], probs=c(0.95))
print_tab("q-0.95:  ", q2, 4, 8)

# print for latex table
a <- c(); digits <- 3
for (i in seq(1, nd-1)) {
    t <- capture.output(cat(print_nd(md[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, " & ["), collapse = "")
    t <- capture.output(cat(print_nd(q1[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, ", "), collapse = "")
    t <- capture.output(cat(print_nd(q2[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, "] & "), collapse = "")
}

#write.table(a, file=fres,sep = "\t", row.names = T)
a <- capture.output(cat(substr(a, 1, nchar(a)-3)))
a <- paste("latex:", a)
a <- paste(a, "\\\\ \n")
cat(a)



if(discrepancy) data[, nd] <-    -data[, nd]
if(truelik)     data[, nd] <- exp(data[, nd])

pname <- paste0(file_path_sans_ext(fname), ".png")
png(pname, width=2500, height=2500, units='px', res=300, pointsize=15, type="cairo", family="times")
kdepairs(data, n_1d=20, n_2d=200, labels=l) #n_2d=200
dev.off()
