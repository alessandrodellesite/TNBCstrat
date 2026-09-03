#Analysis of icluster model tuning

# visualize model tuning results: how many clusters to choose
#BiocManager::install("iClusterPlus")
library(iClusterPlus)

# Read output files
output <- list()
dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster"

for (k in 1:5) {
  path <- file.path(dir, paste0("cv.fit.k", k, ".Rdata"))
  env <- new.env()
  load(path, envir = env)
  output[[k]] <- env$cv.fit
}

nLambda = nrow(output[[1]]$lambda) #nlambda tested: 185
nK = length(output) #nK tested: 5

#selection parameters for the model (each matrix has nlambda rows and nK cols)
BIC = getBIC(output) #Bayesian Information Criterion matrix
devR = getDevR(output) #Deviance Ratio matrix 

#for each K find lambda that minimizes BIC 
#(BIC measures the model's error penalized by its complexity: lowest BIC is the perfect trade-off where the models captures the biological signals without including too much noise)
minBICid = apply(BIC, 2, which.min)

#for each k take the deviance ratio value corresponding to the perfect lambda (lowest BIC value)
#Deviance Ratio: % of explained variation from the clusters, the higher the better the clusters separate patients
devRatMinBIC = rep(NA, nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i], i]
}


#Elbow Plot: look for the point where adding another k doesn't yield an improvement in % explained var.
#x-axis: k clusters; y-axis: % explained variation by each cluster
#pick the K after which the curve plateaus (adding an additional cluster doesn't lead to an improvement in % explained variation) 

plot(1:(nK + 1), c(0, devRatMinBIC), type = "b", pch = 19, col = "blue",
     xlab = "Number of Clusters (K + 1)",
     ylab = "% Explained Variation",
     main = "Model Selection")

#%explained variance table
var_table <- data.frame(
  Clusters = 2:(nK + 1),
  Percent_EV = devRatMinBIC * 100 #%
)
print(var_table)

#between 3 and 4 clusters there is a 3.2% increase in explained variation (not big but not even negligible). a suggestion for noisy data is to plot heatmaps (in this case for 3 and 4 clusters) to see which one shows clearer patterns.

clusters = getClusters(output)

#3 clusters
best_cluster = clusters[, 2]
best_fit = output[[2]]$fit[[minBICid[2]]]
cat("fit for 3 clusters (output[[2]]$fit[[minBICid[2]]]):", best_fit)
lamda_3 <- output[[2]]$lambda[minBICid[2],]
cat("lambda values for 3 clusters:", lambda_3)

#4 clusters
best_cluster_4 = clusters[, 3] 
best_fit_4 = output[[3]]$fit[[minBICid[3]]]
cat("fit for 4 clusters (output[[3]]$fit[[minBICid[3]]])", best_fit_4)
lamda_4 <- output[[3]]$lambda[minBICid[3],]
cat("lambda values for 3 clusters:", lambda_4)

#2 clusters
best_cluster_2 = clusters[, 1] 
best_fit_2 = output[[1]]$fit[[minBICid[1]]]
cat("fit for 2 clusters (output[[1]]$fit[[minBICid[1]]])", best_fit_2)
lamda_2 <- output[[1]]$lambda[minBICid[1],]
cat("lambda values for 2 clusters:", lambda_2)

