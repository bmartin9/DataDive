source('code/sim_gaussian.R')
if (!require(speedglm)) {
    install.packages("speedglm")
    library(speedglm)
}

cv_snps = readRDS('data/cv_snpids.rds')
geno.matrix = readRDS('data/geno.matrix.rds')
in.sample.LD = readRDS('data/ld_matrix.rds')
block_id = '10_30'
chr = '10'
block = '30'

cv_index = which(colnames(geno.matrix) %in% cv_snps)

# 1. Data simulation
print("Simulating data")
pve = c(0.005,0.02,0.1,0.3) # range of percentage variance explained

cvs = colnames(in.sample.LD)
cv = length(cvs)


# To simulate a new dataset for each cluster:
res = list()
for (i in 1:length(cv_index)) {
    print(paste0('Simulating data for cluster ', i))
    cluster = i
    pve.val = sample(pve,1)
    res[[i]] = sim_gaussian(geno.matrix, pve.val, effect_num = 1, cv_idx = cv_index[i], cv_cor = in.sample.LD)
}

# To simulate a new null cluster dataset
null_sim = sim_gaussian(geno.matrix, pve.val, effect_num = 0, cv_idx = NULL, cv_cor = in.sample.LD)
