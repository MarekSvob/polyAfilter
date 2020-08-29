# COPIED function from Ding et al., 2020:
# https://bitbucket.org/jerry00/scumi-dev/src/master/R/filt_low_count_cell.R
FitLogCount = function(x, show.plot=TRUE, loglik=FALSE, filter.ratio=0.1, 
                       xlab="Expression", ylab="Density", take.log=TRUE, ...) {
  # Given an input sparse gene by cell matrix x (raw UMI count matrix or read matrix), 
  # use Student's t mixture model to remove likely low-quality cells 
  # that have small numbers of UMIs or reads 
  # (if x is the raw uncollapsed gene by cell read matrix)
  # 
  # First install the xseq R package (while inside scumi-dev/R) from terminal:
  # R CMD build xseq
  # R CMD install xseq_0.2.2.tar.gz [R CMD INSTALL xseq_0.2.2.tar.gz]
  #
  
  x.cell = Matrix::colSums(x)
  
  if (take.log==TRUE) {
    expr.quantile = log10(x.cell)
  } else {
    expr.quantile = x.cell
  }
  
  th = quantile(expr.quantile, filter.ratio)
  
  mu = c(mean(expr.quantile[expr.quantile<=th]), 
         mean(expr.quantile[expr.quantile>th]))
  
  sigma  = c(sd(expr.quantile[expr.quantile<=th]), 
             sd(expr.quantile[expr.quantile>th]))
  
  lambda = c(filter.ratio, 1 - filter.ratio)
  
  prior = list()
  prior$lambda = lambda
  prior$mu     = mu
  prior$sigma  = sigma
  
  prior$alpha = c(5, 5)
  
  prior$kapp = 0.2
  prior$dof  = 3 # max(length(x.cell) * 0.05, 2)
  
  ## 
  model = xseq:::MixStudentFitEM(expr.quantile, lambda=lambda, prior = prior, 
                                 mu=mu, sigma=sigma, K=2, nu.equal = TRUE)
  if (show.plot == TRUE) {
    xseq:::MixModelPlot(model,  xlab2=xlab, ylab2=ylab, 
                        breaks=40, loglik=loglik, ...)
  }
  
  range.x = range(expr.quantile)
  x.range = seq(range.x[1], range.x[2], diff(range.x) / 1000)
  
  expec = xseq:::MixStudentProb(x.range, K=2, model$lambda, 
                                model$mu, model$sigma, model$nu)
  
  post.prob = expec$tau
  id.cross = which(post.prob[, 1] > post.prob[, 2])
  
  id = which(diff(id.cross) > 1)
  if (length(id) >= 1) {
    id.cross = id.cross[id[1]]
  } else {
    id.cross = id.cross[length(id.cross)]
  }
  
  abline(v = x.range[id.cross], lwd=2.5, col='dodgerblue')
  
  id = which(expr.quantile >= x.range[id.cross])
  if (length(id) > 0) {
    x = x[, id]
  }
  
  x
}

# COPIED function from Ding et al., 2020:
# https://bitbucket.org/jerry00/scumi-dev/src/master/R/sc_compare_util.R
FilterByMito = function(umi.count, return.ratio=FALSE, species='human') {
  num.read = Matrix::colSums(umi.count)
  
  if (species == 'human') {
    mito.genes = grep('MT-', rownames(umi.count), value = TRUE)
  } else if (species == 'mouse') {
    mito.genes = grep('mt-', rownames(umi.count), value = TRUE)
  } else {
    stop('species must be human or mouse!')
  }
  
  mito.ratio = Matrix::colSums(umi.count[mito.genes, ]) / num.read
  mito.th = quantile(mito.ratio, 0.75) + 3 *IQR(mito.ratio)
  
  if (mito.th <= 0.02) {
    print('Mito-th \n')
    mito.th = max(max(mito.ratio) / 2, mito.th)
  }
  
  plot(log10(num.read), mito.ratio)
  abline(h=mito.th)
  
  if (TRUE == return.ratio) {
    return(list(umi.count[, mito.ratio <= mito.th], mito.ratio))
  }
  
  umi.count[, mito.ratio <= mito.th]
}

# MY SCRIPT
library(Seurat)

# Load the datafile
umi.count <- Read10X(
  data.dir = "/Volumes/rc/lab/B/BoscoG/polyA/data/Ding2020/output/Cortex1/10X/count/umi_count",
  gene.column = 1
)

# Filter cells whose posterior probability of UMI count < 0.5 from the high-quality
#  component in bimodal distribution (a mixture of two Student's t distributions)
#> "For each cortex dataset, the number of UMIs per cell barcode across all cell
#>  barcodes returned by scumi followed a bimodal distribution, with some cell
#>  barcodes having few UMIs and others having many. We therefore first used a
#>  mixture of two Student’s t distributions to fit the UMI count distribution
#>  across all of the returned cell barcodes. We considered the mixture component
#>  with a larger mean as the high-quality cell barcode component. We removed from
#>  further analyses the cell barcodes with posterior probabilities <0.5 from the
#>  high-quality component."
umi.count.c.filt <- FitLogCount(umi.count)

# Filter cells with mito ratio > 75th percentile + 3*IQR
#> "We removed cells with a high fraction of reads aligning to mitochondrial
#>  genes (names starting with ‘mt-’ for mouse and ‘MT-’ for human)— greater than
#>  75th percentile + 3 × IQR of the mitochondrial ratios across the top returned
#>  cell barcodes, where IQR stands for interquartile range."
umi.count.cm.filt <- FilterByMito(umi.count.c.filt, species = 'mouse')

#> "For each cell, its UMI counts were divided by the total number of UMIs from
#>  that cell and then scaled by multiplying by 10,000 to get transcripts per
#>  10,000 (TP10K). We then added 1 to these transcripts per 10,000 and log
#>  transformed by the natural log."
seurat.obj <- CreateSeuratObject(counts = umi.count.cm.filt)
seurat.obj <- NormalizeData(object = seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)

