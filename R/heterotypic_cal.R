A_cal <- function(spatial_data, decon_data, n_neigh = 1){
  meta_tmp1 <- spatial_data
  dist_df_1 <- as.matrix(dist(meta_tmp1[, c('row', 'col')], method = "euclidean"))
  dist_df_1[dist_df_1==0] <- 1
  dist_df_1[dist_df_1>2] <- 0 # keep 1, 1.141, 2 values

  meta_tmp2 <- spatial_data
  meta_tmp2$col <- 0
  dist_df_2 <- as.matrix(dist(meta_tmp2[, c('row', 'col')], method = "euclidean"))
  dist_df_2[dist_df_2==0] <- 1
  dist_df_2[dist_df_2>=2] <- 0

  neigh_spots <- dist_df_1 * dist_df_2

  #neigh_spots[neigh_spots==1] <- 0
  neigh_spots[neigh_spots>1] <- 0.5

  neigh_spots_list <- lapply(seq(n_neigh), function(n){
    neigh_spots %^% n
  })

  neigh_spots_mat <- Reduce("+", neigh_spots_list)
  neigh_spots_mat[neigh_spots_mat>0] <- 0.5
  diag(neigh_spots_mat) <- 1

  decon_mat <- t(decon_data)

  neigh_spots_celltypeB <- decon_mat %*% neigh_spots_mat

  A_list <- lapply(rownames(decon_mat), function(celltype){
    spots_celltypeA <- decon_mat[celltype, ,drop=FALSE]
    spots_celltypeA <- matrix(rep(spots_celltypeA, dim(neigh_spots_celltypeB)[1]),
                              nrow = dim(neigh_spots_celltypeB)[1],
                              byrow = TRUE)
    colnames(spots_celltypeA) <- colnames(decon_mat)

    mat <- neigh_spots_celltypeB * spots_celltypeA
    mat <- mat[-which(rownames(mat)==celltype), ]
    rownames(mat) <- paste(celltype, rownames(mat), sep = '_')
    mat <- t(mat)
    mat
  })

  A_df <- do.call(cbind, A_list) %>% as.data.frame(.)

  celltype_pairs <- combn(rownames(decon_mat), m = 2) %>% t(.) %>% as.data.frame(.)
  colnames(celltype_pairs) <- c('celltypeA', 'celltypeB')
  celltype_pairs$pairs1 <- paste(celltype_pairs$celltypeA, celltype_pairs$celltypeB, sep = '_')
  celltype_pairs$pairs2 <- paste(celltype_pairs$celltypeB, celltype_pairs$celltypeA, sep = '_')

  A_aggr_list <- lapply(seq(dim(celltype_pairs)[1]), function(i){
    pairs1 <- celltype_pairs[i, 'pairs1']
    pairs2 <- celltype_pairs[i, 'pairs2']

    mat_sum <- A_df[, pairs1, drop=FALSE] + A_df[, pairs2, drop=FALSE]
    return(mat_sum)
  })

  A_aggr_df <- do.call(cbind, A_aggr_list)

  return(A_aggr_df)
}

heter_aggr_cal <- function(A_df, n_cores = 1){

  A_obs_mat <- A_df[, grepl('heterscore_obs', colnames(A_df))]
  A_exp_mat <- A_df[, grepl('heterscore_random', colnames(A_df))]

  celltype_pairs <- sub('heterscore_random[0-9]+_', '', colnames(A_exp_mat)) %>% unique()

  A_exp_MeanSd_list <- mclapply(celltype_pairs, function(ctp){
    A_exp <- A_exp_mat[, grepl(ctp, colnames(A_exp_mat))]
    A_exp_mean <- rowMeans(A_exp)
    A_exp_sd <- apply(A_exp, 1, sd)
    A_exp_list <- list(mean=A_exp_mean, sd=A_exp_sd)
    A_exp_list
  }, mc.cores = n_cores)

  A_exp_mean_mat <- do.call(cbind, lapply(A_exp_MeanSd_list, function(x){x[[1]]}))
  colnames(A_exp_mean_mat) <- celltype_pairs

  A_exp_sd_mat <- do.call(cbind, lapply(A_exp_MeanSd_list, function(x){x[[2]]}))
  colnames(A_exp_sd_mat) <- celltype_pairs
  A_exp_sd_mat[A_exp_sd_mat==0] <- NA

  heter_df <- (A_obs_mat - A_exp_mean_mat)/A_exp_sd_mat
  colnames(heter_df) <- sub('obs_', '', colnames(heter_df))

  result <- list(A_obs = A_obs_mat, A_exp_mean = A_exp_mean_mat, A_exp_sd = A_exp_sd_mat, heterotypic = heter_df)
  result
}


#' @title Calculate the heterotypic score of spatial transcriptomics
#' @author Jiaxin Luo
#'
#' @param spatial_data data.frame; coordinates of the spots, using 'row' and 'col' (not 'imagerow' and 'imagecol')
#' @param decon_data data.frame; cell type proportions in each spot, with rows as spots and columns as cell types
#' @param n_neigh numeric; number of neighboring rings to consider
#' @param n_permutation numeric; number of permutations, default is 50
#' @param n_cores integer; number of cores for parallel computation, default is 1
#' @param seed integer; random seed
#'
#' @return dataframe with spots as rows and heterotypic scores of cell types as columns
#'
#' @export
#'
#' @import expm
#' @import tidyverse
#' @import parallel
#'
CalHeteroScore <- function(spatial_data, decon_data, n_neigh = 1, n_permutation = 50, n_cores = 1, seed = 123){
  set.seed(seed)

  if (!identical(rownames(spatial_data), rownames(decon_data))) {
    stop('Error: Row names of spatial_data and decon_data are inconsistent.')
  }

  if(!all(c('row', 'col') %in% colnames(spatial_data))){
    stop('Error: Columns "row" and "col" are required in `spatial_data`.')
  }

  if (n_neigh<=0) {
    stop('Error: n_neigh cannot be less than or equal to 0.')
  }


  # step 1: A_obs
  A_obs_df <- A_cal(spatial_data, decon_data, n_neigh)
  colnames(A_obs_df) <- paste0('heterscore_obs_', colnames(A_obs_df))

  # step 2: A_exp
  A_exp_list <- mclapply(seq(n_permutation), function(i){

    shuffle_order <- sample(rownames(spatial_data), length(rownames(spatial_data)), replace = FALSE)
    spatial_data_shuffle <- spatial_data
    rownames(spatial_data_shuffle) <- shuffle_order
    spatial_data_shuffle <- spatial_data_shuffle[rownames(spatial_data), ]

    A_exp_mat <- A_cal(spatial_data_shuffle, decon_data, n_neigh)
    A_exp_mat <- A_exp_mat[rownames(spatial_data), ]
    colnames(A_exp_mat) <- paste0('heterscore_random', i, '_', colnames(A_exp_mat))
    A_exp_mat
  }, mc.cores = n_cores)
  A_exp_df <- do.call(cbind, A_exp_list) %>% as.data.frame(.)

  A_df <- cbind(A_obs_df, A_exp_df)

  # step 3: calculate heterotypic score
  heter_list <- heter_aggr_cal(A_df, n_cores = n_cores)
  heter_df <- heter_list$heterotypic

  return(heter_df)
}

