k_cal <- function(spatial_data, decon_data, n_neigh = 1){

  meta_tmp1 <- spatial_data
  dist_df_1 <- as.matrix(dist(meta_tmp1[, c('row', 'col')], method = "euclidean"))
  dist_df_1[dist_df_1==0] <- 1
  dist_df_1[dist_df_1>2] <- 0 # keep 1, 1.141, 2 values

  meta_tmp2 <- spatial_data
  meta_tmp2$col <- 0
  dist_df_2 <- as.matrix(dist(meta_tmp2[, c('row', 'col')], method = "euclidean"))
  dist_df_2[dist_df_2==0] <- 1
  dist_df_2[dist_df_2>=2] <- 0 # only keep 1

  neigh_spots <- dist_df_1 * dist_df_2

  neigh_spots[neigh_spots==1] <- 0
  neigh_spots[neigh_spots>0] <- 1

  neigh_spots_list <- lapply(seq(n_neigh), function(n){
    neigh_spots %^% n
  })

  neigh_spots_mat <- Reduce("+", neigh_spots_list)
  diag(neigh_spots_mat) <- 0
  neigh_spots_mat[neigh_spots_mat>0] <- 1

  decon_mat <- t(decon_data)

  neigh_spots_celltype <- decon_mat %*% neigh_spots_mat

  result <- neigh_spots_celltype * decon_mat
  result <- t(result)
  return(result)
}


homo_aggr_cal <- function(decon_data, k_df){
  k_obs_mat <- k_df[, grepl('obs', colnames(k_df))]
  k_exp_mat <- k_df[, grepl('random', colnames(k_df))]

  celltypes <- sub('homoscore_obs_', '', colnames(k_obs_mat))

  # mean k_exp value
  k_exp_mean_mat <- lapply(celltypes, function(celltype){
    exp_mat <- k_exp_mat[, grepl(celltype, colnames(k_exp_mat))]
    mean_exp_mat <- rowMeans(exp_mat)
    mean_exp_mat <- data.frame(celltype = mean_exp_mat, row.names = rownames(exp_mat))
    colnames(mean_exp_mat) <- celltype
    mean_exp_mat
  })%>% do.call(cbind, .) %>% as.data.frame()

  # k_obs - k_exp_mean
  homoscore_mat_1 <- as.matrix(k_obs_mat)-as.matrix(k_exp_mean_mat)

  # sum celltype proportion
  celltype_prob <- colSums(decon_data)
  celltype_prob_mat <- matrix(rep(celltype_prob,dim(homoscore_mat_1)[1]), byrow = TRUE,
                              nrow = dim(homoscore_mat_1)[1],
                              ncol = length(celltype_prob))
  rownames(celltype_prob_mat) <- rownames(homoscore_mat_1)
  colnames(celltype_prob_mat) <- colnames(decon_data)


  homoscore_mat_2 <- homoscore_mat_1/celltype_prob_mat
  return(homoscore_mat_2)
}

#' @title Calculate the homotypic score of spatial transcriptomics
#' @author Jiaxin Luo
#'
#' @param spatial_data data.frame; coordinates of the spots, using 'row' and 'col' (not 'imagerow' and 'imagecol')
#' @param decon_data data.frame; cell type proportions in each spot, with rows as spots and columns as cell types
#' @param n_neigh numeric; number of neighboring rings to consider
#' @param n_permutation numeric; number of permutations, default is 50
#' @param n_cores integer; number of cores for parallel computation, default is 1
#' @param seed integer; random seed
#'
#' @return dataframe with spots as rows and homotypic scores of cell types as columns
#'
#' @export
#'
#' @import expm
#' @import tidyverse
#' @import parallel
#'
CalHomoScore <- function(spatial_data, decon_data, n_neigh = 1, n_permutation = 50, n_cores = 1, seed =123){

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

  # step 1: k_obs
  k_obs_df <- k_cal(spatial_data, decon_data, n_neigh = n_neigh)
  colnames(k_obs_df) <- paste0('homoscore_obs_', colnames(k_obs_df))

  # step 2: k_exp
  k_exp_list <- mclapply(seq(n_permutation), function(i){

    shuffle_order <- sample(rownames(spatial_data), length(rownames(spatial_data)), replace = FALSE)
    spatial_data_shuffle <- spatial_data
    rownames(spatial_data_shuffle) <- shuffle_order
    spatial_data_shuffle <- spatial_data_shuffle[rownames(spatial_data), ]


    k_exp_mat <- k_cal(spatial_data_shuffle, decon_data)
    colnames(k_exp_mat) <- paste0('homoscore_random', i, '_', colnames(k_exp_mat))

    return(k_exp_mat)

  }, mc.cores = n_cores)

  k_exp_df <- do.call(cbind, k_exp_list) %>% as.data.frame(.)

  k_df <- cbind(k_obs_df, k_exp_df)

  # step 3: calculate homotypic score
  homo_df <- homo_aggr_cal(decon_data, k_df)
  colnames(homo_df) <- sub('obs_', '', colnames(homo_df))

  return(homo_df)
}



