PPP_sleuth_prep <- function(sample_to_covariates, full_model = NULL, filter_fun = basic_filter, 
    target_mapping = NULL, max_bootstrap = NULL, norm_fun_counts = norm_factors, 
    norm_fun_tpm = norm_factors, aggregation_column = NULL, read_bootstrap_tpm = FALSE, 
    extra_bootstrap_summary = FALSE, transformation_function = log_transform, 
    num_cores = max(1L, parallel::detectCores() - 1L), ...) 
{
    if (!is(sample_to_covariates, "data.frame")) {
        stop(paste0("'", substitute(sample_to_covariates), "' (sample_to_covariates) must be a data.frame"))
    }
    if (!is(full_model, "formula") && !is(full_model, "matrix") && 
        !is.null(full_model)) {
        stop(paste0("'", substitute(full_model), "' (full_model) must be a formula or a matrix"))
    }
    if (!("sample" %in% colnames(sample_to_covariates))) {
        stop(paste0("'", substitute(sample_to_covariates), "' (sample_to_covariates) must contain a column named 'sample'"))
    }
    if (!("path" %in% colnames(sample_to_covariates))) {
        stop(paste0("'", substitute(sample_to_covariates)), "' (sample_to_covariates) must contain a column named 'path'")
    }
    if (!is.null(target_mapping) && !is(target_mapping, "data.frame")) {
        stop(paste0("'", substitute(target_mapping), "' (target_mapping) must be a data.frame or NULL"))
    }
    else if (is(target_mapping, "data.frame")) {
        if (!("target_id" %in% colnames(target_mapping))) {
            stop(paste0("'", substitute(target_mapping), "' (target_mapping) must contain a column named 'target_id'"))
        }
    }
    if (!is.null(max_bootstrap) && max_bootstrap <= 0) {
        stop("max_bootstrap must be > 0")
    }
    if (any(is.na(sample_to_covariates))) {
        warning("Your 'sample_to_covariance' data.frame contains NA values. This will likely cause issues later.")
    }
    if (is(full_model, "matrix") && nrow(full_model) != nrow(sample_to_covariates)) {
        stop("The design matrix number of rows are not equal to the number of rows in the sample_to_covariates argument.")
    }
    if (!is(norm_fun_counts, "function")) {
        stop("norm_fun_counts must be a function")
    }
    if (!is(norm_fun_tpm, "function")) {
        stop("norm_fun_tpm must be a function")
    }
    if (!is.null(aggregation_column) && is.null(target_mapping)) {
        stop(paste("You provided a 'aggregation_column' to aggregate by,", 
            "but not a 'target_mapping'. Please provided a 'target_mapping'."))
    }
    num_cores <- sleuth:::check_num_cores(num_cores)
    sleuth:::msg("reading in kallisto results")
    PPP_stc_sample <<- sample_to_covariates$sample <- as.character(sample_to_covariates$sample)
    if (nrow(sample_to_covariates) == 1 && !is.null(full_model)) {
        warning("There is only one sample present, but you also provided a model. ", 
            "The model will be set to NULL to prevent downstream errors.\n", 
            "The sample can be viewed using sleuth_live after preparation, ", 
            "but you need more than one sample to run the other aspects of Sleuth.")
        full_model <- NULL
    }
    PPP_kal_dirs <<- kal_dirs <- sample_to_covariates$path
    sample_to_covariates$path <- NULL
    sleuth:::msg("dropping unused factor levels")
    sample_to_covariates <- droplevels(sample_to_covariates)
    nsamp <- 0
    PPP_kal_read_kallisto_raw <<- list()
    PPP_kal_abundance <<- list()
    PPP_kal_list <<- kal_list <- lapply(seq_along(kal_dirs), function(i) {
        nsamp <- sleuth:::dot(nsamp)
        path <- kal_dirs[i]
        suppressMessages({
            PPP_kal_read_kallisto_raw[[i]] <<- kal <- sleuth:::read_kallisto(path, read_bootstrap = FALSE, 
                max_bootstrap = max_bootstrap)
        })
        PPP_kal_abundance[[i]] <<- kal$abundance <- dplyr::mutate(kal$abundance, sample = sample_to_covariates$sample[i])
        kal
    })
    sleuth:::msg("")
    check_result <- sleuth:::check_kal_pack(kal_list)
    kal_versions <- check_result$versions
    PPP_obs_raw1 <<- obs_raw <- dplyr::bind_rows(lapply(kal_list, function(k) k$abundance))
    design_matrix <- NULL
    if (is(full_model, "formula")) {
        design_matrix <- model.matrix(full_model, sample_to_covariates)
    }
    else if (is(full_model, "matrix")) {
        if (is.null(colnames(full_model))) {
            stop("If matrix is supplied, column names must also be supplied.")
        }
        design_matrix <- full_model
    }
    if (!is.null(full_model)) {
        rownames(design_matrix) <- sample_to_covariates$sample
    }
    PPP_obs_raw2 <<- obs_raw <- dplyr::arrange(obs_raw, target_id, sample)
    if (!is.null(target_mapping)) {
        tmp_names <- data.frame(target_id = kal_list[[1]]$abundance$target_id, 
            stringsAsFactors = FALSE)
        PPP_target_mapping <<- target_mapping <- sleuth:::check_target_mapping(tmp_names, target_mapping)
        rm(tmp_names)
    }
    ret <- list(kal = kal_list, kal_versions = kal_versions, 
        obs_raw = obs_raw, sample_to_covariates = sample_to_covariates, 
        bootstrap_summary = NA, full_formula = full_model, design_matrix = design_matrix, 
        target_mapping = target_mapping, gene_mode = !is.null(aggregation_column), 
        gene_column = aggregation_column, transform_fun = transformation_function)
    normalize <- TRUE
    if (normalize) {
        sleuth:::msg("normalizing est_counts")
        PPP_est_counts_spread <<- est_counts_spread <- sleuth:::spread_abundance_by(obs_raw, "est_counts", 
            sample_to_covariates$sample)
        filter_bool <- apply(est_counts_spread, 1, filter_fun, 
            ...)
        PPP_filter_true <<- filter_true <- filter_bool[filter_bool]
        sleuth:::msg(paste0(sum(filter_bool), " targets passed the filter"))
        PPP_est_counts_sf <<- est_counts_sf <- norm_fun_counts(est_counts_spread[filter_bool, 
            , drop = FALSE])
        PPP_filter_df <<- filter_df <- sleuth:::adf(target_id = names(filter_true))
        PPP_est_counts_norm1 <<- est_counts_norm <- sleuth:::as_df(t(t(est_counts_spread)/est_counts_sf))
        est_counts_norm$target_id <- rownames(est_counts_norm)
        PPP_est_counts_norm2 <<- est_counts_norm <- tidyr::gather(est_counts_norm, sample, 
            est_counts, -target_id)
        obs_norm <- est_counts_norm
        obs_norm$target_id <- as.character(obs_norm$target_id)
        obs_norm$sample <- as.character(obs_norm$sample)
        rm(est_counts_norm)
        
        
        sleuth:::msg("normalizing tpm")
        PPP_tpm_spread <<- tpm_spread <- sleuth:::spread_abundance_by(obs_raw, "tpm", sample_to_covariates$sample)
        PPP_tpm_sf <<- tpm_sf <- norm_fun_tpm(tpm_spread[filter_bool, , drop = FALSE])
        PPP_tpm_norm1 <<- tpm_norm <- sleuth:::as_df(t(t(tpm_spread)/tpm_sf))
        tpm_norm$target_id <- rownames(tpm_norm)
        PPP_tpm_norm2 <<- tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)
        tpm_norm$sample <- as.character(tpm_norm$sample)

        sleuth:::msg("merging in metadata")
        PPP_merge_obs_norm <<- obs_norm <- dplyr::arrange(obs_norm, target_id, sample)
        PPP_merge_tpm_norm <<- tpm_norm <- dplyr::arrange(tpm_norm, target_id, sample)
        stopifnot(all.equal(obs_raw$target_id, obs_norm$target_id) && 
            all.equal(obs_raw$sample, obs_norm$sample))
        suppressWarnings({
            if (!all.equal(dplyr::select(obs_norm, target_id, 
                sample), dplyr::select(tpm_norm, target_id, sample))) {
                stop("Invalid column rows. In principle, can simply join. Please report error.")
            }
            PPP_bindcol_obs_norm1 <<- obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(tpm_norm, 
                tpm))
        })
        PPP_bindcol_obs_norm2 <<- obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(obs_raw, 
            eff_len, len))
        PPP_bindcol_obs_norm3 <<- obs_norm <- sleuth:::as_df(obs_norm)
        ret$obs_norm <- obs_norm
        ret$est_counts_sf <- est_counts_sf
        ret$filter_bool <- filter_bool
        ret$filter_df <- filter_df
        ret$obs_norm_filt <- dplyr::semi_join(obs_norm, filter_df, 
            by = "target_id")
        ret$tpm_sf <- tpm_sf
        path <- kal_dirs[1]
        kal_path <- sleuth:::get_kallisto_path(path)
        PPP_h5 <<- target_id <- as.character(rhdf5::h5read(kal_path$path, 
            "aux/ids"))
        num_transcripts <- length(target_id)
        ret$bs_quants <- list()
        which_target_id <- ret$filter_df$target_id
        if (ret$gene_mode) { 
            sleuth:::msg(paste0("aggregating by column: ", aggregation_column))
            PPP_target_mapping <<- target_mapping <- data.table::data.table(target_mapping)
            target_mapping[target_mapping[[aggregation_column]] == 
                "", aggregation_column] <- NA
            agg_id <- unique(target_mapping[, aggregation_column, 
                with = FALSE])
            agg_id <- agg_id[[1]]
            PPP_agg_id <<- agg_id <- agg_id[!is.na(agg_id)]
            mappings <- dplyr::select_(target_mapping, "target_id", 
                aggregation_column)
            PPP_mappings <<- mappings <- data.table::as.data.table(mappings)
            PPP_which_tms <<- which_tms <- which(mappings$target_id %in% which_target_id)
            which_agg_id <- unique(mappings[which_tms, aggregation_column, 
                with = FALSE])
            which_agg_id <- which_agg_id[[1]]
            PPP_which_agg_id <<- which_agg_id <- which_agg_id[!is.na(which_agg_id)]
            PPP_filter_df <<- filter_df <- sleuth:::adf(target_id = which_agg_id)
            PPP_filter_bool <<- filter_bool <- agg_id %in% which_agg_id
            sleuth:::msg(paste0(length(which_agg_id), " genes passed the filter"))
            norm_by_length <- TRUE
            tmp <- data.table::as.data.table(ret$obs_raw)
            PPP_tmp1 <<- tmp <- merge(tmp, mappings, by = "target_id", all.x = TRUE)
            PPP_scale_factor <<- scale_factor <- tmp[, `:=`(scale_factor, median(eff_len)), 
                by = list(sample, eval(parse(text = aggregation_column)))]
            PPP_obs_norm_gene <<- obs_norm_gene <- sleuth:::reads_per_base_transform(ret$obs_norm, 
                scale_factor, aggregation_column, mappings, norm_by_length)
            tmp <- data.table::as.data.table(tpm_norm)
            PPP_tmp2 <<- tmp <- merge(tmp, mappings, by = "target_id", all.x = T)
            if (any(is.na(tmp[[aggregation_column]]))) {
                rows_to_remove <- is.na(tmp[[aggregation_column]])
                num_missing <- length(unique(tmp[rows_to_remove, 
                  target_id]))
                warning(num_missing, " target_ids are missing annotations for the aggregation_column: ", 
                  aggregation_column, ".\nThese target_ids will be dropped from the gene-level analysis.", 
                  "\nIf you did not expect this, check your 'target_mapping' table for missing values.")
                PPP_tmp3 <<- tmp <- tmp[!rows_to_remove]
            }
            tpm_norm_gene <- tmp[, j = list(tpm = sum(tpm)), 
                by = list(sample, eval(parse(text = aggregation_column)))]
            data.table::setnames(tpm_norm_gene, "parse", "target_id")
            tpm_norm_gene <- as_df(tpm_norm_gene)
            obs_norm_gene <- dplyr::arrange(obs_norm_gene, target_id, 
                sample)
            PPP_tpm_norm_gene <<- tpm_norm_gene <- dplyr::arrange(tpm_norm_gene, target_id, 
                sample)
            stopifnot(all.equal(dplyr::select(obs_norm_gene, 
                target_id, sample), dplyr::select(tpm_norm_gene, 
                target_id, sample)))
            suppressWarnings({
                if (!all.equal(dplyr::select(obs_norm_gene, target_id, 
                  sample), dplyr::select(tpm_norm_gene, target_id, 
                  sample), check.attributes = FALSE)) {
                  stop("Invalid column rows. In principle, can simply join. Please report error.")
                }
                PPP_obs_norm_gene <<- obs_norm_gene <- dplyr::bind_cols(obs_norm_gene, 
                  dplyr::select(tpm_norm_gene, tpm))
            })
            ret$filter_df <- sleuth:::adf(target_id = which_agg_id)
            ret$filter_bool <- agg_id %in% which_agg_id
            ret$obs_norm <- obs_norm_gene
            ret$obs_norm_filt <- dplyr::semi_join(obs_norm_gene, 
                filter_df, by = "target_id")
            rm(obs_norm, tpm_norm, obs_norm_gene, tpm_norm_gene)
            PPP_all_sample_bootstrap <<- all_sample_bootstrap <- matrix(NA_real_, nrow = length(which_agg_id), 
                ncol = length(ret$kal))
            which_ids <- which_agg_id
        }
        else {
            all_sample_bootstrap <- matrix(NA_real_, nrow = length(which_target_id), 
                ncol = length(ret$kal))
            which_ids <- which_target_id
        }
        sleuth:::msg("summarizing bootstraps")
        apply_function <- if (num_cores == 1) {
            lapply
        }
        else {
            function(x, y) parallel::mclapply(x, y, mc.cores = num_cores)
        }
        PPP_bs_results <<- bs_results <- apply_function(seq_along(kal_dirs), function(i) {
            samp_name <- sample_to_covariates$sample[i]
            kal_path <- sleuth:::get_kallisto_path(kal_dirs[i])
            sleuth:::process_bootstrap(i, samp_name, kal_path, num_transcripts, 
                est_counts_sf[[i]], read_bootstrap_tpm, ret$gene_mode, 
                extra_bootstrap_summary, target_id, mappings, 
                which_ids, ret$gene_column, ret$transform_fun)
        })
        error_status <- sapply(bs_results, function(x) is(x, 
            "try-error"))
        if (any(error_status)) {
            print(attributes(bs_results[error_status])$condition)
            stop("At least one core from mclapply had an error. See the above error message(s) for more details.")
        }
        PPP_indices <<- indices <- sapply(bs_results, function(result) result$index)
        stopifnot(identical(indices, order(indices)))
        if (read_bootstrap_tpm | extra_bootstrap_summary) {
            ret$bs_quants <- lapply(bs_results, function(result) result$bs_quants)
            names(ret$bs_quants) <- sample_to_covariates$sample
        }
        PPP_all_sample_bootstrap <<- all_sample_bootstrap <- sapply(bs_results, function(result) result$bootstrap_result)
        rownames(all_sample_bootstrap) <- which_ids
        sleuth:::msg("")
        PPP_sigma_q_sq <<- sigma_q_sq <- rowMeans(all_sample_bootstrap)
        if (ret$gene_mode) {
            names(sigma_q_sq) <- which_agg_id
            obs_counts <- sleuth:::obs_to_matrix(ret, "scaled_reads_per_base")[which_agg_id, 
                , drop = FALSE]
        }
        else {
            names(sigma_q_sq) <- which_target_id
            obs_counts <- sleuth:::obs_to_matrix(ret, "est_counts")[which_target_id, 
                , drop = FALSE]
        }
        PPP_final_sigma_q_sq <<- sigma_q_sq <- sigma_q_sq[order(names(sigma_q_sq))]
        obs_counts <- ret$transform_fun(obs_counts)
        PPP_funal_obs_counts <<- obs_counts <- obs_counts[order(rownames(obs_counts)), 
            ]
        ret$bs_summary <- list(obs_counts = obs_counts, sigma_q_sq = sigma_q_sq)
    }
    class(ret) <- "sleuth"
    ret
}



####################variation of the code############
FFF_sleuth_fit <- function (obj, formula = NULL, fit_name = NULL, ...) 
{
    stopifnot(is(obj, "sleuth"))
    if (is.null(formula)) {
        FFF_formula <<- formula <- obj$full_formula
        if (is.null(formula)) {
            stop("'formula' was not specified and the 'full' model was not specified in `sleuth_prep`.", 
                " Please specify a formula and a label.")
        }
    }
    else if (!is(formula, "formula") && !is(formula, "matrix")) {
        stop("'", substitute(formula), "' is not a valid 'formula' or 'matrix'")
    }
    if (is.null(fit_name)) {
        fit_name <- "full"
    }
    else if (!is(fit_name, "character")) {
        stop("'", substitute(fit_name), "' is not a valid 'character'")
    }
    if (length(fit_name) > 1) {
        stop("'", substitute(fit_name), "' is of length greater than one.", 
            " Please only supply one string.")
    }
    X <- NULL
    if (is(formula, "formula")) {
        X <- model.matrix(formula, obj$sample_to_covariates)
    }
    else {
        if (is.null(colnames(formula))) {
            stop("If matrix is supplied, column names must also be supplied.")
        }
        FFF_X <<- X <- formula
    }
    FFF_rnamesX <<- rownames(X) <- obj$sample_to_covariates$sample
    FFF_A <<- A <- solve(t(X) %*% X)
    sleuth:::msg("fitting measurement error models")
    FFF_mes <<- mes <- sleuth:::me_model_by_row(obj, X, obj$bs_summary)
    tid <- names(mes)
    mes_df <- dplyr::bind_rows(lapply(mes, function(x) {
        data.frame(rss = x$rss, sigma_sq = x$sigma_sq, sigma_q_sq = x$sigma_q_sq, 
            mean_obs = x$mean_obs, var_obs = x$var_obs)
    }))
    mes_df$target_id <- tid
    FFF_tid <<- mes_df
    rm(tid)
    FFF_mes_df <<- mes_df <- dplyr::mutate(mes_df, sigma_sq_pmax = pmax(sigma_sq, 
        0))
    sleuth:::msg("shrinkage estimation")
    FFF_swg <<- swg <- sleuth:::sliding_window_grouping(mes_df, "mean_obs", "sigma_sq_pmax", 
        ignore_zeroes = TRUE, ...)
    FFF_shrink1 <<- l_smooth <- sleuth:::shrink_df(swg, sqrt(sqrt(sigma_sq_pmax)) ~ mean_obs, 
        "iqr", ...)
    FFF_shrink_select <<- l_smooth <- dplyr::select(dplyr::mutate(l_smooth, smooth_sigma_sq = shrink^4), 
        -shrink)
    FFF_shrink_mutate <<- l_smooth <- dplyr::mutate(l_smooth, smooth_sigma_sq_pmax = pmax(smooth_sigma_sq, 
        sigma_sq))
    sleuth:::msg("computing variance of betas")
    beta_covars <- lapply(1:nrow(l_smooth), function(i) {
        row <- l_smooth[i, ]
        with(row, sleuth:::covar_beta(smooth_sigma_sq_pmax + sigma_q_sq, 
            X, A))
    })
    names(beta_covars) <- l_smooth$target_id
    FFF_beta_covars <<- beta_covars
    if (is.null(obj$fits)) {
        obj$fits <- list()
    }
    obj$fits[[fit_name]] <- list(models = mes, summary = l_smooth, 
        beta_covars = beta_covars, formula = formula, design_matrix = X, 
        transform_synced = TRUE)
    class(obj$fits[[fit_name]]) <- "sleuth_model"
    obj
}



RRR_sleuth_fit <- function (obj, formula = NULL, fit_name = NULL, ...) 
{
    stopifnot(is(obj, "sleuth"))
    if (is.null(formula)) {
        RRR_formula <<- formula <- obj$full_formula
        if (is.null(formula)) {
            stop("'formula' was not specified and the 'full' model was not specified in `sleuth_prep`.", 
                " Please specify a formula and a label.")
        }
    }
    else if (!is(formula, "formula") && !is(formula, "matrix")) {
        stop("'", substitute(formula), "' is not a valid 'formula' or 'matrix'")
    }
    if (is.null(fit_name)) {
        fit_name <- "full"
    }
    else if (!is(fit_name, "character")) {
        stop("'", substitute(fit_name), "' is not a valid 'character'")
    }
    if (length(fit_name) > 1) {
        stop("'", substitute(fit_name), "' is of length greater than one.", 
            " Please only supply one string.")
    }
    X <- NULL
    if (is(formula, "formula")) {
        X <- model.matrix(formula, obj$sample_to_covariates)
    }
    else {
        if (is.null(colnames(formula))) {
            stop("If matrix is supplied, column names must also be supplied.")
        }
        RRR_X <<- X <- formula
    }
    RRR_rnamesX <<- rownames(X) <- obj$sample_to_covariates$sample
    RRR_A <<- A <- solve(t(X) %*% X)
    sleuth:::msg("fitting measurement error models")
    RRR_mes <<- mes <- sleuth:::me_model_by_row(obj, X, obj$bs_summary)
    tid <- names(mes)
    mes_df <- dplyr::bind_rows(lapply(mes, function(x) {
        data.frame(rss = x$rss, sigma_sq = x$sigma_sq, sigma_q_sq = x$sigma_q_sq, 
            mean_obs = x$mean_obs, var_obs = x$var_obs)
    }))
    mes_df$target_id <- tid
    RRR_tid <<- mes_df
    rm(tid)
    RRR_mes_df <<- mes_df <- dplyr::mutate(mes_df, sigma_sq_pmax = pmax(sigma_sq, 
        0))
    sleuth:::msg("shrinkage estimation")
    RRR_swg <<- swg <- sleuth:::sliding_window_grouping(mes_df, "mean_obs", "sigma_sq_pmax", 
        ignore_zeroes = TRUE, ...)
    RRR_shrink1 <<- l_smooth <- sleuth:::shrink_df(swg, sqrt(sqrt(sigma_sq_pmax)) ~ mean_obs, 
        "iqr", ...)
    RRR_shrink_select <<- l_smooth <- dplyr::select(dplyr::mutate(l_smooth, smooth_sigma_sq = shrink^4), 
        -shrink)
    RRR_shrink_mutate <<- l_smooth <- dplyr::mutate(l_smooth, smooth_sigma_sq_pmax = pmax(smooth_sigma_sq, 
        sigma_sq))
    sleuth:::msg("computing variance of betas")
    beta_covars <- lapply(1:nrow(l_smooth), function(i) {
        row <- l_smooth[i, ]
        with(row, sleuth:::covar_beta(smooth_sigma_sq_pmax + sigma_q_sq, 
            X, A))
    })
    names(beta_covars) <- l_smooth$target_id
    RRR_beta_covars <<- beta_covars
    if (is.null(obj$fits)) {
        obj$fits <- list()
    }
    obj$fits[[fit_name]] <- list(models = mes, summary = l_smooth, 
        beta_covars = beta_covars, formula = formula, design_matrix = X, 
        transform_synced = TRUE)
    class(obj$fits[[fit_name]]) <- "sleuth_model"
    obj
}






LLL_sleuth_lrt <- function (obj, null_model, alt_model) 
{
    stopifnot(is(obj, "sleuth"))
    sleuth:::model_exists(obj, null_model)
    sleuth:::model_exists(obj, alt_model)
    if (!obj$fits[[alt_model]]$transform_synced) {
        stop("Model '", alt_model, "' was not computed using the sleuth object's", 
             " current transform function. Please rerun sleuth_fit for this model.")
    }
    if (!obj$fits[[null_model]]$transform_synced) {
        stop("Model '", null_model, "' was not computed using the sleuth object's", 
             " current transform function. Please rerun sleuth_fit for this model.")
    }
    if (!sleuth:::likelihood_exists(obj, null_model)) {
        LLL_null_obj <<- obj <- sleuth:::compute_likelihood(obj, null_model)
    }
    if (!sleuth:::likelihood_exists(obj, alt_model)) {
        LLL_alt_obj <<- obj <- sleuth:::compute_likelihood(obj, alt_model)
    }
    LLL_n_ll <<- n_ll <- sleuth:::get_likelihood(obj, null_model)
    LLL_a_ll <<- a_ll <- sleuth:::get_likelihood(obj, alt_model)
    obj$fits[[null_model]]$likelihood <- n_ll
    obj$fits[[alt_model]]$likelihood <- a_ll
    LLL_test_statistic <<- test_statistic <- 2 * (a_ll - n_ll)
    LLL_degrees_free <<- degrees_free <- obj$fits[[null_model]]$models[[1]]$ols_fit$df.residual - 
        obj$fits[[alt_model]]$models[[1]]$ols_fit$df.residual
    LLL_p_value <<- p_value <- pchisq(test_statistic, degrees_free, lower.tail = FALSE)
    LLL_result_adf <<-result <- sleuth:::adf(target_id = names(obj$fits[[alt_model]]$likelihood), 
                  test_stat = test_statistic, pval = p_value)
    LLL_result_adf_post <<- result <- dplyr::mutate(result, qval = p.adjust(pval, method = "BH"))
    LLL_model_info1 <<- model_info <- data.table::data.table(obj$fits[[null_model]]$summary)
    LLL_model_info2 <<- model_info <- dplyr::select(model_info, -c(iqr))
    LLL_result_join <<- result <- dplyr::left_join(data.table::data.table(result), 
                               model_info, by = "target_id")
    LLL_result_mutate <<- result <- dplyr::mutate(result, degrees_free = degrees_free)
    test_name <- paste0(null_model, ":", alt_model)
    obj <- sleuth:::add_test(obj, result, test_name, "lrt")
    obj
}
