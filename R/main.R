#' Encode DNA sequences into factor columns
#'
#' This function encodes fixed-length DNA sequences ( 5-mers)
#' into a data frame where each nucleotide position becomes a separate column.
#'
#' @param dna_strings A character vector of DNA sequences (all must have the same length).
#'
#' @return A data frame with columns named nt_pos1, nt_pos2, etc.,
#' where each column is a factor with levels A, T, C, and G.
#'
#' @keywords internal
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Multiple Sample Prediction for m6A Sites
#'
#' This function predicts m6A probability and status for multiple input samples
#' stored in a data frame containing the required features.
#'
#' @param ml_fit A trained machine learning model (random forest).
#' @param feature_df A data frame containing the input features:
#' gc_content, RNA_type, RNA_region, exon_length, distance_to_junction,
#' evolutionary_conservation, and DNA_5mer.
#' @param positive_threshold A numeric threshold (0–1) for classifying as Positive.
#'
#' @return The same data frame with two additional columns:
#' predicted_m6A_prob and predicted_m6A_status.
#'
#' @import randomForest
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length",
                  "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
                %in% colnames(feature_df)))

  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  dna_encoded <- dna_encoding(feature_df$DNA_5mer)

  predict_data <- data.frame(
    gc_content = feature_df$gc_content,
    RNA_type = feature_df$RNA_type,
    RNA_region = feature_df$RNA_region,
    exon_length = feature_df$exon_length,
    distance_to_junction = feature_df$distance_to_junction,
    evolutionary_conservation = feature_df$evolutionary_conservation
  )

  predict_data <- cbind(predict_data, dna_encoded)

  predicted_probs <- predict(ml_fit, newdata = predict_data, type = "prob")
  predicted_prob_positive <- predicted_probs[, "Positive"]
  predicted_status <- ifelse(predicted_prob_positive >= positive_threshold,
                             "Positive", "Negative")

  feature_df$predicted_m6A_prob <- predicted_prob_positive
  feature_df$predicted_m6A_status <- predicted_status
  return(feature_df)
}


#' Single Sample Prediction for m6A Sites
#'
#' This function takes individual feature values and uses a trained model
#' to predict the probability and classification status for a single m6A site.
#'
#' @param ml_fit A trained machine learning model (random forest).
#' @param gc_content A numeric value for the GC content.
#' @param RNA_type A character string specifying the RNA type (mRNA, lincRNA, lncRNA, pseudogene).
#' @param RNA_region A character string specifying the RNA region (CDS, intron, 3'UTR, 5'UTR).
#' @param exon_length A numeric value for the exon length.
#' @param distance_to_junction A numeric value for the distance to the nearest junction.
#' @param evolutionary_conservation A numeric value for the evolutionary conservation score.
#' @param DNA_5mer A character string representing the 5-mer DNA sequence (e.g. "GGACA").
#' @param positive_threshold A numeric threshold (0–1) for classifying as Positive.
#'
#' @return A named vector with two elements: predicted_m6A_prob and predicted_m6A_status.
#'
#' @import randomForest
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' prediction_single(
#'   ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "ATGCAT",
#'   positive_threshold = 0.5
#' )
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

    encoded_dna <- dna_encoding(feature_df$DNA_5mer)


   feature_df_encoded <- cbind(
    feature_df[c("gc_content", "RNA_type", "RNA_region", "exon_length",
                 "distance_to_junction", "evolutionary_conservation")],
    encoded_dna
  )


  feature_df_encoded$RNA_type <- factor(
    feature_df_encoded$RNA_type,
    levels = if (!is.null(ml_fit$forest$xlevels$RNA_type)) ml_fit$forest$xlevels$RNA_type else c("mRNA","lincRNA","lncRNA","pseudogene")
  )
  feature_df_encoded$RNA_region <- factor(
    feature_df_encoded$RNA_region,
    levels = if (!is.null(ml_fit$forest$xlevels$RNA_region)) ml_fit$forest$xlevels$RNA_region else c("CDS","intron","3'UTR","5'UTR")
  )


  prob_mat <- predict(ml_fit, newdata = feature_df_encoded, type = "prob")
  pred_prob <- if ("Positive" %in% colnames(prob_mat)) as.numeric(prob_mat[, "Positive"]) else as.numeric(prob_mat[, 2])
  pred_status <- ifelse(pred_prob >= positive_threshold, "Positive", "Negative")

  returned_vector <- c(
    predicted_m6A_prob = pred_prob,
    predicted_m6A_status = pred_status
  )
  return(returned_vector)
}
