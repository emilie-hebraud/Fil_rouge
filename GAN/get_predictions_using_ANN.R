get_predictions_using_ANN_batched <- function(raw_profile, profile_reading_ANN, batch_size = 32){
  
  source('getmode.R')
  
  # ==================== ADAPTATION POUR MODÈLE 6 DYES ====================
  # Le modèle attend 6 dyes, mais nous avons 4 channels
  # Solution : Ajouter 2 channels vides (remplis de zéros)
  
  number_dyes_in_data <- 4      # Vos vraies données
  number_dyes_in_model <- 6     # Ce que le modèle attend
  expected_input_scans <- 5200  # Taille attendue par le modèle
  context_window <- 100
  saturation <- 33000
  
  # =======================================================================
  
  cat("Configuration détectée:\n")
  cat("  - Dyes dans les données:", number_dyes_in_data, "\n")
  cat("  - Dyes attendus par le modèle:", number_dyes_in_model, "\n")
  cat("  - Scans attendus:", expected_input_scans, "\n")
  cat("  - Context window:", context_window, "\n\n")
  
  # Calculer les scans analysables
  startScan <- 1 + context_window
  endScan <- nrow(raw_profile) - context_window
  total_scans <- endScan - startScan + 1
  
  if (total_scans <= 0) {
    stop(sprintf("ERREUR: Profil trop court (%d scans). Minimum requis: %d", 
                 nrow(raw_profile), 2 * context_window + 1))
  }
  
  # Initialiser les prédictions (seulement pour les 4 vrais channels)
  predictions <- array("X", dim = c(nrow(raw_profile), number_dyes_in_data))
  
  # Catégories de sortie par dye (du modèle original à 6 dyes)
  outputs_info_per_dye <-  list(
    list(categories = c("B", "P", "S")),
    list(categories = c("B", "P", "S")),
    list(categories = c("B", "P", "S")),
    list(categories = c("B", "P", "S")))
    
    
  
  # Pré-calculer les modes pour les 4 vrais channels
  dye_modes <- sapply(1:number_dyes_in_data, function(idye) {
    getmode(raw_profile[, idye + 1])
  })
  
  # Traiter par batch
  scan_positions <- startScan:endScan
  n_batches <- ceiling(total_scans / batch_size)
  
  cat(sprintf("Traitement de %d scans en %d batches de %d...\n\n", 
              total_scans, n_batches, batch_size))
  
  for (iBatch in 1:n_batches){
    
    # Indices du batch
    batch_start <- (iBatch - 1) * batch_size + 1
    batch_end <- min(iBatch * batch_size, total_scans)
    batch_scans <- scan_positions[batch_start:batch_end]
    current_batch_size <- length(batch_scans)
    
    if (iBatch %% 10 == 0 || iBatch == 1 || iBatch == n_batches) {
      cat(sprintf("  Batch %d/%d (scans %d-%d)\r", 
                  iBatch, n_batches, batch_scans[1], batch_scans[current_batch_size]))
    }
    
    # ============ CONSTRUIRE L'INPUT POUR LE MODÈLE ============
    # Forme attendue: (batch_size, 6, 5200, 1)
    
    batch_data <- array(0, dim = c(current_batch_size, 
                                   number_dyes_in_model, 
                                   expected_input_scans, 
                                   1))
    
    for (i in 1:current_batch_size){
      iScan <- batch_scans[i]
      window_start <- iScan - context_window
      window_end <- iScan + context_window
      window_size <- window_end - window_start + 1  # Devrait être 201
      
      # Remplir les 4 premiers channels avec vos vraies données
      for (idye in 1:number_dyes_in_data){
        signal_window <- raw_profile[window_start:window_end, idye + 1]
        signal_normalized <- (signal_window - dye_modes[idye]) / saturation
        
        # Centrer la fenêtre dans l'input de taille 5200
        # Position centrale
        center_pos <- floor(expected_input_scans / 2)
        start_pos <- center_pos - floor(window_size / 2) + 1
        end_pos <- start_pos + window_size - 1
        
        batch_data[i, idye, start_pos:end_pos, 1] <- signal_normalized
      }
      
      # Les channels 5 et 6 restent à zéro (pas de données)
    }
    
    # ============ PRÉDICTION ============
    # Predict entire batch at once
    batch_data_tensor <- tensorflow::tf$constant(batch_data, dtype = tensorflow::tf$float32)
    batch_predictions <- profile_reading_ANN$predict(batch_data_tensor, verbose = 0L)
    
    # ============ DÉCODER LES PRÉDICTIONS ============
    # On garde seulement les prédictions des 4 premiers channels
    for (i in 1:current_batch_size){
      iScan <- batch_scans[i]
      
      for (idye in 1:number_dyes_in_data){
        # Extraire les probabilités pour ce dye
        # batch_predictions est une liste de sorties, une par dye
        probs <- batch_predictions[[idye]][i, 1, 1, ]
        
        # Trouver la catégorie avec la probabilité maximale
        best_index <- which.max(probs)
        predictions[iScan, idye] <- outputs_info_per_dye[[idye]]$categories[best_index]
      }
    }
    
  } # end batch loop
  
  cat("\n✓ Traitement terminé\n")
  return(predictions)
}