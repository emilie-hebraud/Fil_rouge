library(here)
read_labels_safe <- function(filepath) {
  raw <- readLines(filepath, warn = FALSE)
  
  # Detect separator and clean \r
  raw <- gsub("\r", "", raw)
  
  # Get header columns
  header <- strsplit(raw[[1]], ";")[[1]]
  
  # Check required columns exist
  if (!all(c("plant_id","multiplex","channel","peak_kind","mu_pb") %in% header)) {
    stop(paste("Missing required columns in:", filepath))
  }
  
  # Parse the file
  df <- read.table(
    text = paste(raw, collapse = "\n"),
    sep = ";",
    header = TRUE,
    dec = ".",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    row.names = NULL,
    fill = TRUE
  )
  
  # Keep only plant_id and mu_pb
  df <- df[, c("plant_id","multiplex","channel","peak_kind","mu_pb")]
  return(df)
}

# Read both files
file1 <- read_labels_safe(paste(here(), "synthetic_data/labels_map/peak_positions_detailed.csv", sep="/"))
file2 <- read_labels_safe(paste(here(), "input_data/peak_position_detailed_raw.csv", sep="/"))

# Concatenate
result <- rbind(file1, file2)

# Save
write.table(
  result,
  file =paste(here(), "input_data/labels.csv", sep="/"),
  sep = ";",
  row.names = FALSE,
  quote = FALSE,
  dec = "."
)

cat("Done:", nrow(result), "rows written\n")