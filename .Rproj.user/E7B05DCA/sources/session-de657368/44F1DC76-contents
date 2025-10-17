
required_packages <- c(
  # Core
  "dplyr", "tibble", "rioja", "tidypaleo", "randomForest", "gbm",
  # Data analysis
  "stringr", "reshape2",
  # Visualization
  "ggplot2", "ggnewscale", "RColorBrewer", "patchwork", "ggthemes", "ggpmisc",
  # Reproducibility
  "renv"
)

install_and_load <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("Installing missing package: ", pkg))
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# Run it
install_and_load(required_packages)

message("All required packages are installed and loaded.")

