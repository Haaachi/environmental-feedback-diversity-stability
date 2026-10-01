packages <- c("tidyverse", "readxl", "cowplot", "patchwork", "scales")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
cat("R dependencies are available.\n")
