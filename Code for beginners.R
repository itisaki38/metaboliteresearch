# ───────────────────────────────────────────────────────────────────────────────
# Affinity Analysis Script (Commented for Beginners)
# ───────────────────────────────────────────────────────────────────────────────

# 0) — Install & load required packages
# --------------------------------------
# Core dependencies:
packages <- c("tidyverse", "httr", "jsonlite", "knitr")
# For Excel export: writexl (binary on Windows) or alternative openxlsx
excel_pkgs <- c("writexl")
installed <- c(packages, excel_pkgs) %in% rownames(installed.packages())
if (any(!installed)) {
  # Attempt to install core packages first
  install.packages(packages[!installed], repos = "https://cloud.r-project.org")
  # Install writexl; on Windows, force binary to avoid compilation
  install.packages("writexl", repos = "https://cloud.r-project.org", type = "binary")
  # If writexl still fails, you can install openxlsx instead:
  # install.packages("openxlsx")
}
library(tidyverse)  # data wrangling, string ops, ggplot2
library(httr)       # HTTP requests
library(jsonlite)   # parse JSON
library(knitr)      # kable()
# Try loading writexl; if unavailable, uncomment openxlsx below
library(writexl)
# library(openxlsx)  # alternative Excel export, no compilation needed

# 1) — Define your input affinity files
# -------------------------------------
files <- list(
  Glucose   = "GLC.pdb_affinities.txt",
  Metformin = "Conformer3D_COMPOUND_CID_4091.pdb_affinities.txt",
  Berberine = "CID_2353_berberine1.pdb_affinities.txt"
)
missing <- names(files)[!file.exists(unlist(files))]
if (length(missing) > 0) {
  stop("ERROR: Missing files: ", paste(missing, collapse = ", "))
}

# 2) — Function to parse a single affinity file
# ---------------------------------------------
parse_affinities <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- grep("log\\.log:", lines, value = TRUE)
  pdb_ids <- str_extract(lines, "(?<=/)[^/]+(?=/log\\.log)")
  affinities <- as.numeric(str_match(lines,
                     "log\\.log:\\s*\\d+\\s*(-?[0-9]+\\.?[0-9]*)")[,2])
  tibble(PDB_ID = pdb_ids, affinity = affinities)
}

# 3) — Read & rename each ligand’s data
# -------------------------------------
parsed_list <- imap(files, function(path, ligand) {
  parse_affinities(path) %>%
    filter(!is.na(affinity)) %>%
    rename(!!ligand := affinity)
})

# 4) — Find shared PDB IDs and clean IDs
# --------------------------------------
df_shared <- reduce(parsed_list, inner_join, by = "PDB_ID") %>%
  mutate(PDB_ID = str_remove(PDB_ID, "\\.pdb$"))

# 5) — Top N targets by strongest binding
# ---------------------------------------
TOP_N <- 20

df_top <- df_shared %>%
  rowwise() %>%
  mutate(Average = mean(c_across(-PDB_ID))) %>%
  ungroup() %>%
  arrange(Average) %>%
  slice_head(n = TOP_N)

# 6) — Fetch PDB titles
# ---------------------
get_pdb_title <- function(code4) {
  url <- paste0("https://data.rcsb.org/rest/v1/core/entry/", code4)
  res <- GET(url)
  if (status_code(res) != 200) return(NA_character_)
  content(res, "parsed")$struct$title
}

final_table <- df_top %>%
  mutate(EntryCode = str_remove(PDB_ID, "_.*$")) %>%
  rowwise() %>%
  mutate(Identification = get_pdb_title(EntryCode)) %>%
  ungroup() %>%
  select(
    Target         = PDB_ID,
    Identification,
    !!!syms(names(files))
  )

# 7) — Display in console or RMarkdown
# ------------------------------------
kable(final_table, caption = paste0("Top ", TOP_N, 
                                 " shared targets (strongest average binding)"))

# 8) — Plot per-ligand affinities
# -------------------------------
df_long <- df_top %>%
  pivot_longer(cols = all_of(names(files)), names_to = "Ligand", values_to = "Affinity")

ggplot(df_long, aes(x = PDB_ID, y = Affinity, color = Ligand, group = Ligand)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = paste0("Binding Affinities of Top ", TOP_N, " Shared Targets"),
       x = "Target (PDB ID)", y = "Affinity (kcal/mol)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")

# 9) — Save tables to files
# -------------------------
# CSV exports
write_csv(df_shared,   "df_shared.csv")
write_csv(df_top,      "df_top.csv")
write_csv(final_table, "final_table.csv")

# Excel exports
# Option A: writexl
write_xlsx(list(Shared = df_shared, Top20 = df_top, Final = final_table),
           path = "affinity_results.xlsx")
# Option B: openxlsx (if writexl not available)
# write.xlsx(final_table, file = "final_table_openxlsx.xlsx")

# Tips:
# - If writexl fails to install, uncomment openxlsx lines and comment out writexl.
# - Use write.csv() or saveRDS() for alternative formats.


