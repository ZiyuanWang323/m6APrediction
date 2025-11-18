# m6APrediction

### Overview  

m6APrediction is an open-source R package that delivers **fast, accurate, and ready-to-use** prediction of *N6-methyladenosine* (m6A) modification sites in RNA sequences. Built on a rigorously trained **random-forest** model, the package provides both **batch** and **single-sequence** functions, enabling researchers to screen thousands of candidate sites or analyze individual transcripts in seconds. All core assets—pre-trained model, example data, and documentation—are bundled inside, so no external downloads or bioinformatics pipelines are required.

### Purpose  
The rapid accumulation of high-throughput sequencing data has made m6A the most intensively studied RNA modification, yet experimental mapping (miCLIP, MeRIP-seq) remains costly and time-consuming. m6APrediction fills this gap by offering **reliable in-silico** predictions using only **basic sequence and genomic features** (5-mer frequency, GC content, RNA region, evolutionary conservation, etc.). By converting complex machine-learning steps into two simple R commands, the package helps biologists:

* prioritize candidate sites for downstream validation  
* explore m6A distribution in novel transcriptomes or species  
* integrate methylation potential into broader gene-regulation analyses  

without requiring prior expertise in bioinformatics or model building.
## Installation

```r
# Option 1: with devtools
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yourusername/m6APrediction")

# Option 2: with remotes (lighter dependency)
if (!require("remotes")) install.packages("remotes")
remotes::install_github("yourusername/m6APrediction")

Quickly Start
library(m6APrediction)

# 1. Pull in the bundled model and demo file
model   <- readRDS(system.file("extdata", "rf_fit.rds",        package = "m6APrediction"))
demo_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))

# 2. Score the entire demo set at once
batch_out <- prediction_multiple(
  model             = model,
  feature_table     = demo_df,
  cutoff            = 0.5               # decision threshold
)
print(head(batch_out[, c("DNA_5mer", "m6A_prob", "m6A_call")]))

# 3. Rapid check on a single transcript
single_out <- prediction_single(
  model           = model,
  gc_content      = 0.51,
  RNA_type        = "lncRNA",
  RNA_region      = "3'UTR",
  exon_length     = 320,
  junction_dist   = 75,
  conserv_score   = 0.82,
  fivemer         = "GCTAG",
  cutoff          = 0.5
)
str(single_out)   # compact view of the returned list


