# Analysis of fragments results from protein.py
# Load libraries and set variables
setwd("/home/antortjim/MEGA/Master/SB/exam")
library("ggplot2")
library("data.table")
library("dplyr")
library("cowplot")
theme_set(theme_bw(base_size = 20))
flen <- 5
out_folder = "../shiny/out/"
plot_folder = "../..plots/"

# check_intraprotein <- function(my_row) {
#   if (my_row[2] == my_row[6]) {
#     if(my_row[3] == my_row[7] & my_row[4] == my_row[8]) {
#       return(NA)
#     } else {
#       return(TRUE)
#     }
#   } else {
#     return(FALSE) 
#   }
# }


# read data
sequence_fn <- paste(out_folder, flen, "-mers_rmsd.txt", sep = "")
random_fn <- paste(out_folder, flen, "-mers_random.txt", sep = "")
fragments_fn <- paste(out_folder, flen, "-mer_fragments.csv", sep = "")

print("Reading data")
sequence <- read.csv(file = sequence_fn, header = T, stringsAsFactors = F)
sequence$pair <- "sequence"
random <- read.csv(file = random_fn, header = T, stringsAsFactors = F)
random$pair <- "random"
fragments <- read.csv(file = fragments_fn, header = T, stringsAsFactors = F)



# Pymol
# Select residues in pdb4ptp and pdb1mct
# pdb1 <- "pdb4ptp.ent"
# pdb2 <- "pdb1mct.ent"
# 
#

# (my_data %>% filter(f1.seq_id == pdb1 &
#                      f2.seq_id == pdb2) %>%
#   arrange(f1.start))$f1.start
# 
# (my_data %>% filter(f1.seq_id == pdb1 &
#                       f2.seq_id == pdb2) %>%
#   arrange(f2.start))$f2.start
# 
# 
#pdb1 <- "pdb1cpc.ent"
#pdb2 <- "pdb2rhe.ent"
####################################################################
# Mark duplicates
# Duplicate = two or more entries where f1 and f2.seq_id
# are identical or
# entries where f1 and f2.seq_id are the same
sequence$id <- 1:nrow(sequence)
sequence$duplicated <- FALSE
sequence[sequence$f1.seq_id == sequence$f2.seq_id, 12] <- TRUE
DT <- as.data.table(sequence)
DT[DT[,.SD[2,.N], by=.(f1.seq_id, f2.seq_id)]$id, 12] <- TRUE
sequence <- as.data.frame(DT)
# mark the same number of random entries
random$duplicated <- c(rep(T, sum(sequence$duplicated)),
                       rep(F, nrow(random) - sum(sequence$duplicated)))




# Merge sequence and random data
my_data <- full_join(sequence, as.data.table(random))

# Use formatted labels in facetting
pair_names <- c(
  'random' = "Random pairs",
  'sequence' = "Sequence pairs")

print("Plotting histograms")
A <- ggplot(data = my_data, aes(rmsd, fill = pair)) +
  geom_histogram(bins = 60) +
  facet_wrap( ~ pair, labeller = as_labeller(pair_names)) +
  guides(fill = FALSE) +
  labs(x = "RMSD", y = "Count")
# Histogram including duplicated entries. Crazy amount of sequence
# pairs with low RMSD
A

# Only non duplicated
ndDT <- my_data %>% filter(duplicated == F)


B <- ggplot(data = ndDT, aes(rmsd, fill = pair)) +
  geom_histogram(bins = 60) +
  facet_wrap( ~ pair, labeller = as_labeller(pair_names)) +
  guides(fill = FALSE) +
  labs(x = "RMSD", y = "Count")

B

# check that there are as many random pairs as sequence pairs
(my_data %>% filter(duplicated == F))$pair %>% table
# TRUE :)

histogram <- plot_grid(A, B, labels = c("A", "B"), nrow = 2, align = "v")


fn <- paste(plot_folder, "histogram.png", sep = "")
# Save exam figure
save_plot(fn, histogram,
          ncol = 1, # we're saving a grid plot of 1 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 0.5
          base_aspect_ratio = 3
)


# Obtain statistics
ndDT %>% group_by(pair) %>% summarise(x = sum(rmsd < 2))
ndDT %>% group_by(pair) %>% summarise(x = sum(rmsd > 2))
ndDT %>% group_by(pair) %>% summarise(x = median(rmsd))
ndDT %>% group_by(pair) %>% summarise(x = mean(rmsd))
ndDT %>% group_by(pair) %>% summarise(x = table(pair))
sum(ndDT$pair == "sequence")
ndDT %>% group_by(pair) %>% summarise(x = mean(rmsd))



C <- ggplot(data = ndDT, aes(rmsd, col = pair, fill = pair)) +
  geom_density(alpha = 0.5) +
#  facet_wrap( ~ pair, labeller = as_labeller(pair_names)) +
  guides(color = FALSE) +
  labs(x = "RMSD", y = "Density") +
  scale_fill_discrete(name = "Pairs",
                      breaks = c("random", "sequence"),
                      labels = c("Random", "Sequence"))


D <- ggplot(data = ndDT, aes(x = pair, y = rmsd, fill = pair)) +
  geom_boxplot() +
  labs(x = "", y = "RMSD") +
  guides(fill = FALSE) +
  scale_x_discrete(labels = pair_names)
D

wilcox.test(ndDT %>% filter(pair == "random") %>% select(rmsd) %>% .$rmsd,
            ndDT %>% filter(pair == "sequence") %>% select(rmsd) %>% .$rmsd)$p.value
