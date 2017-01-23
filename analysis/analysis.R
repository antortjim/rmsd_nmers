# Analysis of fragments results from protein.py
setwd("/home/antortjim/MEGA/Master/SB/exam")
library("ggplot2")
library("dplyr")
theme_set(theme_bw(base_size = 20))
flen <- 5
out_folder = "shiny/sb_exam/out/"
plot_folder = "plots/"

signal_fn <- paste(out_folder, flen, "-mers_rmsd.txt", sep = "")
random_fn <- paste(out_folder, flen, "-mers_random.txt", sep = "")
fragments_fn <- paste(out_folder, flen, "-mer_fragments.csv", sep = "")

print("Reading data")
signal <- read.csv(file = signal_fn, header = T)
signal$pair <- "signal"
random <- read.csv(file = random_fn, header = T)
random$pair <- "random"
fragments <- read.csv(file = fragments_fn, header = T)



check_intraprotein <- function(my_row) {
  if (my_row[2] == my_row[6]) {
    if(my_row[3] == my_row[7] & my_row[4] == my_row[8]) {
      return(NA)
    } else {
      return(TRUE)
    }
  } else {
    return(FALSE) 
  }
}

signal$intraprotein <- apply(signal, 1, check_intraprotein)
#Working with interprotein only
signal <- filter(signal, intraprotein == FALSE) # keep only interprotein

# Should duplicates be excluded?


idx <- sample(x = 1:nrow(random), size = nrow(signal))
random <- random[idx, ]

my_data <- full_join(signal, random)


print("Plotting histograms")
ggplot(data = my_data, aes(rmsd, fill = pair)) +
  geom_histogram(bins = 60) +
  facet_wrap( ~ pair)
fn <- paste(plot_folder, "duplicates_histogram.png", sep = "")
ggsave(fn)
# 
# 
# 
# ggplot(data = my_data, aes(rmsd, col = pair)) +
#   geom_density(stat = "density", lwd = 1.2)
# 
# fn <- paste(plot_folder, "duplicates_density.png", sep = "")
# ggsave(fn)

# Select residues in pdb4ptp and pdb1mct
pdb1 <- "pdb4ptp.ent"
pdb2 <- "pdb1mct.ent"


(my_data %>% filter(f1.seq_id == pdb1 &
                     f2.seq_id == pdb2) %>%
  arrange(f1.start))$f1.start

(my_data %>% filter(f1.seq_id == pdb1 &
                      f2.seq_id == pdb2) %>%
  arrange(f2.start))$f2.start


pdb1 <- "pdb1aru.ent"
pdb2 <- "pdb1nif.ent"
####################################################################

# Reload data and repeat analysis without excluding duplicates as above
print("Reading data")
signal <- read.csv(file = signal_fn, header = T)
signal$pair <- "signal"
random <- read.csv(file = random_fn, header = T)
random$pair <- "random"
fragments <- read.csv(file = fragments_fn, header = T)
#my_data <- full_join(signal, random)


signal$intraprotein <- apply(signal, 1, check_intraprotein)
#Working with interprotein only
signal <- filter(signal, intraprotein == FALSE) # keep only interprotein

# Remove duplicates
# Duplicate = two or more entries where f1 and f2.seq_id
# are identical
signal <- signal %>%
  group_by(f1.seq_id, f2.seq_id) %>%
  filter(row_number() == 1)


idx <- sample(x = 1:nrow(random), size = nrow(signal))
random <- random[idx, ]

my_data <- full_join(signal, random)


print("Plotting histograms")
ggplot(data = my_data, aes(rmsd, fill = pair)) +
  geom_histogram(bins = 60) +
  facet_wrap( ~ pair)
fn <- paste(plot_folder, "no_duplicates_histogram.png", sep = "")
ggsave(fn)



ggplot(data = my_data, aes(rmsd, col = pair)) +
  geom_density(stat = "density", lwd = 1.2)

fn <- paste(plot_folder, "no_duplicates_density.png", sep = "")
ggsave(fn)



