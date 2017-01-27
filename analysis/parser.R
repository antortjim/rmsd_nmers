rm(list = ls())
setwd("~/MEGA/Master/SB/exam/rna/")
suppressMessages(library("dplyr"))
library("ggplot2")
theme_set(theme_bw(base_size = 20))

bp_probs <- data.frame()
seqs <- character(2)

for(id in 1:2) {

 message(paste("Processing sequence", id))
 out <- paste("out", id, sep = "_")
 
 pfn <- paste(out, "bp_probabilities", sep = "/")
 sfn <- paste(out, "sequence.txt", sep = "/")
 
                   
 df <- read.table(file = pfn, col.names = c("start", "end", "p"))
 df$id <- id %>% as.factor
 bp_probs <- rbind(bp_probs, df)
 seqs[id] <- read.table(file = sfn)
}


ggplot(data = bp_probs, aes(x = p, fill = id)) +
  geom_density(alpha = 0.5)

ggplot(data = bp_probs, aes(x = id, y = p, fill = id)) +
  geom_violin(alpha = 0.5) +
  labs(y = "Probability", x = "Sequence id") +
  guides(fill = FALSE) +
  geom_hline(mapping = aes(yintercept = 0.8), linetype = 2)
ggsave("latex/figures/violin.png", height = 6, width = 10)


ggplot(data = bp_probs, aes(x = id, y = p, col = id)) +
  geom_jitter() +
  labs(y = "Probability", x = "Sequence id") +
  guides(col = FALSE) +
  geom_hline(mapping = aes(yintercept = 0.8), linetype = 2)
ggsave("latex/figures/jitter.png", height = 6, width = 10)

ggplot(data = bp_probs) +
  stat_ecdf(aes(x = p, col = id)) +
  labs(y = "Rank", x = "Probability") +
  scale_color_discrete(name = "Sequence",
                      breaks = c(1, 2),
                      labels = c("Seq1", "Seq2"))
ggsave("latex/figures/rank_plot.png", height = 6, width = 10)

g80 <- bp_probs %>% arrange(-p) %>% summarize(g80 = sum(p > 0.8))

message(g80)
  

