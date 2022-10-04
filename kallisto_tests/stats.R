library(tidyverse)
library(kableExtra)
library(faraway)
library(broom)

top = c(read.csv("final_res_top100.tsv", sep = "\t")['tpm'])
TP = c(read.csv("final_res_TP100.tsv", sep = "\t")['tpm'])


dully <-data.frame(tpm =c(TP),
                   Results =rep(c('100 True Positives'), each = 97))

gully <-data.frame(tpm =c(top),
                   Results =rep(c('Top 100 MinCE'), each = 97))

total <- rbind(gully, dully)

total$Results = factor(total$Results, levels = c("Top 100 MinCE",
                                               "100 True Positives"))
ggplot(total, aes(x=Results, y=tpm)) + geom_boxplot()

group_by(total, Results) %>%
  summarise(
    min = min(tpm, na.rm = TRUE),
    med = median(tpm, na.rm = TRUE),
    mean = mean(tpm, na.rm = TRUE),
    max = max(tpm, na.rm = TRUE),
    sd = sd(tpm, na.rm = TRUE) ) %>% kbl(align = 'c') %>%
  kable_styling(latex_options = "HOLD_position")
