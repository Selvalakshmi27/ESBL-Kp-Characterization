library(hrbrthemes)
library(tidyverse)
library(ggplot2)

#data
df<- read.csv("recurrence_plot.csv")
attach(df)
df$CRE
#for recurrent isolates
ggplot(df) +
    geom_segment(aes(x = 0, xend = max(Days), y = Patient, yend = Patient),
                 linetype = "solid", size = 0.5, color="black") +
  geom_point(aes(x = Days, y = Patient, color = as.factor(CRE)),size = 7.5) + theme_bw() +
  scale_color_manual(values = c("ST307" = "#61a5c2",
                                "ST29"="#31572c",
                                "ST258"="#000000",
                                "ST395"="#fcb9b2",
                                "Other" = "#c8553d",
                                "non-CP-CRE" = "red",
                                "CP-CRE" = "orange")) +
  scale_y_continuous(breaks = seq(from = 1, to = 21, by = 1)) +
  labs(
    x = "Number of Days", y = "Patient",
    title = "Recurrent Isolates",
    color="ST",size=1.5
  ) +
  theme(legend.position = "top") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=16),
                                        axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16),
                                        legend.text = element_text(size = 16))


#for recurrent isolates with non-CP-CRE
ggplot(df) +
  geom_segment(aes(x = 0, xend = max(Days), y = Patient, yend = Patient),
               linetype = "solid", size = 0.5, color="black") +
  geom_point(aes(x = Days, y = Patient, color = as.factor(ST)),size = 5.5) + theme_bw() +
  scale_color_manual(values = c("ST307" = "#023e8a",
                                "ST29"="#c1121f",
                                "ST258"="#000000",
                                "ST395"="#fcb9b2",
                                "Other" = "grey")) +
  scale_y_continuous(breaks = seq(from = 1, to = 27, by = 1)) +
  scale_color_ipsum(name = "ST") +
  labs(
    x = "Number of Days", y = "Patient",
    title = "Recurrent Isolates",size=1.4
  ) +
  theme(legend.position = "top")   





ggsave("test.png",width = 9.36, height = 8.14,dpi=900)


               