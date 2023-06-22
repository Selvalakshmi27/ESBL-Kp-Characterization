library("ggpubr")

CG307 <- c(1,
           8,
           9,
           8,
           6,
           2,
           4)
CG307_rates <- c(7.692307692,
                 34.7826087,
                 47.36842105,
                 33.33333333,
                 28.57142857,
                 9.523809524,
                 40
)

rare <- c(7,
          12,
          5,
          11,
          10,
          8,
          4)
rare_rates <- c(53.84615385,
                52.17391304,
                26.31578947,
                45.83333333,
                47.61904762,
                38.0952381,
                40
)
ESBL_numbers <- c(16,
                24,
                25,
                23,
                20,
                19,
                12)

ESBL_rates <- c(19.51219512,
                22.64150943,
                26.04166667,
                29.48717949,
                33.89830508,
                21.34831461,
                25.6666
)

non_ESBL_rates <- c(80.48780488,
                    77.35849057,
                    73.95833333,
                    70.51282051,
                    66.10169492,
                    78.65168539,
                    75.55555556
)

year <- c(2016,2017,2018,2019,2020,2021,2022)

df <- data.frame(CG307,ESBL_numbers,year)
df
ggscatter(df, x = "CG307", y = "ESBL_numbers",
          color="black", size=3,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",cor.coef.size = 10,cor.coef.coord =c(1.5,27),
          add.params = list(color = "blue", fill = "lightgray"),
          xlab = "Number of CG307s per year", ylab ="Number of ESBL-Kp per year", label = df$year) + font("xlab", size = 22, face="italic") + font("ylab", size = 22, face ="italic")+
  theme(axis.text = element_text(size = 23))  

ggsave("correlation_CG307.png",width = 10.36, height = 8.14, dpi = 900)
