rm(list = ls())


# Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,
  rio,
  stringr,
  compareGroups,
  hrbrthemes,
  lmtest,
  sandwich,
  missForest,
  stats,
  MASS,
  gdata,
  Matrix,
  cowplot,
  genscore,
  plyr,
  ggplot2)


set.seed(1221)

### load data 
#### change here to path of your data

path <- "/Users/kgoebler/Desktop/Lektorat/Sabrinas Doktorarbeit"  

data <- import(paste(path, "Tabelle ParoVeg_Imputation_kurz.xlsx", sep = "/"),
               which = 1, na = "NA")

arcsinevars <- colnames(dplyr::select(data, contains("arcsin")))
renamevars <- colnames(dplyr::select(data, contains("...")))
dropvars <- c("id", renamevars, arcsinevars)

data.raw <- dplyr::select(data, -all_of(dropvars)) 


changevars <- c("crp_t01", "crp_t02", "crp_t03", "holo_t02", "ins_t03",
                "sf_t01", "sf_t02", "sf_t03")

data.raw <- data.raw %>% 
  mutate(score = ifelse(
            SCORE == "niedrig",1, ifelse(
            str_detect(SCORE, "mo"), 2, ifelse(
            SCORE == "hoch", 3, NA))),
         crp_t01 = as.numeric(ifelse(crp_t01 == "<0,3", 0, crp_t01)),
         crp_t02 = as.numeric(ifelse(crp_t02 == "<0,3", 0, crp_t02)),
         crp_t03 = as.numeric(ifelse(crp_t03 == "<0,3", 0, crp_t03)),
         ins_t03 = as.numeric(ifelse(ins_t03 == "folgt", NA, ins_t03)),
         sf_t01 = as.numeric(ifelse(sf_t01 == "<20", 0, sf_t01)),
         sf_t02 = as.numeric(ifelse(sf_t02 == "<20", 0, sf_t02)),
         sf_t03 = as.numeric(ifelse(sf_t03 == "<20", 0, sf_t03)),
         )

data.raw <- dplyr::select(data.raw, -all_of(c("SCORE", "holo_t02", "hy")))


data.final <- filter(data.raw, !(is.na(t02) & is.na(t03)))

##### imputation #### 

misstoohigh <- which(sapply(1:ncol(data.final), function(i) sum(is.na(data.final[,i]))/nrow(data.final)) > .2)

data.imp.pre <- dplyr::select(data.final, -all_of(c("t01", "t02","t03", colnames(data.final[,misstoohigh]))))
data.imputed <- missForest(data.imp.pre)

data.complete <- cbind(data.final[,c("t01", "t02","t03")], data.imputed$ximp)


#### Descriptive Analysis ####

desc.data <- data.complete

p1 <- desc.data %>%
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control"))) %>%
  ggplot(aes(x=age, fill=group)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', bins = 10) +
  labs(fill="") +
  facet_grid(~group)



p2 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control"))) %>% 
  ggplot(aes(sexe-1, group = as.factor(group))) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "Gender", labels = c("male", "female")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

p3 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         medics = ifelse(m == 2,0,1)) %>% 
  ggplot(aes(medics, group = as.factor(group))) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "presricption drugs", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)


p4 <- desc.data %>%
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control"))) %>%
  ggplot(aes(x=py, fill=group)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', bins = 10) +
  labs(fill="", x = "Pack years") +
  facet_grid(~group)


p5 <- desc.data %>%
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control"))) %>%
  ggplot(aes(x=nu, fill=group)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', bins = 10) +
  labs(fill="", x = "Number of Teeth") +
  facet_grid(~group)


p6 <- desc.data %>%
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control"))) %>%
  ggplot(aes(x=bmi_t01, fill=group)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', bins = 10) +
  labs(fill="", x = "Body Mass Index") +
  facet_grid(~group)


p7 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         therapy = ifelse(th == 2,0,1)) %>% 
  ggplot(aes(therapy, group = as.factor(group))) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "antihypertensive therapy", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)


p8 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         acohol = ifelse(th == 2,0,1)) %>% 
  ggplot(aes(acohol, group = as.factor(group))) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "regular drinking", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

p9 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         metsyndrom = ifelse(th == 2,0,1)) %>% 
  ggplot(aes(metsyndrom, group = as.factor(group))) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "Regular Drinking", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

p10 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         metsyndrom = ifelse(th == 2,0,1)) %>% 
  ggplot(aes(metsyndrom, group = as.factor(group))) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "Metabolic Syndrom", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

p11 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         medzu = ifelse(medzu == 2,0,1)) %>% 
  ggplot(aes(medzu, group = as.factor(group))) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "Change in BP medication additionally", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

p12 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         medred = ifelse(medred == 2,0,1)) %>% 
  ggplot(aes(medred, group = as.factor(group))) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "Change in BP medication reduced", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

### could potentially remove those two variables

p13 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         fa = ifelse(fa == 2,0,1)) %>% 
  ggplot(aes(fa, group = group)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "Family History", labels = c("negative", "positive")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

p14 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         fa = ifelse(path == 2,0,1)) %>% 
  ggplot(aes(path, group = group)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "PATherapy w.o. year", labels = c("negative", "positive")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

p15 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         zvp = ifelse(zvp == 2,0,1)) %>% 
  ggplot(aes(zvp, group = group)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "Loss of Teeth (paradontologic)", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)

p16 <- desc.data %>% 
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control")),
         muhy1 = ifelse(muhy1 == 2,0,1)) %>% 
  ggplot(aes(muhy1, group = group)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", alpha = .6) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  scale_fill_discrete(name = "", labels = c("no", "yes")) +
  labs(x="") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~group)


p16 <- desc.data %>%
  mutate(group = factor(g, levels = c(1,2), labels = c("Treat", "Control"))) %>%
  ggplot(aes(x=muhy1, fill=group)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', bins = 3) +
  labs(fill="", x = "Oral hygiene (# per day)") +
  facet_grid(~group)


desc.data <- desc.data %>% 
  mutate(group = factor(desc.data$g, levels = c(1,2), labels = c("Treat", "Control")),
         sex = factor(desc.data$sexe, levels = c(1,2), labels = c("Male", "Female")))

table <- compareGroups(group ~ sex + age + py + zig + bop_t01 + bop_t02 + bop_t03 + 
                              pcr_t01 + pcr_t02 + pcr_t03 + st_t01 + st_t02 + st_t03 + 
                              hb1_t01 + hb1_t02 + hb1_t03 + glu_t01 + glu_t03 + glu_t02 + 
                              bmi_t01 + bmi_t02 + bmi_t03 + whr_t01 + whr_t02 + whr_t03 + 
                              `pisa-score_t01` + `pisa-score_t02` + `pisa-score_t03`, 
                              data = desc.data, method = NA)
descr.table <- createTable(table)
export2word(descr.table, file='table1.docx')
#export2xls(restab, file='table1.xlsx')
#export2csv(restab, file='table1.csv')
### even more detailed

summary(table)

#### make nice boxplots of the two significant variables

bopdata <- as.data.frame(cbind(
                 c(desc.data$bop_t01,desc.data$bop_t02,desc.data$bop_t03),
                 rep(factor(desc.data$g, levels = c(1,2), labels = c("Treat", "Control")),3),
                 factor(c(rep(1,length(desc.data$bop_t01)), rep(2,length(desc.data$bop_t01)), rep(3,length(desc.data$bop_t01))),
                          levels = c(1,2,3), labels = c("Baseline", "T1", "T2"))
                 ))

bopdata <- bopdata %>% 
  mutate(BOP = V1,
         Group = factor(V2, levels = c(1,2), labels = c("Treat", "Control")),
         Time = factor(V3, levels = c(1,2,3), labels = c("Baseline", "T1", "T2"))
  )

p17 <- ggplot(bopdata, aes(x=Time, y=BOP, fill=Group)) + 
  geom_boxplot() +
  facet_wrap(~Group)

##### pisa #####

pisadata <- as.data.frame(cbind(
  c(desc.data$`pisa-score_t01`,desc.data$`pisa-score_t02`,desc.data$`pisa-score_t03`),
  rep(factor(desc.data$g, levels = c(1,2), labels = c("Treat", "Control")),3),
  factor(c(rep(1,length(desc.data$`pisa-score_t01`)), rep(2,length(desc.data$`pisa-score_t01`)), rep(3,length(desc.data$`pisa-score_t01`))),
         levels = c(1,2,3), labels = c("Baseline", "T1", "T2"))
))
pisadata <- pisadata %>% 
  mutate(PisaScore = V1,
         Group = factor(V2, levels = c(1,2), labels = c("Treat", "Control")),
         Time = factor(V3, levels = c(1,2,3), labels = c("Baseline", "T1", "T2"))
  )

p18 <- ggplot(pisadata, aes(x=Time, y=PisaScore, fill=Group)) + 
  geom_boxplot() +
  facet_wrap(~Group)



#### Regression Analysis ####
reg.data <- data.complete
reg.data <- reg.data %>% 
  mutate(group = factor(reg.data$g, levels = c(2,1), labels = c("Control", "Treat")),
         sex = factor(reg.data$sexe, levels = c(1,2), labels = c("Male", "Female")))
reg.data <- dplyr::rename(reg.data, pisa_score_t01 = `pisa-score_t01`,
       pisa_score_t02 = `pisa-score_t02`,
       pisa_score_t03 = `pisa-score_t03`)


frml1 <- as.formula(bop_t01 ~ group + sex + age + py + zig + 
          pcr_t01 + st_t01 +
          hb1_t01 + glu_t01 + 
          bmi_t01 + whr_t01 + 
          pisa_score_t01)

frml2 <- as.formula(bop_t02 ~ group + sex + age + py + zig + 
                      bop_t01 + pcr_t01 + st_t01 +
                      hb1_t01 + glu_t01 + 
                      bmi_t01 + whr_t01 + 
                      pisa_score_t01 + pcr_t02 + st_t02 +
                      hb1_t02 + glu_t02 + 
                      bmi_t02 + whr_t02 + 
                      pisa_score_t02)

frml3 <- as.formula(bop_t03 ~ group + sex + age + py + zig + 
                      bop_t02 + bop_t01 + pcr_t01 + st_t01 +
                      hb1_t01 + glu_t01 + 
                      bmi_t01 + whr_t01 + 
                      pisa_score_t01 + pcr_t02 + st_t02 +
                      hb1_t02 + glu_t02 + 
                      bmi_t02 + whr_t02 + 
                      pisa_score_t02 + pcr_t03 + st_t03 +
                      hb1_t03 + glu_t03 + 
                      bmi_t03 + whr_t03 + 
                      pisa_score_t03)

lm_mod1 <- lm(frml1, data = reg.data)
summary(lm_mod1)

stargazer(lm_mod1, out = "reg1.txt", type = "text")

lm_mod2 <- lm(frml2, data = reg.data)
summary(lm_mod2)

stargazer(lm_mod2, out = "reg2.txt", type = "text")

lm_mod3 <- lm(frml3, data = reg.data)
summary(lm_mod3)


stargazer(lm_mod3, out = "reg3.txt", type = "text")


frml.new1 <- as.formula(bop_t01 ~ group + age + py + nu + whr_t01 + 
                          pcr_t01 + sf_t01 + hb1c_t01 + tag_t01 + crp_t01 +aasi_t01)

lm_mod.new1 <- lm(frml.new1, data = reg.data)
summary(lm_mod.new1)
coeftest(lm_mod.new1, vcov = vcovHC(lm_mod.new1, type = "HC5")) 

#frml.new2 <- as.formula(bop_t02 ~ group + age + py + nu + whr_t01 + whr_t02+ 
#                          pcr_t01 + pcr_t02 + sf_t01 + sf_t02 + hb1c_t01 + hb1c_t02 +
#                          tag_t01 + tag_t02 + crp_t01 + crp_t02 +
#                          aasi_t01 + aasi_t02)
frml.new2 <- as.formula(bop_t02 ~ group + age + py + nu + 
                          whr_t02 + 
                          pcr_t02 + 
                          sf_t02 + 
                          hb1c_t02 +
                          tag_t02 +
                          crp_t02 + 
                          aasi_t02)

lm_mod.new2 <- lm(frml.new2, data = reg.data)
summary(lm_mod.new2)
# robust version due to likely present heteroskedasticity
coeftest(lm_mod.new2, vcov = vcovHC(lm_mod.new2, type = "HC5")) 


#frml.new3 <- as.formula(bop_t03 ~ group + age + py + #nu + 
#                          whr_t01 + whr_t02 + whr_t03 + 
#                          pcr_t01 + pcr_t02 + pcr_t03+ 
#                          sf_t01 + sf_t02 + sf_t03 + 
#                          hb1c_t01 + hb1c_t02 + hb1c_t03 +
#                          tag_t01 + tag_t02 + tag_t03 +
#                          crp_t01 + crp_t02 + crp_t03 + 
#                          aasi_t01 + aasi_t02 + aasi_t03)

frml.new3 <- as.formula(bop_t03 ~ group + age + py  + nu +
                          whr_t03 + 
                          pcr_t03 + 
                          sf_t03 + 
                          hb1c_t03 +
                          tag_t03 +
                          crp_t03 + 
                          aasi_t03)

lm_mod.new3 <- lm(frml.new3, data = reg.data)
summary(lm_mod.new3)
# robust version due to likely present heteroskedasticity and leverage points
coeftest(lm_mod.new3, vcov = vcovHC(lm_mod.new3, type = "HC5")) 


### turns out: nu alters treatment effect insignificant and smaller but still consistently negative
### couldn't use bu (Bauchumfang) as this is not in the dataset. I used the whr instead.



