
# Project: Spatial Analysis of COVID-19 Mortality in Bariloche
# Manuscript title: Social inequalities and COVID-19 mortality between neighborhoods of Bariloche city, Argentina
# Manuscript author: Perner et al, 2023
# Code Author: Binod Acharya

#==============================================================================##
# Section 1: Spatial Adjacency stuffs 
#==============================================================================##

rm(list=ls())

library(rgdal)
library(tidyverse)
library(tmap)
library(readxl)
library(R2WinBUGS)
library(mitml)
library(rtf)
library(mcmcplots)
library(biscale)
library(cowplot)
library(sf)


# Import shapefile
shape_bar<-readOGR("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Data From Serena/mapa con  datos nse v1.shp")


# Adjacency matrix stuff 
shape_bar@data$orderid  <- as.numeric(rownames(shape_bar@data))+1
# I am adding snap parameter. Otherwise, one region would not have a neighbor.
tract_fid_link_bar     <- shape_bar@data[, c("orderid", "Link")]
my_adj                 <- spdep::poly2nb(shape_bar, queen = T,snap = 1/10^5) 
NumCells               <- length(my_adj)
num                    <- sapply(my_adj,length)
adj                    <- unlist(my_adj)
sumNumNeigh            <- length(unlist(my_adj))

# Check if any hole in adjacency structure etc.
summary(num) # Min number of adjacent areas is 1. Not zero.

#save.image("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Data From Serena/adjacency_stuff_bariloche.RData")



#==============================================================================##
# Section 2: Analytic Data Preparation
#==============================================================================##

rm(list=ls())

setwd("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Output")

# Bring in spatial stuffs, mortality, population stuffs
load("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Data From Serena/adjacency_stuff_bariloche.RData")

# Import mortality
mort<- read_delim("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Data From Serena/covid_bariloche.csv", delim=";")
mort$date_death<- as.Date(mort$date_death, format="%d/%m/%Y")

# Import population
pop <- read_csv("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/PopulationData/pop_data_bariloche.csv")

# Import standard population from WHO
standard_pop<-read_excel("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/PopulationData/WHO_StandardPop_2000-2025_v2.xlsx")


# Recode sex:
table(mort$sex, useNA="a")
mort$sex<- ifelse(mort$sex==1,"Females", "Males")

# check deaths that are not georeferenced: 10 deaths
mort[mort$radio_censal ==".",]

# check age that is missing: 1 death
mort[mort$age ==".",]

# Remove records with missing age and radio censales
mort<-mort %>% 
  mutate(L3LOCALID= as.numeric(radio_censal),
         age= as.numeric(age)) %>%
  drop_na(L3LOCALID, age) %>%
  mutate(age_100 = ifelse(age>=100, 100, age)) %>%
  select(-age, -radio_censal)

dim(mort) # number of rows (deaths) and columns in the data

# Check how many radio censales have no death at all.
setdiff(unique(mort$L3LOCALID),shape_bar@data$Link) # zero, all L3 ids in  death data are in the shapefile
setdiff(shape_bar@data$Link, unique(mort$L3LOCALID))# 9 L3s have not deaths at all.


#Work on cleaning standard population and population stuffs. For population clean up, we just need to subset for core Bariloche.

# scale the standard weight percent
standard_pop$wt <- standard_pop$avg2000_2025*100/sum(standard_pop$avg2000_2025)

# some wranglings..
standard_pop <- standard_pop %>%
  mutate(age_lowlim = str_extract(agecat, "[^-]+")) %>%
  mutate(age_lowlim = as.numeric(gsub("\\+","",age_lowlim))) %>%
  mutate(age= ifelse(age_lowlim >= 100, 100, age_lowlim)) %>%
  group_by(age) %>%
  summarize(wt=sum(wt))

# re-categorize standard population with 4 age groups
standard_pop_4 <- standard_pop %>%
 mutate(age4= case_when(age <= 35 ~ 1, # up to 39
                         age %in% c(40,45,50,55) ~ 2,
                         age %in% c(60,65,70,75) ~ 3,
                         age >=80 ~ 4)) %>%
  group_by(age4)%>%
  summarize(wt=sum(wt)/100)

# Subset population for core Bariloche geography
pop <- pop %>% filter(L3LOCALID %in% tract_fid_link_bar$Link)


# Now proceed in creating analytic dataset. Exposure population is currently approximated as 2 times the estimated population of 2018. We were only able to estimate population upto 2018 because of the lack of population projection dataset post 2018.

# Grid
template <- expand.grid(age_85= 0:85,
            sex = c("Females", "Males"),
            L3LOCALID = as.numeric(unique(tract_fid_link_bar$Link)))

# pop
pop_24month <- pop %>%
  select(SALID3, L3LOCALID, age_85, sex, pop18) %>%
  mutate(pop_24mo= pop18*2)

# mort
mort_24month <- mort %>%
  mutate(age_85 = ifelse(age_100 >=85, 85, age_100))%>%
  group_by(L3LOCALID,age_85,sex)%>%
  summarize(dth=n())%>%
  full_join(template, by = c("L3LOCALID", "age_85", "sex"))%>%
  mutate(dth=replace_na(dth,0))
  
#combine pop and mort
dat_24mo <- inner_join(mort_24month, pop_24month, by = c("L3LOCALID", "age_85", "sex"))


# Make the analytic data
ad4<- dat_24mo %>%
  mutate(age4= case_when(age_85 %in% 0:39 ~ 1,
                         age_85 %in% 40:59 ~ 2,
                         age_85 %in% 60:79 ~ 3,
                         age_85 >=80 ~ 4))%>%
  group_by(L3LOCALID,age4) %>%
  summarize(dth=sum(dth), pop= sum(pop_24mo))%>%
  mutate(pop=ifelse(pop==0,1,pop))
  

#==============================================================================##
# Section 3: Modeling part in WinBUGS
#==============================================================================##

# Fully Bayesian model in WinBUGS. 
# Have 4 age groups plus the spatial random effects for each L3s. 
# Put Besag CAR priors on spatial effects. Also have iid non-spatial component to capture heterogeneity for each area and age. 


# Save the model to use as text file
sink("bariloche_age4_car_v2.txt")        
cat("
model
{
  for(a in 1:N_age){ 
      for(i in 1:N_neighborhoods){ 
        Y[i,a] ~ dpois(nlambda[i,a])
        nlambda[i,a] <- n[i,a] * lambda[i,a]  
        lambda[i,a] <- exp(beta0[a] + w[i]+ nonspat[i,a])
        nonspat[i,a] ~ dnorm(0,tauns[a])
        
        lambda_who_wt[i,a] <- lambda[i,a] * who_wt[a]
        lambda_city_wt[i,a]<- lambda[i,a] * city_wt[a]
      
      }
      beta0[a] ~ dunif(-20,20)
      tauns[a] ~ dgamma(1, 0.01)
    }
  
      w[1:N_neighborhoods] ~ car.normal(adj[], weights[], num[], tau)
      tau~dgamma(1,0.142) 
      
  
  for(h in 1:sumNumNeigh){
    weights[h] <- 1
  }
  
  # obtain age-adjusted rate
   for(i in 1:N_neighborhoods){ 
   who_aa_rate[i] <- sum(lambda_who_wt[i,])*10000
   city_aa_rate[i]<- sum(lambda_city_wt[i,])*10000
   }
   
}
")
sink()


# Now prepare the analytic dataset in the form needed by WinBUGS.
dta<-ad4 %>% 
  inner_join(tract_fid_link_bar %>% 
               mutate(L3LOCALID= as.numeric(Link)),
             by = "L3LOCALID") %>%
  ungroup()%>% 
  group_by(age4, orderid, L3LOCALID) %>% 
  summarise(deaths=sum(dth),
            total=sum(pop)) %>% 
  mutate(rate=deaths/total) %>%
  arrange(age4, orderid) # this order is crucial


# Make the arrays  
N_neighborhoods = nrow(tract_fid_link_bar)
N_age = length(unique(dta$age4))
Y = array(dta$deaths, dim = c(N_neighborhoods,N_age))
n = array(dta$total, dim = c(N_neighborhoods, N_age))

who_wt <- standard_pop_4 %>% pull(wt)

city_wt <- dta %>% 
  group_by(age4) %>% 
  summarize(total=sum(total)) %>%
  mutate(wt=total/sum(total)) %>%
  pull(wt)


# Give initial values
inits <- function() {
  list(tau= 2,
       tauns=rep(2, 4),
       beta0=rep(-10, 4),
       nonspat=array(0, dim=c(N_neighborhoods,N_age)),
       w = array(0,dim=c(N_neighborhoods,N_age)))
}

# Define parameters to save
parameters= c("w","lambda","tau","tauns","beta0", "nonspat", "who_aa_rate", "city_aa_rate")

# Define iterations numbers, burn-ins, and thin factor
n.iter   <- 600000
n.burnin <- 500000
n.thin   <- 50


#Call WinBUGS and save the BUGS object.
bugs_v3 = bugs(model.file= "bariloche_age4_car_v2.txt",
               data = list(Y = Y, 
                           n = n,
                           N_neighborhoods=N_neighborhoods,
                           N_age= N_age,
                           who_wt=who_wt,
                           city_wt=city_wt,
                           sumNumNeigh=sumNumNeigh, num=num, adj=adj),
               inits= inits,
               parameters= parameters,
               n.chains = 2,
               n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin,
               debug=F, clearWD=TRUE,
               DIC= TRUE,
               working.directory = getwd(),
               bugs.directory = "C:/Users/ba525/Documents/WB/winbugs14_full_patched/WinBUGS14")

saveRDS(bugs_v3, file = 'bugs_bariloche_a4_v3.Rds')


#==============================================================================##
# Section 4: Post-WinBUGS processing
#==============================================================================##

bugs<- readRDS('//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Output/bugs_bariloche_a4_v3.Rds')

# Look at the trace plots 
traplot(bugs, "tau")   #precision of spatial random effect
traplot(bugs, "tauns") # precision of non-spatial random effect
traplot(bugs, "beta0") # coefficients for age effect


# Now work on the mortality rates. WHO age-adjusted rates are automatically saved in the BUGS object. 
# The "True" L3 level rates, not adjusted to external population but aggregated up from the age-specific rates using the L3 level population composition needs some extra work.


# Age specific rates
lambda_quant <- apply(bugs$sims.list$lambda, 2:3, quantile, probs= c(0.025, 0.5, 0.975))
lambda_quant <- data.table::as.data.table(lambda_quant)

lambda_quant <- lambda_quant %>% 
  setNames(c("quant", "orderid", "age4","rate")) %>%
  pivot_wider(names_from=quant, values_from=rate) %>% 
  setNames(c('orderid', 'age4', 'lci', 'median', 'uci'))

# True aggregated rates (not standardized to external pop)
lambda_df <- lambda_quant %>% 
  inner_join(dta, by = c("orderid", "age4")) %>% 
  mutate(Link=as.character(L3LOCALID))%>%
  inner_join(tract_fid_link_bar) %>%
  mutate(pop_final=total)%>%
  mutate(pop_rate = pop_final * median ,
         pop_rate_lci = pop_final * lci,
         pop_rate_uci = pop_final * uci) %>%
  group_by(L3LOCALID, Link)%>%
  summarise(numerator = sum(pop_rate),
            numerator_lci = sum(pop_rate_lci),
            numerator_uci = sum(pop_rate_uci),
            denom = sum(pop_final)) %>%
  mutate(rate_per10k_car = numerator *10000/denom,
         lci_per10k_car = numerator_lci *10000/denom,
         uci_per10k_car = numerator_uci *10000/denom)

# Data frame of WHO-age adjusted rate, pulled from BUGS
aa_who <- apply(bugs$sims.list$who_aa_rate,2, quantile,probs=c(0.025, 0.5, 0.975))%>%
  as.data.frame()%>%
  rownames_to_column()%>%
  pivot_longer(2:160)%>%
  mutate(orderid= as.numeric(substr(name,2,4))) %>%
  pivot_wider(id_cols=orderid, names_from=rowname, values_from=value)%>%
  setNames(c('orderid', 'who_aa_lci', 'who_aa_median', 'who_aa_uci'))


# Combine 2 sets of model-based rates: true, WHO adjusted
rates_from_model <-aa_who %>%
  inner_join(tract_fid_link_bar)%>%
  inner_join(lambda_df)%>%
  select(L3LOCALID, rate_per10k_car, lci_per10k_car, uci_per10k_car,
         who_aa_median, who_aa_lci, who_aa_uci)


# export data sets of age-adjusted rates
# write_csv(rates_from_model, "rates_from_model.csv")

# Convert this to long form that we might need later
rates_from_model_long<- rates_from_model %>%
  pivot_longer(!L3LOCALID)%>%
  transmute(L3LOCALID=L3LOCALID,
            value=value,
            method=case_when(grepl("car",name)~ "True",
                             grepl("who",name) ~"WHO adjusted"),
            
            quant=case_when(grepl("rate|median",name)~ "median",
                            grepl("lci",name) ~"lci",
                            grepl("uci", name) ~ "uci")) %>%
  pivot_wider(id_cols=c(L3LOCALID, method), names_from=quant, values_from=value)



# Make the map  of rates
shape2 <- shape_bar
shape2@data <- shape2@data %>% left_join(rates_from_model %>% mutate(Link=as.character(L3LOCALID)))


# Make the map for manuscript - WHO adjusted
# Map ultimately created in ArcGIS Pro.

map_who_manuscript<- tm_shape(shape2)+
  tm_layout(legend.show = T)+
  tm_polygons(col='who_aa_median' ,
              breaks= quantile(rates_from_model$who_aa_median, na.rm=T),
              title= "Rate per 10,000",palette = "seq",lwd=0.01)+   
  tm_layout(frame = FALSE,
            legend.title.size = 1.1,
            legend.text.size=.9)+
  tm_layout(aes.palette = list(seq = "-RdBu"),
            legend.position = c(0.05, -0.15))+
  tm_compass(type="4star", position=c(.7, .7), size=2.5, show.labels = 1)+
  tm_scale_bar(breaks = c(0,5,10), position=c(.7, 0), just=c(0,0.2))


cairo_ps(filename='Figure1_map_who_rate.eps', width=5, height=5)
print(map_who_manuscript)
dev.off()



# Make the map without age adjustment 
map_true_manuscript<- tm_shape(shape2)+
  tm_layout(legend.show = T)+
  tm_polygons(col='rate_per10k_car' ,
              breaks= quantile(rates_from_model$rate_per10k_car, na.rm=T),
              title= "Rate per 10,000",palette = "seq",lwd=0.01)+   
  tm_layout(
    frame = FALSE,
    main.title.size = 1.1,
    legend.text.size=.9)+
  tm_layout(aes.palette = list(seq = "-RdBu"),
            legend.position = c(0.05, -0.15))+
  tm_compass(type="4star", position=c(.7, .7), size=2.5, show.labels = 1)+
  tm_scale_bar(breaks = c(0,5,10), position=c(.7, 0), just=c(0,0.2))


cairo_ps(filename='eFigure1_map_true_modeled_rate.eps', width=5, height=5)
print(map_true_manuscript)
dev.off()


# Scatter-plot of unadjusted and adjusted rates for manuscript supplement
adjust_vs_crude_plt=ggplot(rates_from_model)+
  geom_point(aes(x=rate_per10k_car, y=who_aa_median), 
             pch=21, size=1.5,  fill="gray", color="black")+
  ylim(c(5,42))+
  xlim(c(5,42))+
  ylab("Age-adjusted mortality rate")+
  xlab("Crude mortality rate")+
  theme_bw()+
  theme(axis.title.y = element_text(size = 14, 
                                    margin = margin(t = 0, r = 15, b = 0, l = 0)))+
  theme(axis.title.x=element_text(size=14))


# Export plot
cairo_ps(filename='eFigure2_adjust_vs_crude_plt.eps', width=5, height=5)
plot(adjust_vs_crude_plt)
dev.off()


#==============================================================================##
# Section 5 : Association between mortality rate and SES/exposures
#==============================================================================##

#The L3-level covariates for study would be: 
# overcrowding (CNSCROWD3RML3)
# school attendance (CNSST1517L3), 
# education (CNSMINHS_L3),
# water inside home (CNSWATINL3), 
# sewage connection (CNSSEWANYL3),
# unemployment(CNSUNEMPL3)
# Unmet basic needs (NBI)

# Bring in co-variate data
covdata <- read_excel("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Data Request/2022_03_16/SEC_Census_L3_AR2010_08272020_localid_v2.xlsx")

#Unmet Basic Needs (Necesidades Básicas Insatisfechas): percentage of households with at least one unmet #basic need
nbi<-read_delim("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Data From Serena/nbi_l3_bche.csv",delim=";", locale=locale(decimal_mark = ","))
names(nbi) <- c("L3LOCALID", "NBI")
nbi$L3LOCALID <- as.character(nbi$L3LOCALID)

# subset the covariates data for Bariloche L3s
covdata <- covdata %>% filter(L3LOCALID %in% shape_bar@data$Link)%>% inner_join(nbi)

# Specify the potential covariates
need_vars <- c('CNSCROWD3RML3','CNSST1517L3','CNSMINHS_L3',
               'CNSWATINL3','CNSSEWANYL3','CNSSEWNETL3','CNSUNEMPL3', 'NBI')

# Pair-wise Pearson correlation
round(cor(covdata[all_of(need_vars)]),1)

# Only keep selected variables / no parameterization
covdata_sel <- covdata %>% 
  select(L3LOCALID, all_of(need_vars)) %>%
  mutate(L3LOCALID= as.numeric(L3LOCALID))

#write_csv(covdata_sel, "covdata_sel.csv")



# no parameterization, reversed
covdata_sel_reversed <- covdata_sel 
to_reverse= c("CNSMINHS_L3", "CNSST1517L3", "CNSWATINL3", "CNSSEWANYL3", "CNSSEWNETL3")
to_reverse_vec=which(colnames(covdata_sel_reversed) %in% to_reverse)
covdata_sel_reversed[,to_reverse_vec] <- (covdata_sel_reversed[,to_reverse_vec])*-1

# decile-based parm and reversed
covdata_sel_minmaxparm_decile_reversed<-covdata_sel_reversed
covdata_sel_minmaxparm_decile_reversed[,2:9] <-data.frame(sapply(covdata_sel_reversed[,2:9], ntile,10))
covdata_sel_minmaxparm_decile_reversed[,2:9] <-(covdata_sel_minmaxparm_decile_reversed[,2:9] -1)/9

# Combine mortality rates and covariate
rate_with_cov <- rates_from_model %>% inner_join(covdata_sel, by = "L3LOCALID")

# Needed covs df
need <- data.frame(name=c('NBI',
                          'CNSCROWD3RML3',
                          'CNSMINHS_L3',
                          'CNSST1517L3',
                          'CNSUNEMPL3',
                          'CNSWATINL3',
                          'CNSSEWANYL3'),
                   label=c('Unmet basic needs',
                           'Overcrowding',
                           'Completed high school',
                           'School attendance',
                           'Unemployment',
                           'Piped water inside dwelling',
                           'Sewage connection'))

# Make the map of SES (univariable maps)
# These maps are ultimately replaced with bi-variate maps during revision 
map_ses=tm_shape(shape3)+
  tm_polygons(col= c( "NBI", 
                      "CNSCROWD3RML3",
                      "CNSMINHS_L3", 
                      "CNSST1517L3",
                      "CNSUNEMPL3",
                      "CNSWATINL3",
                      "CNSSEWANYL3"),
              title="",
              lwd=0.01,
              border.alpha = 0.0,
              ncol=2,
              breaks= list(
                quantile(covdata_sel$NBI),
                quantile(covdata_sel$CNSCROWD3RML3),
                quantile(covdata_sel$CNSMINHS_L3),
                quantile(covdata_sel$CNSST1517L3),
                quantile(covdata_sel$CNSUNEMPL3),
                quantile(covdata_sel$CNSWATINL3),
                quantile(covdata_sel$CNSSEWANYL3)),
              palette=list("Oranges", "Oranges", "Purples", "Purples", "Reds", "Purples", "Purples"))+
  
  tm_layout(title.position = c("right", "top"),
            legend.position=  c(0.05, 0),
            legend.text.size = 1, 
            frame = T,
            panel.show=TRUE,
            panel.labels=c("Unmet basic needs",
                           "Overcrowding",
                           "Completed high school",
                           "School attendance", 
                           "Unemployment", 
                           "Piped water inside dwelling",
                           "Sewage connection"),
            panel.label.size=1.3,
            asp=0)


cairo_ps(filename='eFigure3_map_ses.eps', width=8.5, height=11)
print(map_ses)
dev.off()


# Create scatter plot between mortality rate and SES
plt_ass= rate_with_cov %>%
  select(L3LOCALID, who_aa_median, all_of(need_vars))%>%
  pivot_longer(cols= !c('L3LOCALID', 'who_aa_median')) %>%
  inner_join(need)%>%
  mutate(label=factor(label, levels=need$label))%>%
  ggplot(aes(x=value, y= who_aa_median))+
  geom_point(pch=21)+
  facet_wrap(~label, scales="free_x", ncol=4)+
  #scale_y_continuous(trans="log10", breaks=c(10,20,30,40,50))+
  scale_y_continuous(breaks=c(10,20,30,40,50))+
  theme_bw()+
  geom_smooth(method="lm", se=F, color="gray")+
  labs(y="Age-adjusted mortality rate per 10,000", x="")+
  theme(panel.spacing = unit(.75, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside")+
  theme(axis.title.y = element_text(size = 14, color="black",
                                    margin = margin(t = 0, r = 15, b = 0, l = 0)))+
  theme(axis.title.x=element_text(size=13, color="black"))+
  theme(axis.text=element_text(size=12, color="black"))+
  theme(strip.text = element_text(size = 13, color="black"))


# Print and export
cairo_ps(filename='Figure2_bariloche_plt_assv2.eps', width=10, height=8)
plot(plt_ass)
dev.off()



#Remove variable that are correlated with each other above certain threshold
tmp <- cor(covdata[all_of(need_vars)])
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0
data.old= covdata[all_of(need_vars)]
data.new <- data.old[, !apply(tmp, 2, function(x) any(abs(x) > 0.6, na.rm = TRUE))]
head(data.new)


# Data frame for association modeling; incorporates samples of posterior distribution.
# Reversed
aa_who_nits_minmaxparm_decile_reversed <-as.data.frame(bugs$sims.list$who_aa_rate)%>%
  mutate(nits=1:n())%>%
  pivot_longer(1:159)%>%
  mutate(orderid= as.numeric(substr(name,2,4))) %>%
  select(-name)%>%
  inner_join(tract_fid_link_bar, by = "orderid")%>%
  transmute(nits=nits, L3LOCALID= as.numeric(Link), mort_rate_per10k= value) %>%
  inner_join(covdata_sel_minmaxparm_decile_reversed,by = "L3LOCALID")

#write_csv(aa_who_nits_minmaxparm_decile_reversed, "aa_who_nitsv2_minmaxparm_decile_reversed.csv")



### (A) Run model with log(rate) as outcome and some predictors are reversed. Gives RII in min-max parameterization upon exponentiation of coefficients.

aa_who_nits_minmaxparm_decile_reversed$log_mort = log(aa_who_nits_minmaxparm_decile_reversed$mort_rate_per10k/10000)
aa_who_nits_minmaxparm_list_decile_reversed  <- aa_who_nits_minmaxparm_decile_reversed %>% group_split(nits)
aa_who_nits_minmaxparm_list2_decile_reversed <- as.mitml.list(aa_who_nits_minmaxparm_list_decile_reversed)


#1. with edu
dfit.m1 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(log_mort ~ CNSMINHS_L3))
dfit.m1_ = testEstimates(dfit.m1)


#2. with school attendance
dfit.m2 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(log_mort ~ CNSST1517L3))
dfit.m2_ = testEstimates(dfit.m2)

#3. with overcrowding
dfit.m3 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(log_mort ~ CNSCROWD3RML3))
dfit.m3_ = testEstimates(dfit.m3)

#4. with unemployment
dfit.m4 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(log_mort ~ CNSUNEMPL3))
dfit.m4_ = testEstimates(dfit.m4)

#5. with NBI
dfit.m5 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(log_mort ~ NBI))
dfit.m5_ = testEstimates(dfit.m5)


#6. with sewage
dfit.m6 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(log_mort ~ CNSSEWANYL3))
dfit.m6_ = testEstimates(dfit.m6)


#7. water inside dwelling
dfit.m7 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(log_mort ~ CNSWATINL3))
dfit.m7_ = testEstimates(dfit.m7)


# Extract coefficients from decile based min max parameterized model
var=rbind(rownames(dfit.m1_$estimates)[2],
          rownames(dfit.m2_$estimates)[2],
          rownames(dfit.m3_$estimates)[2],
          rownames(dfit.m4_$estimates)[2],
          rownames(dfit.m5_$estimates)[2],
          rownames(dfit.m6_$estimates)[2],
          rownames(dfit.m7_$estimates)[2])

val=rbind(dfit.m1_$estimates[2, c(1,2)],
          dfit.m2_$estimates[2, c(1,2)],
          dfit.m3_$estimates[2, c(1,2)],
          dfit.m4_$estimates[2, c(1,2)],
          dfit.m5_$estimates[2, c(1,2)],
          dfit.m6_$estimates[2, c(1,2)],
          dfit.m7_$estimates[2, c(1,2)])

tab_d=cbind(data.frame(var),val) %>%
  mutate(d_beta_ci = paste0(round(Estimate,3),
                            "(",
                            round(Estimate - 1.96*Std.Error, 3),
                            ",",
                            round(Estimate + 1.96*Std.Error, 3),
                            ")")) %>%
  mutate(d_exp_beta_ci= paste0(round(exp(Estimate),3),
                               "(",
                               round(exp(Estimate - 1.96*Std.Error), 3),
                               ",",
                               round(exp(Estimate + 1.96*Std.Error), 3),
                               ")")) %>%
  select(-Estimate,-Std.Error)


min= sapply(covdata_sel_reversed[,2:9], min)
max= sapply(covdata_sel_reversed[,2:9], max)

# Combine estimates from all model forms
outtab_reversed = tab_d %>% 
  inner_join(data.frame(min) %>% rownames_to_column("var")) %>%
  inner_join(data.frame(max) %>% rownames_to_column("var")) 

sddf=covdata_sel %>%
  select(all_of(outtab_reversed$var))%>%
  summarize_all(sd) %>%
  pivot_longer(cols=1:7,names_to="var", values_to ="SD")

logmodel_reversed_outtab<- outtab_reversed%>%
  inner_join(sddf)%>%
  select(var, min, max, SD, everything())

#write.csv(logmodel_reversed_outtab, "logmodel_coefficients_reversed.csv", row.names = F)



##B. Run model with rate per 10,000 as outcome. Some Predictors are reversed. Gives SII in min-max parameterization.


#1. with edu
dfit.m1 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(mort_rate_per10k ~ CNSMINHS_L3))
dfit.m1_ = testEstimates(dfit.m1)


#2. with school attendance
dfit.m2 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(mort_rate_per10k ~ CNSST1517L3))
dfit.m2_ = testEstimates(dfit.m2)


#3. with overcrowding
dfit.m3 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(mort_rate_per10k ~ CNSCROWD3RML3))
dfit.m3_ = testEstimates(dfit.m3)


#4. with unemployment
dfit.m4 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(mort_rate_per10k ~ CNSUNEMPL3))
dfit.m4_ = testEstimates(dfit.m4)


#5. with NBI
dfit.m5 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(mort_rate_per10k ~ NBI))
dfit.m5_ = testEstimates(dfit.m5)


#6. with sewage
dfit.m6 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(mort_rate_per10k ~ CNSSEWANYL3))
dfit.m6_ = testEstimates(dfit.m6)


#7. water inside dwelling
dfit.m7 <- with(aa_who_nits_minmaxparm_list2_decile_reversed, lm(mort_rate_per10k ~ CNSWATINL3))
dfit.m7_ = testEstimates(dfit.m7)


# Extract coefficients from decile based min max parameterized model
var=rbind(rownames(dfit.m1_$estimates)[2],
          rownames(dfit.m2_$estimates)[2],
          rownames(dfit.m3_$estimates)[2],
          rownames(dfit.m4_$estimates)[2],
          rownames(dfit.m5_$estimates)[2],
          rownames(dfit.m6_$estimates)[2],
          rownames(dfit.m7_$estimates)[2])

val=rbind(dfit.m1_$estimates[2, c(1,2)],
          dfit.m2_$estimates[2, c(1,2)],
          dfit.m3_$estimates[2, c(1,2)],
          dfit.m4_$estimates[2, c(1,2)],
          dfit.m5_$estimates[2, c(1,2)],
          dfit.m6_$estimates[2, c(1,2)],
          dfit.m7_$estimates[2, c(1,2)])

tab_d=cbind(data.frame(var),val) %>%
  mutate(d_beta_ci = paste0(round(Estimate,3),
                            "(",
                            round(Estimate - 1.96*Std.Error, 3),
                            ",",
                            round(Estimate + 1.96*Std.Error, 3),
                            ")")) %>%
  select(-Estimate,-Std.Error)

min= sapply(covdata_sel_reversed[,2:9], min)
max= sapply(covdata_sel_reversed[,2:9], max)

# Combine estimates from all model forms
outtab_reversed = tab_d %>%
  inner_join(data.frame(min) %>% rownames_to_column("var")) %>%
  inner_join(data.frame(max) %>% rownames_to_column("var")) 

nologmodel_reversed_outtab<- outtab_reversed%>%
  inner_join(sddf)%>%
  select(var, min, max, SD, everything())

#write.csv(nologmodel_reversed_outtab, "linear_model_coefficients_reversed.csv", row.names = F)



#Clean-up the tables a little bit--based on the some variables reversed 
var_label=data.frame(
  var=c("NBI",
        "CNSCROWD3RML3",
        "CNSMINHS_L3",
        "CNSST1517L3",
        "CNSUNEMPL3",
        "CNSWATINL3",
        "CNSSEWANYL3"),
  label=c(
    "% of population with at least one unmeet basic needs",
    "% of households with overcrowding",
    "% of population over 25 years who completed high school",
    "% of population aged 15-17 attending school",
    "Unemployment", 
    "% of households with piped water access inside the dwelling",
    "% of households connected to a sewage system of any type"))
var_label$label= factor(var_label$label, levels=unique(var_label$label))


table_reversed = nologmodel_reversed_outtab %>%
  transmute(type="4 var reversed",
            var=var, 
            sii_decile= d_beta_ci) %>%
  inner_join(
    (logmodel_reversed_outtab %>%
       transmute(var=var, 
                 rii_decile= d_exp_beta_ci))) %>%
  inner_join(var_label) %>%
  select(-var) %>%
  select(label, everything()) %>%
  arrange(label)

# Save regression results
rtffile <- RTF("reg_result_table_reversed.doc") 
addTable(rtffile, table_reversed)
done(rtffile)






#==============================================================================##
# Bivariate maps of age-adjusted mortality and SES
#==============================================================================##

## Bi-variate maps of age-adjusted mortality rate and SES
## These maps are ultimately created in ArcGIS Pro but as an example R code is provided below.
shape3 <- read_sf("//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL/Manuscripts/MS191_Perner/Data From Serena/mapa con  datos nse v1.shp")

shape3 <- shape3 %>%
  left_join(covdata_sel, by=c("Link"="L3LOCALID")) %>%
  left_join(rates_from_model,by=c("Link"="L3LOCALID"))

# create classes
data <- bi_class(shape3, y = CNSMINHS_L3, x = who_aa_median, style = "quantile", dim = 3)
map <- ggplot() +
  geom_sf(data = shape3, mapping = aes(fill = data$bi_class), color = "white", size = 0.1, 
          show.legend = FALSE) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  bi_theme()
legend <- bi_legend(pal = "DkViolet", dim = 3, xlab = "Mortality", ylab = "Education", size = 10)

finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.2, 0.65, 0.2, 0.2)

print(finalPlot)





