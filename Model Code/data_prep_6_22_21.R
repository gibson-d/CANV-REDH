require(sp)
library(RColorBrewer)
library(maps)
library(rgdal)
library(mapproj)
library(rgeos)
library(maptools)
library(raster)
library(dplyr)
library(readr)
library(ggplot2)
library(dplyr)
library(reshape)
library(Imap)
library(sf)
library(spdep)
library(abind)

file1 <- "https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data/releases.csv"
raw.band<-read_csv(file1)  #reading in CSV
tail(raw.band) #Unknown number records
dim(raw.band)  #71193 records 43 variables
summary(raw.band)

########################################################################
#cleaning data
########################################################################
#only use status 3 birds
clean.re_ca <-subset(raw.band,Status==3)                #69,363 note:records refer to "count of birds" not n banded
clean.re_ca$Count.of.Birds <- as.numeric(clean.re_ca$Count.of.Birds)

#only use B.Year from 1951 onwards
clean.re_ca<-subset(clean.re_ca,B.Year>=1951)          #65,199 bandings 1951-present

#remove unknown sex
clean.re_ca<-subset(clean.re_ca,Sex!=0)                #64,649 records (418 excluded)

#remove sex determined subsequent
clean.re_ca<-subset(clean.re_ca,Sex!=6)                #64,644 records   #!remove (5) males from subs encounter
clean.re_ca<-subset(clean.re_ca,Sex!=7)                #64,639 records   #!remove (5) females from subs encounter

#remove unknown Banding Flyway
clean.re_ca<-subset(clean.re_ca,B.Flyway!=0)           # No one removed

# Add Info include ==
#0-Fed band, 4-control band, 14-oral swab, 16-trach swab, 18-blood sample, 70-spotlight capture
clean.re_ca<-subset(clean.re_ca, Add.Info==0 | Add.Info==4 | Add.Info==6 |Add.Info==7 |
                      Add.Info==8 | Add.Info==14 |Add.Info==15 | Add.Info==16 | 
                      Add.Info==18 | Add.Info==19 | Add.Info==25 |
                      Add.Info==39 | Add.Info==51 |
                      Add.Info==70| Add.Info==71| Add.Info==89)

# restrict B.month 1 Jun-30 Sept
clean.re_ca$Regime <- ifelse(clean.re_ca$B.Month==6|clean.re_ca$B.Month==7| clean.re_ca$B.Month==8| clean.re_ca$B.Month==9,"B",
                             ifelse(clean.re_ca$B.Month==10| clean.re_ca$B.Month==11| clean.re_ca$B.Month==12,"F",
                                    ifelse(clean.re_ca$B.Month==1| clean.re_ca$B.Month==2| clean.re_ca$B.Month==3,"W",'E'))) 

#remove unknown age
clean.re_ca<-subset(clean.re_ca,Age!=0)           #60,261 records #!remove 51 unknown age

# create new variable, banding age: 0 = local, 1=juv, 2=adult (0 TY and ATY) 
clean.re_ca$b_age   <- ifelse(clean.re_ca$Age_VAGE=="After Hatch Year", 2,
                              ifelse(clean.re_ca$Age_VAGE=="After Second Year", 2,
                                     ifelse(clean.re_ca$Age_VAGE=="After Third Year", 2,
                                            ifelse(clean.re_ca$Age_VAGE=="Third Year", 2,
                                                   ifelse(clean.re_ca$Age_VAGE=="Hatch Year", 1,
                                                          ifelse(clean.re_ca$Age_VAGE=="Juvenile (obsolete)",1,
                                                                 ifelse(clean.re_ca$Age_VAGE=="Local",1,
                                                                        ifelse(clean.re_ca$Age_VAGE=="Second Year" & clean.re_ca$Regime == 'W',1, 2))))))))
xy <- clean.re_ca[,c(33,32)]

clean.re_ca$TYPE <- ifelse(clean.re_ca$`Band Type` == '00' | clean.re_ca$`Band Type` == '11'| clean.re_ca$`Band Type` == '12' | clean.re_ca$`Band Type` == '13'| clean.re_ca$`Band Type` == '18'|clean.re_ca$`Band Type` == '21', 1,
                           ifelse(clean.re_ca$`Band Type` == '01' | clean.re_ca$`Band Type` == '04'| clean.re_ca$`Band Type` == '08' | clean.re_ca$`Band Type` == '51' | clean.re_ca$`Band Type` == '53', 2,
                                  3))

#------------------------------------------------------------------------------
# Data manipulation of recovery data 
#------------------------------------------------------------------------------
raw.recovs<-read.csv("https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data/encounters_minimal.csv")  #reading in CSV
tail(raw.recovs) #Unknown number records
dim(raw.recovs)  #88,411 records 71 variables
summary(raw.recovs)

########################################################################
#cleaning data
########################################################################
#only use status 3 birds
re_cai.recov <-subset(raw.recovs,Status==3)              # 86,428
nrow(re_cai.recov)
# only use HSS = 1
re_cai.recov <-subset(re_cai.recov, Hunt..Season.Surv.>= 1)            
# Exclude unlikely 99's

re_cai.recov <-subset(re_cai.recov, Hunt..Season.Surv.!= 99 | Hunt..Season.Surv.== 99 & (R.Year - B.Year) < 20) 
re_cai.recov <- subset(re_cai.recov, Band.Type.Orig != 98)

#only use B.Year from 1951 onwards
re_cai.recov<-subset(re_cai.recov,B.Year>=1951)          # 82,107 
nrow(re_cai.recov)
#remove unknown sex
re_cai.recov<-subset(re_cai.recov,Sex!=0)                # 81,350
nrow(re_cai.recov)
#remove sex determined subsequent
re_cai.recov<-subset(re_cai.recov,Sex!=6)           # 81,345 records   #!remove (1) males from subs encounter
re_cai.recov<-subset(re_cai.recov,Sex!=7)           # 81,339 records   #!remove (0) females from subs encounter

re_cai.recov$Sex<-factor(re_cai.recov$Sex)

nrow(re_cai.recov)
#remove unknown Banding Flyway 0 and Central America 7

re_cai.recov$TYPE  <-   ifelse(  re_cai.recov$Band.Type.Orig == 11 | 
                                   re_cai.recov$Band.Type.Orig == 12 |
                                   re_cai.recov$Band.Type.Orig == 13 | 
                                   re_cai.recov$Band.Type.Orig == 18 |
                                   re_cai.recov$Band.Type.Orig == 21 , 1,
                                 ifelse( re_cai.recov$Band.Type.Orig ==  1 | 
                                           re_cai.recov$Band.Type.Orig == 4|
                                           re_cai.recov$Band.Type.Orig == 8 |
                                           re_cai.recov$Band.Type.Orig == 51 | 
                                           re_cai.recov$Band.Type.Orig == 53, 2, 3))


#remove Unknown Recov Flyway 0, Central America 7, South America 8, Caribbean 9
re_cai.recov <-subset(re_cai.recov, R.Flyway!=0)
re_cai.recov <-subset(re_cai.recov, R.Flyway!=7)
re_cai.recov <-subset(re_cai.recov,!grepl('X',R.Flyway))
re_cai.recov <-subset(re_cai.recov, R.Flyway!=8)
re_cai.recov <-subset(re_cai.recov, R.Flyway!=9)
re_cai.recov$R.Flyway<-factor(re_cai.recov$R.Flyway)
table (re_cai.recov$R.Flyway)
nrow(re_cai.recov)


re_cai.recaps     <- subset(re_cai.recov,How.Obt==66)  # 29020 # recaps
re_cai.found.dead <- subset(re_cai.recov,How.Obt==0)   # 2552  # found dead

re_cai.other.dead <- subset(re_cai.recov, How.Obt==2|How.Obt==3|How.Obt==4|How.Obt==5|How.Obt==6|How.Obt==7|
                              How.Obt==9|How.Obt==11|How.Obt==12|How.Obt==13|How.Obt==14|How.Obt==15|
                              How.Obt==17|How.Obt==18|How.Obt==20|How.Obt==23|How.Obt==26|How.Obt==30|How.Obt==64) # 862

re_cai.other.dead <- subset(re_cai.other.dead , Pres..Cond. <= 5)  # 800 Misc deaths

re_cai.recov <- subset(re_cai.recov,How.Obt==1)        # 46,522 harvest reports
nrow(re_cai.recov)
# create new variable, banding age: 1=juv, 2=adult (0 TY and ATY) 
re_cai.recov$b_age[(re_cai.recov$Age..VAGE=="After Hatch Year")] <- 2
re_cai.recov$b_age[(re_cai.recov$Age..VAGE=="After Second Year")] <- 2
re_cai.recov$b_age[(re_cai.recov$Age..VAGE=="Hatch Year")] <- 1
re_cai.recov$b_age[(re_cai.recov$Age..VAGE=="Local")] <- 1


#view n shot birds not dead -condition
table(re_cai.recov$Condition..VCondition)

re_cai.recaps <- subset(re_cai.recaps, grepl('Alive-Released', Condition..VCondition))
re_cai.recov  <- subset(re_cai.recov, grepl('Dead', Condition..VCondition))


# Double Check
# Add Info include ==
#0-Fed band, 4-control band, 6- State/Prov Band, 7 - double marked, 8 - temp mark, 14-oral swab, 16-trach swab, 18-blood sample, 70-spotlight capture     
re_cai.recov<-subset(re_cai.recov,Add.Info==0 | Add.Info==4 | Add.Info==6 |Add.Info==7 |
                       Add.Info==8 | Add.Info==14 |Add.Info==15 | Add.Info==16 | 
                       Add.Info==18 | Add.Info==19 | Add.Info==25 |
                       Add.Info==39 | Add.Info==51 |
                       Add.Info==70| Add.Info==71| Add.Info==89)   # 47317

nrow(re_cai.recov)
# Currently not using recaps
re_cai.recaps<-subset(re_cai.recaps,Add.Info==0 | Add.Info==4 | Add.Info==6 |Add.Info==7 |
                        Add.Info==8 | Add.Info==14 | Add.Info==15 |Add.Info==16 | 
                        Add.Info==18 | Add.Info==19 | Add.Info==25 |
                        Add.Info==39 | Add.Info==51 |
                        Add.Info==70| Add.Info==71| Add.Info==89)   # 25875

nrow(re_cai.recaps)

#create new variable 'n_recov' or 'n_recov'
re_cai.recov$n_recov <- 1
re_cai.recaps$n_recap <- 1

# Verify Seasonality of banding-regime (Breeding, Winter, E (Early Season), F (Late Season)): currently not using E and F releases
#restrict B.month 1 Jun-31 Aug
re_cai.recov$Regime <- ifelse(re_cai.recov$B.Month==6| re_cai.recov$B.Month==7| re_cai.recov$B.Month==8| re_cai.recov$B.Month==9,"B",
                              ifelse(re_cai.recov$B.Month==10| re_cai.recov$B.Month==11| re_cai.recov$B.Month==12,"F",
                                     ifelse(re_cai.recov$B.Month==1| re_cai.recov$B.Month==2| re_cai.recov$B.Month==3,"W",'E'))) 


re_cai.recaps$Regime <- ifelse(re_cai.recaps$B.Month==6| re_cai.recaps$B.Month==7| re_cai.recaps$B.Month==8| re_cai.recaps$B.Month==9,"B",
                               ifelse(re_cai.recaps$B.Month==10| re_cai.recaps$B.Month==11| re_cai.recaps$B.Month==12,"F",
                                      ifelse(re_cai.recaps$B.Month==1| re_cai.recaps$B.Month==2| re_cai.recaps$B.Month==3,"W",'E'))) 


re_cai.recaps$Rec.Regime <- ifelse(re_cai.recaps$R.Month==6| re_cai.recaps$R.Month==7| re_cai.recaps$R.Month==8| re_cai.recaps$R.Month==9,"B",
                                   ifelse(re_cai.recaps$R.Month==10| re_cai.recaps$R.Month==11| re_cai.recaps$R.Month==12,"F",
                                          ifelse(re_cai.recaps$R.Month==1| re_cai.recaps$R.Month==2| re_cai.recaps$R.Month==3,"W",'E'))) 



re_cai.recov$b_age   <- ifelse(re_cai.recov$Age..VAGE=="After Hatch Year", 2,
                               ifelse(re_cai.recov$Age..VAGE=="After Second Year", 2,
                                      ifelse(re_cai.recov$Age..VAGE=="After Third Year", 2,
                                             ifelse(re_cai.recov$Age..VAGE=="Third Year", 2,
                                                    ifelse(re_cai.recov$Age..VAGE=="Hatch Year", 1,
                                                           ifelse(re_cai.recov$Age..VAGE=="Juvenile (obsolete)",1,
                                                                  ifelse(re_cai.recov$Age..VAGE=="Local",1,
                                                                         ifelse(re_cai.recov$Age..VAGE=="Second Year" & re_cai.recov$Regime == 'W',1, 2))))))))


re_cai.recaps$b_age   <- ifelse(re_cai.recaps$Age..VAGE=="After Hatch Year", 2,
                                ifelse(re_cai.recaps$Age..VAGE=="After Second Year", 2,
                                       ifelse(re_cai.recaps$Age..VAGE=="After Third Year", 2,
                                              ifelse(re_cai.recaps$Age..VAGE=="Third Year", 2,
                                                     ifelse(re_cai.recaps$Age..VAGE=="Hatch Year", 1,
                                                            ifelse(re_cai.recaps$Age..VAGE=="Juvenile (obsolete)",1,
                                                                   ifelse(re_cai.recaps$Age..VAGE=="Local",1,
                                                                          ifelse(re_cai.recaps$Age..VAGE=="Second Year" & re_cai.recaps$Regime == 'W',1, 2))))))))
table(re_cai.recov$b_age)




#restrict R.month 1 Sep -31 Jan,  Fall 93, Winter 92, Hunting Season 94   #14,172
re_cai.recov <-subset(re_cai.recov, R.Month==9|R.Month==10| R.Month==11|R.Month==12|R.Month==1| R.Month==92|R.Month== 93|R.Month==94) 
nrow(re_cai.recov)

#remove unknown age
re_cai.recov<-subset(re_cai.recov,Age!=0)        

re_cai.recaps <-subset(re_cai.recaps,Age!=0)  #1347
nrow(re_cai.recov)

clean.re_ca_sub <- clean.re_ca
re_cai.recov_sub <- re_cai.recov
re_cai.recaps_sub <- re_cai.recaps


clean.re_ca_sub$PPH <-    ifelse(clean.re_ca_sub$Region..State == 'Alberta'| 
                                   clean.re_ca_sub$Region..State == 'Saskatchewan'|
                                   clean.re_ca_sub$Region..State == 'Manitoba'|
                                   clean.re_ca_sub$Region..State == 'Montana'|
                                   clean.re_ca_sub$Region..State == 'North Dakota'|
                                   clean.re_ca_sub$Region..State == 'South Dakota'|
                                   clean.re_ca_sub$Region..State == 'Minnesota',1,0)

re_cai.recov_sub$PPH <-      ifelse(re_cai.recov_sub$BRegion..State == 'Alberta'| 
                                      re_cai.recov_sub$BRegion..State == 'Saskatchewan'|
                                      re_cai.recov_sub$BRegion..State == 'Manitoba'|
                                      re_cai.recov_sub$BRegion..State == 'Montana'|
                                      re_cai.recov_sub$BRegion..State == 'North Dakota'|
                                      re_cai.recov_sub$BRegion..State == 'South Dakota'|
                                      re_cai.recov_sub$BRegion..State == 'Minnesota',1,0)


re_cai.recaps$PPH  <-         ifelse(re_cai.recaps$BRegion..State == 'Alberta'| 
                                       re_cai.recaps$BRegion..State == 'Saskatchewan'|
                                       re_cai.recaps$BRegion..State == 'Manitoba'|
                                       re_cai.recaps$BRegion..State == 'Montana'|
                                       re_cai.recaps$BRegion..State == 'North Dakota'|
                                       re_cai.recaps$BRegion..State == 'South Dakota'|
                                       re_cai.recaps$BRegion..State == 'Minnesota',1,0)

clean.re_ca_sub$sex    <- recode(clean.re_ca_sub$Sex..VSEX, "Female" = 1, 'Male' = 2)
clean.re_ca_sub$species<- recode(clean.re_ca_sub$`Species Game Birds::Species` , "Canvasback" = 1, 'Redhead' = 2)
clean.re_ca_sub$age    <- clean.re_ca_sub$b_age
clean.re_ca_sub$rel.type <- ifelse(clean.re_ca_sub$Regime == 'B' & clean.re_ca_sub$b_age == 1 & clean.re_ca_sub$B.Month <8, 'Early', 
                                   ifelse(clean.re_ca_sub$Regime == 'B' & clean.re_ca_sub$b_age == 1 & clean.re_ca_sub$B.Month == 8 & clean.re_ca_sub$`B Day Code` < 4, 'Early',
                                          ifelse(clean.re_ca_sub$Regime == 'B' & clean.re_ca_sub$b_age == 1 & clean.re_ca_sub$B.Month == 8 & clean.re_ca_sub$`B Day Code` > 3, 'Late',
                                                 ifelse(clean.re_ca_sub$Regime == 'B' & clean.re_ca_sub$b_age == 1 & clean.re_ca_sub$B.Month > 8, 'Late',
                                                        ifelse(clean.re_ca_sub$Regime == 'B' & clean.re_ca_sub$b_age == 2 & clean.re_ca_sub$B.Month < 8, 'Early',
                                                               ifelse(clean.re_ca_sub$Regime == 'B' & clean.re_ca_sub$b_age == 2 & clean.re_ca_sub$B.Month == 8 & clean.re_ca_sub$`B Day Code` < 4, 'Early',
                                                                      ifelse(clean.re_ca_sub$Regime == 'B' & clean.re_ca_sub$b_age == 2 & clean.re_ca_sub$B.Month == 8 & clean.re_ca_sub$`B Day Code` > 3, 'Late',
                                                                             ifelse(clean.re_ca_sub$Regime == 'B' & clean.re_ca_sub$b_age == 2 & clean.re_ca_sub$B.Month > 8, 'Late', 
                                                                                    clean.re_ca_sub$B.Month  ))))))))


re_cai.recov_sub$sex    <- recode(re_cai.recov_sub$Sex..VSEX, "Female" = 1, 'Male' = 2)
re_cai.recov_sub$species<- recode(re_cai.recov_sub$Species.Game.Birds..Species , "Canvasback" = 1, 'Redhead' = 2)
re_cai.recov_sub$age    <- re_cai.recov_sub$b_age 

re_cai.recov_sub$rel.type <- ifelse(re_cai.recov_sub$Regime == 'B' & re_cai.recov_sub$b_age == 1 & re_cai.recov_sub$B.Month <8, 'Early', 
                                    ifelse(re_cai.recov_sub$Regime == 'B' & re_cai.recov_sub$b_age == 1 & re_cai.recov_sub$B.Month == 8 & re_cai.recov_sub$B.Day < 16, 'Early',
                                           ifelse(re_cai.recov_sub$Regime == 'B' & re_cai.recov_sub$b_age == 1 & re_cai.recov_sub$B.Month == 8 & re_cai.recov_sub$B.Day > 15, 'Late',
                                                  ifelse(re_cai.recov_sub$Regime == 'B' & re_cai.recov_sub$b_age == 1 & re_cai.recov_sub$B.Month >8, 'Late',
                                                         ifelse(re_cai.recov_sub$Regime == 'B' & re_cai.recov_sub$b_age == 2 & re_cai.recov_sub$B.Month < 8, 'Early',
                                                                ifelse(re_cai.recov_sub$Regime == 'B' & re_cai.recov_sub$b_age == 2 & re_cai.recov_sub$B.Month == 8 & re_cai.recov_sub$B.Day < 16, 'Early',
                                                                       ifelse(re_cai.recov_sub$Regime == 'B' & re_cai.recov_sub$b_age == 2 & re_cai.recov_sub$B.Month == 8 & re_cai.recov_sub$B.Day > 15, 'Late',
                                                                              ifelse(re_cai.recov_sub$Regime == 'B' & re_cai.recov_sub$b_age == 2 & re_cai.recov_sub$B.Month > 8, 'Late', 
                                                                                     re_cai.recov_sub$B.Month  ))))))))



re_cai.recaps$sex    <- recode(re_cai.recaps$Sex..VSEX, "Female" = 1, 'Male' = 2)
re_cai.recaps$species<- recode(re_cai.recaps$Species.Game.Birds..Species , "Canvasback" = 1, 'Redhead' = 2)
re_cai.recaps$age    <- re_cai.recaps$b_age 

re_cai.recaps$rel.type <- ifelse(re_cai.recaps$Regime == 'B' & re_cai.recaps$b_age == 1 & re_cai.recaps$B.Month <8, 'Early', 
                                 ifelse(re_cai.recaps$Regime == 'B' & re_cai.recaps$b_age == 1 & re_cai.recaps$B.Month == 8 & re_cai.recaps$B.Day < 16, 'Early',
                                        ifelse(re_cai.recaps$Regime == 'B' & re_cai.recaps$b_age == 1 & re_cai.recaps$B.Month == 8 & re_cai.recaps$B.Day > 15, 'Late',
                                               ifelse(re_cai.recaps$Regime == 'B' & re_cai.recaps$b_age == 1 & re_cai.recaps$B.Month >8, 'Late',
                                                      ifelse(re_cai.recaps$Regime == 'B' & re_cai.recaps$b_age == 2 & re_cai.recaps$B.Month < 8, 'Early',
                                                             ifelse(re_cai.recaps$Regime == 'B' & re_cai.recaps$b_age == 2 & re_cai.recaps$B.Month == 8 & re_cai.recaps$B.Day < 16, 'Early',
                                                                    ifelse(re_cai.recaps$Regime == 'B' & re_cai.recaps$b_age == 2 & re_cai.recaps$B.Month == 8 & re_cai.recaps$B.Day > 15, 'Late',
                                                                           ifelse(re_cai.recaps$Regime == 'B' & re_cai.recaps$b_age == 2 & re_cai.recaps$B.Month > 8, 'Late', 
                                                                                  re_cai.recaps$B.Month  ))))))))



start <- 1961

clean_releases_b <- clean.re_ca_sub %>%
  subset(Regime == "B" & PPH == 1) %>%
  subset(B.Year >= start) %>%
  group_by(sex,species,age) %>% #TYPE
  group_split()

#clean.re_ca_sub$Flyway <-recode(clean.re_ca_sub$B.Flyway, `1` = 1, `2` = 2, `3`= 3, `4` = 3,`5` =5,  `6` = 6)
clean_releases_w <-  clean.re_ca_sub %>%
  subset(Regime == "W") %>%
  subset(B.Flyway < 4) %>%
  subset(Region..State == 'Maryland'| Region..State == 'New York'|  Region..State == 'Ohio'|  Region..State == 'Illinois'| Region..State == 'Virginia'| Region..State == 'Texas')%>%
  #subset(Region..State == 'Maryland'| Region..State == 'New York')%>%
  subset(B.Year >= start+1) %>%
  group_by(sex,species) %>% #,age) %>% #TYPE
  group_split()


start <- 1961
end <- 2020
n.years <- end - start

re_cai.recov_sub <- subset(re_cai.recov_sub, R.Month < 4 | R.Month > 8)
re_cai.recov_sub$R.Year[which(re_cai.recov_sub$R.Month < 4)] <- re_cai.recov_sub$R.Year[which(re_cai.recov_sub$R.Month < 4)] - 1
re_cai.recov_sub <- subset(re_cai.recov_sub, R.Year < end)
re_cai.recov_sub <- subset(re_cai.recov_sub,R.Year != B.Year - 1)


clean_recovs_b <- re_cai.recov_sub %>%
  subset(Regime == "B" & PPH == 1) %>%
  subset(B.Year >= start) %>%
  group_by(sex,species,age) %>%  #TYPE
  group_split()

clean_recovs_w <- re_cai.recov_sub %>%
  subset(Regime == "W") %>%
  subset(B.Flyway < 4) %>%
  subset(BRegion..State == 'Maryland'| BRegion..State == 'New York'|  BRegion..State == 'Ohio'|  BRegion..State == 'Illinois'|  BRegion..State == 'Virginia'| BRegion..State == 'Texas')%>%
  #subset(BRegion..State == 'Maryland'| BRegion..State == 'New York')%>%
  subset(B.Year >= start+1) %>%
  group_by(sex,species) %>% #,age) %>%  #TYPE
  group_split()
start   <- 1961
end     <- 2020
n.years <- end - start

# re_cai.recaps <- subset(re_cai.recaps, R.Month < 4 | R.Month > 8)
# re_cai.recaps$R.Year[which(re_cai.recaps$R.Month < 4)] <- re_cai.recaps$R.Year[which(re_cai.recaps$R.Month < 4)] - 1
# re_cai.recaps <- subset(re_cai.recaps, R.Year < end)
# re_cai.recaps <- subset(re_cai.recaps, R.Year != B.Year)
# 
# 
# clean_recaps_w <- re_cai.recaps %>%
#   subset(Regime == "W") %>%
#   subset(B.Flyway < 4) %>%
#   subset(B.Year > start) %>%
#   group_by(Band) %>%
#   slice(which.min(R.Year)) %>%
#   group_by(sex,species,age) %>%  
#   group_split()



ordering.b <- expand.grid( age = c('HY','AHY'), species = c('CANV','REDH'), sex = c('female','male'))
ordering.w <- expand.grid( species = c('CANV','REDH'), sex = c('female','male')) # age = c("SY",'ASY'), 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Breeding Season 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
marr.main <- array(0, dim = c(n.years, n.years+1, nrow(ordering.b)))
rel.main <- array(0, dim = c(n.years,nrow(ordering.b)))
dim(marr.main)


for(j in 1:nrow(ordering.b)){
  for (i in 1:nrow(clean_releases_b[[j]] )){
    rel.main[clean_releases_b[[j]]$B.Year[i] - start + 1,j] <- rel.main[clean_releases_b[[j]]$B.Year[i]- start + 1,j] + clean_releases_b[[j]]$Count.of.Birds[i]
  }
  
  for (i in 1:nrow(clean_recovs_b[[j]])){
    marr.main[(clean_recovs_b[[j]]$B.Year[i] - (start - 1)), (clean_recovs_b[[j]]$R.Year[i] - (start - 1)),j] <- marr.main[(clean_recovs_b [[j]]$B.Year[i] - (start - 1)), (clean_recovs_b[[j]]$R.Year[i] - (start - 1)),j] + 1
  }
  marr.main[,dim(marr.main)[2],j] <- rel.main[,j] - rowSums(marr.main[,,j])
  
  for (t in 2:n.years){
    marr.main[t,1:(t-1),j] <- 0
  }
}  
colSums(rel.main)


ordering.b <- cbind.data.frame(ordering.b, colSums(rel.main))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Winter Season 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
start <- 1962
end   <- 2021
n.years <- end - start

marr.off <- array(0, dim = c(n.years,  n.years+1, nrow(ordering.w)))
rel.off <- array(0, dim = c(n.years,nrow(ordering.w)))
dim(marr.off)


for(j in 1:nrow(ordering.w)){
  for (i in 1:nrow(clean_releases_w[[j]] )){
    rel.off[clean_releases_w[[j]]$B.Year[i] - start + 1,j] <- rel.off[clean_releases_w[[j]]$B.Year[i]- start + 1,j] + clean_releases_w[[j]]$Count.of.Birds[i]
  }
  
  for (i in 1:nrow(clean_recovs_w[[j]])){
    marr.off[(clean_recovs_w[[j]]$B.Year[i] - (start -1)), (clean_recovs_w[[j]]$R.Year[i] - (start -1)),j] <-  marr.off[(clean_recovs_w [[j]]$B.Year[i] - (start -1)), (clean_recovs_w[[j]]$R.Year[i] - (start -1)),j] + 1
    #marr.off[(clean_recaps_w[[j]]$B.Year[i] - (start -1)), n.years + (clean_recaps_w[[j]]$R.Year[i] - (start-1)),j] <- marr.off[(clean_recaps_w [[j]]$B.Year[i] - (start-1)), n.years + (clean_recaps_w[[j]]$R.Year[i] - (start-1)),j] + 1
  }
  marr.off[,dim(marr.off)[2],j] <- rel.off[,j] - rowSums(marr.off[,,j])
  
  for (t in 2:n.years){
    marr.off[t,1:(t-1),j] <- 0
    
  }
}  
colSums(rel.off)

ordering.w <- cbind.data.frame(ordering.w, colSums(rel.off))
get.first <- function(x) min(which(x!=0))
get.last <- function(x) max(which(x!=0))

f <- apply(rel.main, 2, get.first)
l <- apply(rel.main, 2, get.last)

f.off <- apply(rel.off, 2, get.first)
l.off <- apply(rel.off, 2, get.last)
###########################################################################################################################

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}
GeoDistanceInMetresMatrix <- function(df.geopoints){
  GeoDistanceInMetres <- function(g1, g2){
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  df.geopoints$index <- 1:n.geopoints
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  return(mat.distances)
}

start <- 1961

pond_data <- read.csv('https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data/Ponds-stratum_bpop.csv')
duck_data <- read.csv('https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data/CANV_REDH-stratum_bpop.csv')

pond_data <- subset(pond_data, YEAR >= start)
duck_data <- subset(duck_data, YEAR >= start)

pond_data <- subset(pond_data, STRATUM > 18)
pond_data <- subset(pond_data, !grepl('50', STRATUM))

duck_data <- subset(duck_data, STRATUM > 18)
duck_data <- subset(duck_data, !grepl('50', STRATUM))

#create a couple temp files
temp <- tempfile()
temp2 <- tempfile()
download.file("https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data//stratum.zip",temp)
unzip(zipfile = temp, exdir = temp2)

bpop <-list.files(temp2, pattern = ".shp$",full.names=TRUE)
bpop <- read_sf(bpop)

code <- c(11,11,12,10,11,11,7,6,4,5,13,14,15,9,2,3,2,16,8,17,18,19,78,22,24,68,21,77,70,67,25,71,60,76,23,59,27,72,31,51,26,71,32,28,69,58,73,37,33,29,65,38,35,34,64,30,66,52,41,57,36,63,42,39,40,1,46,43,48,44,1,1,47,53,54,56,1,45,55,49,50)
code <- code - 1
bpop$code <- code

keep <- unique(pond_data$STRATUM)

bpop_sub <- bpop[bpop$code %in% keep,]
plot(bpop_sub )

keep2 <- unique(duck_data$STRATUM)
bpop_sub2 <- bpop[bpop$code %in% keep2,]
plot(bpop_sub2 )
W     = st_touches(bpop_sub$geometry , sparse=FALSE)
listW = mat2listw(W)

bpop_coords = bpop_sub %>% st_centroid() %>% st_coordinates()

bpop_coords.1 <- cbind.data.frame(strata = 1:26, lat = bpop_coords[,2], lon = bpop_coords[,1])

dist.pond <- GeoDistanceInMetresMatrix(bpop_coords.1)/100000

dist.prime <- dist.pond[-c(1:2),-c(1:2)] 
strata <- cbind.data.frame(STRATUM = bpop_sub$code, STRATA =  1:length(bpop_sub$code))

pond_data <- left_join(pond_data, strata)

ponds <- cast(pond_data,  YEAR ~ STRATA, value = 'POP')
ponds_se <- cast(pond_data,  YEAR ~ STRATA, value = 'SEPOP')

lat.p <- scale(bpop_coords.1$lat)
ponds <- ponds[,-c(1:3)]
ponds_se <- ponds_se[,-c(1:3)]

all_obs <- rowSums(ponds)
all_error <- sqrt(rowSums(ponds_se^2))
ponds $YEAR <- NULL
ponds  <- (data.matrix(ponds/1000))

all_obs <- rbind.data.frame(cbind.data.frame(Year = 1961:2015, Ponds = all_obs/1000, SE = all_error/1000),
                            cbind.data.frame(Year = 2016:2019, Ponds = c(5012.50,6096.00,5227.40,4990.30), SE = c(156.40,182.70,173.00,172.10)))

# Read in dataset
code <- c(11,11,12,10,11,11,7,6,4,5,13,14,15,9,2,3,2,16,8,17,18,19,78,22,24,68,21,77,70,67,25,71,60,76,23,59,27,72,31,51,26,71,32,28,69,58,73,37,33,29,65,38,35,34,64,30,66,52,41,57,36,63,42,39,40,1,46,43,48,44,1,1,47,53,54,56,1,45,55,49,50)
code <- code - 1
bpop$code <- code
bpop$area <- st_area(bpop$geometry)/(1000 * 1000)

keep <- unique(pond_data$STRATUM)

bpop_sub  <- bpop[bpop$code %in% keep,]
keep2     <- unique(duck_data$STRATUM)
bpop_sub2 <- bpop[bpop$code %in% keep2,]
area <- scale(bpop_sub$area [-c(1:2)])
W     = st_touches(bpop_sub2$geometry , sparse=FALSE)
listW = mat2listw(W)

adj  = W * 1L
nadj = (rowSums(W))

bpop_coords = bpop_sub2 %>% st_centroid() %>% st_coordinates()
plot(st_geometry(bpop_sub2), col = code)
plot(listW, bpop_coords, add=TRUE, col='blue', pch=16)

bpop_coords.1 <- cbind.data.frame(strata = 1:33, lat = bpop_coords[,2], lon = bpop_coords[,1])

dist <- GeoDistanceInMetresMatrix(bpop_coords.1)/100000

strata.d <- cbind.data.frame(STRATUM = bpop_sub2$code, STRATA =  1:length(bpop_sub2$code), Area = bpop_sub2$area)

duck_data <- left_join(duck_data, strata.d)
duck_data$DEN <- duck_data$POP/duck_data$Area

STRATA_KEEP <- paste(strata$STRATUM[-c(1:2)], collapse = '|')

ducks <- cast(duck_data, SPECIES ~ YEAR ~ STRATA, value = 'POP')
duck.p <- cast(duck_data,  YEAR+SPECIES ~ STRATA, value = 'POP')
ducks.r <- round(ducks,0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
US_Data <- read.csv('https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data/US_CAN_RED_HARV_RAW.csv')

wide <- cast(US_Data, Season ~ Flyway ~ Species ~ Sex ~ Age, value = 'Number_parts')

CANV_FA <-wide[,-2,1,1,1]
CANV_FJ <-wide[,-2,1,1,2]
CANV_FU <-wide[,-2,1,1,3]

CANV_MA <-wide[,-2,1,2,1]
CANV_MJ <-wide[,-2,1,2,2]
CANV_MU <-wide[,-2,1,2,3]

CANV_UA <-wide[,-2,1,3,1]
CANV_UJ <-wide[,-2,1,3,2]
CANV_UU <-wide[,-2,1,3,3]

REDH_FA <-wide[,-2,2,1,1]
REDH_FJ <-wide[,-2,2,1,2]
REDH_FU <-wide[,-2,2,1,3]

REDH_MA <-wide[,-2,2,2,1]
REDH_MJ <-wide[,-2,2,2,2]
REDH_MU <-wide[,-2,2,2,3]

REDH_UA <-wide[,-2,2,3,1]
REDH_UJ <-wide[,-2,2,3,2]
REDH_UU <-wide[,-2,2,3,3]

FJ <- abind(CANV_FJ,REDH_FJ, along = 3)
FA <- abind(CANV_FA,REDH_FA, along = 3)
MJ <- abind(CANV_MJ,REDH_MJ, along = 3)
MA <- abind(CANV_MA,REDH_MA, along = 3)

FJ[is.na(FJ)] <- 0
FA[is.na(FA)] <- 0
MJ[is.na(MJ)] <- 0
MA[is.na(MA)] <- 0


t_rr <- read.csv('https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data/trophy_duck_reporting_rates.csv')
t_rr_sub <- subset(t_rr, ï..Year >= start)


h.wide <- cast(US_Data, Season ~ Species ~ Sex ~ Age, value = 'Harvest_estimate', fun.aggregate = 'sum')

H.CANV_FA <-h.wide[,1,1,1]
H.CANV_FJ <-h.wide[,1,1,2]
H.CANV_FU <-h.wide[,1,1,3]

H.CANV_MA <-h.wide[,1,2,1]
H.CANV_MJ <-h.wide[,1,2,2]
H.CANV_MU <-h.wide[,1,2,3]

H.CANV_UA <-h.wide[,1,3,1]
H.CANV_UJ <-h.wide[,1,3,2]
H.CANV_UU <-h.wide[,1,3,3]

H.REDH_FA <-h.wide[,2,1,1]
H.REDH_FJ <-h.wide[,2,1,2]
H.REDH_FU <-h.wide[,2,1,3]

H.REDH_MA <-h.wide[,2,2,1]
H.REDH_MJ <-h.wide[,2,2,2]
H.REDH_MU <-h.wide[,2,2,3]

H.REDH_UA <-h.wide[,2,3,1]
H.REDH_UJ <-h.wide[,2,3,2]
H.REDH_UU <-h.wide[,2,3,3]

library(abind)
H.FJ <- cbind(H.CANV_FJ,H.REDH_FJ)
H.FA <- cbind(H.CANV_FA,H.REDH_FA)
H.MJ <- cbind(H.CANV_MJ,H.REDH_MJ)
H.MA <- cbind(H.CANV_MA,H.REDH_MA)

HARVEST.J <- abind(H.FJ, H.MJ, along = 3)
HARVEST.A <- abind(H.FA, H.MA, along = 3)

# time - species - sex - age
HARVEST <- abind(HARVEST.J, HARVEST.A, along = 4)

FJ[is.na(FJ)] <- 0

cut <- 15

marr.main.abr <- marr.main

for(i in 1:(dim(marr.main.abr )[1]-cut)){
  for(j in (i + cut):dim(marr.main.abr )[1]){
    for(k in 1:8){
      marr.main.abr[i,j,k] <- 0
    }
  }
}
marr.main.abr[4,,2] <- 0

marr.off.edit <- abind(array(0, dim = c(dim( marr.off)[1],1,4)), marr.off, along = 2)
marr.off.edit <- marr.off.edit[,-(dim(marr.off.edit)[2] - 1),]

for(i in 1:(dim(marr.off.edit)[1]-cut)){
  for(j in (i + cut):dim(marr.off.edit )[1]){
    for(k in 1:4){
      marr.off.edit[i,j,k] <- 0
    }
  }
}

#########################################################

# Average (Nov-Apr)  Pacific North American Index 1961 - 2020
pna <- c(0.393333333,-0.605,-0.078333333,0.121666667,-0.778333333,-0.431666667,-1,-0.243333333,-0.348333333,0.658333333,-0.461666667,-1.028333333,0.191666667,
         -0.685,-0.051666667,-0.348333333,0.938333333,0.221666667,-0.541666667,0.538333333,0.743333333,-0.615,1.065,1.063333333,-0.341666667,0.318333333,0.881666667,
         0.988333333,-0.476666667,-0.35,-0.27,0.5,0.388333333,-0.07,0.171666667,0.025,-0.206666667,0.971666667,0.28,0.32,0.591666667,-0.295,0.87,0.236666667,0.481666667,
         0.186666667,0.421666667,-0.046666667,-0.241666667,0.99,-0.635,0.196666667,-0.438333333,-0.33,0.361666667,0.948333333,0.383333333,-0.6,0.081666667,-0.55)

load(url("https://raw.githubusercontent.com/gibson-d/CANV-REDH/main//Data/CHB_SAL.rdata"))
df1 <- CHB_SAL
df1$SampleDepth <- round_any( df1$SampleDepth, 0.5)
df1$lat_round <- round_any(df1$Latitude,.08333)
df1$long_round <- round_any(df1$Longitude,.08333)
df1$Observation <- paste(df1$lat_round,"-",df1$long_round,'-',df1$SampleDate,'-', df1$SampleDepth)
df1$Unique1 <-  paste(df1$lat_round,"-",df1$long_round)
df1 <- df1%>% group_by(Observation) %>% sample_n(size = 1)

df1 <- subset(df1, Year >= 1961 & Year < 2020 & Longitude < -76.15 &  Longitude > -76.5)
df1 <- df1[complete.cases(df1$ReportedValue),]
df1 <- subset(df1, SampleDepth <= 1 )
df1 <- subset(df1,ReportedValue > .1 & ReportedValue<30 )

df_unique1 <- df1[!duplicated(df1$Unique1),]
study_points1 <- cbind.data.frame(long = df_unique1$long_round, lat = df_unique1$lat_round)
blockcoord1 <- cbind(study_points1$long,study_points1$lat)
neigh1      <- dnearneigh(blockcoord1, d1 = 0, d2 = 20, longlat = TRUE)
winnb1 <- nb2WB(neigh1)          

distance <- as.matrix(dist(study_points1))

covs1 <- df1[,c('ACE','AMO','NAO','PREC','SampleDepth')]
covs1$SampleDepth <- round_any(covs1$SampleDepth,0.5)
covs_scale1 <- scale(covs1)
covs_scale1[is.na(covs_scale1)] <- 0

covs_scale1[,5] <- as.factor(covs_scale1[,5])

year1       <- df1$Year - 1960
month1      <- df1$Month
system.ind  <- ifelse(year1 < 24,2,3)

df <- read.csv('https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data/texas_salinity.csv')
df <- subset(df, Depth == 0.3)
df <- subset(df, !grepl('RIVER', SiteName))
df <- subset(df, !grepl('CREEK', SiteName))

df$lat_round <- round_any(df$Latitude,.08333)
df$long_round <- round_any(df$Longitude,.08333)
df$Unique1 <-  paste(df$lat_round,"-",df$long_round)
df$Observation <- paste(df$lat_round,"-",df$long_round,'-',df$SiteName,'-',df$Year,'-',df$Month)
df <- df%>% group_by(Observation) %>% sample_n(size = 1)

df_unique <- df[!duplicated(df$Unique1),]
study_points <- cbind.data.frame(long = df_unique$long_round, lat = df_unique$lat_round)
blockcoord <- cbind(study_points$long,study_points$lat)
neigh      <- dnearneigh(blockcoord, d1 = 0, d2 = 20, longlat = TRUE)
winnb <- nb2WB(neigh)          

covs <- df[,c('ACE','AMO6','NAO6','PREC6','FLOW6')]
covs_scale <- scale(covs)
covs_scale[is.na(covs_scale )] <- 0

year  <- df$Year - 1960
month <- df$Month

preds <- read.csv("https://raw.githubusercontent.com/gibson-d/CANV-REDH/main/Data/fall_predictions.csv")
preds <- subset(preds, ï..Year > 1960)
preds_lgm <- cbind.data.frame(Ace = preds$Ace, AMO = preds$AMO, NAO = preds$NAO, PREC = preds$TEX.PREC, FLOW = preds$TEX.FLOW)
preds_chb <- cbind.data.frame(Ace = preds$Ace, AMO = preds$AMO, NAO = preds$NAO, PREC = preds$CHB.PREC, FLOW = rep(NA))

preds_lgm_sc <-  preds_chb_sc <- matrix(NA, 59,5)

for(i in 1:dim(preds_lgm)[1]){
  for(j in 1:dim(preds_lgm)[2]){
    preds_lgm_sc[i,j] <- (preds_lgm[i,j] - attr(covs_scale,"scaled:center")[j])/ attr(covs_scale,"scaled:scale")[j]
  }
  for(j in 1:(dim(preds_chb)[2] - 1)){
    preds_chb_sc[i,j] <- (preds_chb[i,j] - attr(covs_scale1,"scaled:center")[j])/ attr(covs_scale1,"scaled:scale")[j]
  }
}

preds_chb_sc[is.na(preds_chb_sc)] <- 0
preds_lgm_sc[is.na(preds_lgm_sc)] <- 0

pred_covs <- abind(preds_lgm_sc,preds_chb_sc, along = 3)


##########################################################################################################
# Clean up
##########################################################################################################

#Duck census data from 2015-2019 (traditional - AK)
end <- matrix(c(658000,	647000,	603000,	574000, 1289000,	1115000,	1000000,	732000), 2,4, byrow = TRUE)

y <- apply(ducks,c(1:2), sum)
y <- abind(y,end,along = 2)
lny <- log(y)
lny[is.infinite(lny)] <- NA
ponds[1,20] <- 0.5
y_scale <- scale(t(y))
all.ponds <- c(rowSums(ponds), 5012.500,6096.000,5227.400,4990.30)
p_scale <- scale(c(rowSums(ponds), 5012.500,6096.000,5227.400,4990.30))
pond.std <- p_scale[,1]

library(splines)

CANS <- y[1,]
REDS <- y[2,]

bs_bbase <- function(x, xl = min(x), xr = max(x), nseg = 5, deg = 3) {
  # Compute the length of the partitions
  dx <- (xr - xl) / nseg
  # Create equally spaced knots
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # Use bs() function to generate the B-spline basis
  get_bs_matrix <- matrix(bs(x, knots = knots, degree = deg, Boundary.knots = c(knots[1], knots[length(knots)])), nrow = length(x))
  # Remove columns that contain zero only
  bs_matrix <- get_bs_matrix[, -c(1:deg, ncol(get_bs_matrix):(ncol(get_bs_matrix) - deg))]
  
  return(bs_matrix)
}

B2 <- bs_bbase(CANS, nseg = 20)
B1 <- bs_bbase(REDS, nseg = 20)

D  <- diff(diag(ncol(B2)), diff = 1)
Q  <- t(D) %*% solve(D %*% t(D))
Z1 <- B1 %*% Q
Z2 <- B2 %*% Q

BM1 <- abind(Z1,Z2,along = 3)

projN <- cbind.data.frame(CANS = scale(CANS), REDS = scale(REDS))
rho.sex <- array(c(1,.5,.5,1), dim = c(2,2,2))
N.std    <- projN

mean_sal <- c(9.942515,61.816014)
std_sal  <- c(2.234635, 17.407757)

# Bundle data 
dat <- list(
  N.std = N.std, 
  pond.std = pond.std, 
  
  lower.pond  =  -1.81511, upper.pond = 2.001904,
  RR = c(t_rr_sub$Mean, rep(t_rr_sub$Mean[50], 9)),
  SD = c(t_rr_sub$SD, rep(t_rr_sub$SD[50], 9)),
  
  marr.wint = marr.off.edit,                       #
  rel.wint = apply( marr.off.edit, c(1,3), sum),   #
  marr.fall = marr.main.abr,                      #
  rel.fall = apply( marr.main.abr, c(1,3), sum),  #
  IF = FJ , F = FJ + FA,
  IM = MJ,  M = MJ + MA,
  FA = FA,  T = FA + MA,
  BM = BM1,
  start.n.std = y_scale[1,],
  y = (lny),#
  pond = log(rbind(ponds,matrix(NA, 4,24))),                               # annual number of ponds in each stratum
  dist = dist.prime,                          # distance matrix b/w stratums
  pna = scale(pna)[,1],
  lat = lat.p[-c(1:2)],                        # latitude of each stratum
  area = area[,1],
  all.ponds = all.ponds,
  
  salinity.lgm = log(df$Salinity), pred.covs = pred_covs,
  salinity.chb = log(df1$ReportedValue),
  covs.lgm = covs_scale, 
  covs.chb = covs_scale1[,c(1:4)]
)

y_inits <- cbind(lny)
y_inits <- ifelse(is.na(y_inits) == TRUE, mean(lny, na.rm = TRUE), NA)

constants <- list(  sex    = as.numeric(ordering.b$sex), cut = cut,
                    spec   = as.numeric(ordering.b$species),
                    age    = as.numeric(ordering.b$age),
                    season = as.numeric(ordering.b$type),
                    flyway = 4, 
                    ages = c(1,2,1,2),
                    sex1 = as.numeric(ordering.w$sex),
                    spec1 = as.numeric(ordering.w$species),
                    age1 = as.numeric(ordering.w$age),
                    season1 = as.numeric(ordering.w$type),
                    n.sex = 2, n.species = 2, n.age = 2, n.season= 2,n.year = 61,
                    K = ncol(BM1),
                    
                    N.lgm = length(table(df$Unique1)), nobs.lgm = length(df$Salinity), lgm_breaks = as.numeric(as.factor(df$Unique1)),
                    L.lgm = length( winnb$weights),  year.lgm = year, month.lgm = month,
                    weights.lgm = winnb$weights, num.lgm = winnb$num, adj.lgm = winnb$adj,
                    
                    N.chb = length(table(df1$Unique1)), nobs.chb = length(df1$ReportedValue), chb_breaks = as.numeric(as.factor(df1$Unique1)),
                    L.chb = length( winnb1$weights),  year.chb = year1, month.chb = month1,
                    weights.chb = winnb1$weights, num.chb = winnb1$num, adj.chb = winnb1$adj,
                    sys.est = c(rep(2,22),rep(3,37)),system.ind=system.ind, depths = covs_scale1[,5],
                    
                    diff = 5, p.scale = sd(rowSums(ponds)),
                    n.occasions = dim(marr.off.edit)[2]-1, n.stratum = 24, n.duck.strat = 33)

pond.inits <- apply(ponds, 2, max)
rho        <- array(c(1,.8,.8,1), dim = c(2,2,3))
rho.prime  <- abind(rho,rho, along = 4)

# Initial values
initsFunction <- function()list(
  
  n.start = matrix(2.5, 2, 33),
  rho.duck = array(c(1,.5,.5,1),dim = c(2,2,33)), 
  sig.duck = matrix(1,2, 33),
  sigma.obs.strata= matrix(1,2, 33),
  beta.indirect = matrix(0,2,2),
  beta.winter = c(0,0),
  
  rho = diag(6), sig.dem = runif(6,0,2),
  sig.sex = matrix(runif(4, 0, 1), 2,2),
  mu.kappa  = array(runif(8,-4,-1.5),dim = c(2,2,2)),
  mu.sur    = array(runif(12,-4,-3),dim = c(2,2,3)),
  beta.surv = array(rnorm(2*2*4, 0, .1), dim = c(2,2,4)), 
  beta.surv3 = matrix(rnorm(4,0,.1),2,2),
  cr        = rbeta(1,25,100),
  beta.np   = matrix(rnorm(6,0,1),2,3),
  nu.theta  = matrix(rlogis(4,0,1),2,2),
  nu.omega  = rlogis(2,0,1),
  sig.beta  = runif(2,0,1),
  global.mort = c(-3, -2, -5), sig.gm = c(2,1,2),
  global.kappa = c(-1.5, - 3), sig.gk = c(.5,.5),
  mu.trend  = rnorm(1, 0, .1),
  sig.trend = runif(1,0,.1),
  b0        = rnorm(1,0, 1), 
  sig.b0 = runif(1,0,.1),
  beta.trend.area=  rnorm(1,0, .1),
  beta.trend=  rnorm(24,0, .1),
  beta.sst  = rnorm(24,0, 1),
  beta.trend.area = rnorm(1,0, 1),
  beta.area = rnorm(1,0, 1),
  phi = rbeta(24, 10,10),
  sig = rexp(1,1),
  tau = rexp(1,1),
  decay.pond = rexp(1,1),
  sig.obs.chb = runif(2, 0, 2),
  
  mu.sst = 0, sig.sst =1,
  tau.car = c(.1,.1), sig.sal = c(1,1), car_lgm = rep(0, length(table(df$Unique1))),car_chb = rep(0, length(table(df1$Unique1))),
  eps.year = matrix(0, nrow = 59,2), eps.month = matrix(0,2,12),
  beta.sal = matrix(0,2,4), beta.flow = 0, beta.depth = c(0,0,0), sig.month = c(1,1), sig.year = c(1,1),beta.inter= c(0,0),
  rho.sal = diag(2),
  
  beta.spline = matrix(0, 2, ncol(BM1)), 
  lN = cbind(lny, matrix(rowMeans(lny), 2, 4)),
  y = y_inits,
  N = exp(lny),N.start = c(13,13),
  rho.sex = rho.sex,rho.sal = matrix(c(1,.5,.5,1),2,2),
  pond = log(rbind(ponds,pond.inits ,pond.inits ,pond.inits ,pond.inits)),
  mu.pond = t(log(rbind(ponds,pond.inits ,pond.inits ,pond.inits ,pond.inits ))),
  
  p.sig = 500,
  sigma.obs = runif(2,0,.5))

# Parameters monitored
pars <- c("mean.sw","mean.sb",'mean.sj', 'SJ','SA')

inits <- initsFunction() 

rm(list=setdiff(ls(), c('dat','constants','pars','inits')))
