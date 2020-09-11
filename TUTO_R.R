## ----message=FALSE, warning=FALSE, echo=FALSE-------------
library(tidyverse)
theme_set(theme_classic())

## ----include=FALSE,eval=FALSE-----------------------------
#  # automatically create a bib database for R packages
#  knitr::write_bib(c(
#    .packages(), 'bookdown', 'knitr', 'rmarkdown'
#  ), 'packages.bib')

## ---- eval=FALSE, include=TRUE----------------------------
#  install.packages(package.name)
#  library(packages.name)

## ----eval=FALSE,include=TRUE------------------------------
#  iris %>% summarize(mean_Petal=mean(Petal.Length))

## ----eval=FALSE,echo=correct------------------------------
#  install.packages("tidyverse")

## ----teacher=correct--------------------------------------
library(tidyverse)
iris %>% summarize(mean_Petal=mean(Petal.Length))

## ----echo=TRUE,eval=FALSE---------------------------------
#  x <- seq(-2*pi,2*pi,by=0.01)
#  y <- cos(x)
#  plot(x,y,type="l")

## ----eval=FALSE-------------------------------------------
#  getwd()

## ---------------------------------------------------------
b <- 41.3  # assigne la valeur 41.3 à l'objet b
x <- b     # b est assigné à x
x = b      # b est assigné à x
b -> x     # b est assigné à x
is.numeric(b)
mode(b)

## ---------------------------------------------------------
x <- "La mort"
y <- "aux trousses"
paste(x,y)
is.character(x)

## ---------------------------------------------------------
V1 <- factor(c("less20years","more50years","less20years","more50years","less20years"))
V1
levels(V1)
levels(V1) <- c("Young","Old")
V1

## ---------------------------------------------------------
x <- TRUE
is.logical(x)
mode(x)
a <- 1
a==1
a!=1
a<0
a>0

## ---------------------------------------------------------
x <- c(1.2,5,9,11)
x

## ---------------------------------------------------------
1:5

## ---------------------------------------------------------
seq(1,10,by=2)
seq(0,1,length=10)

## ---------------------------------------------------------
rep(1,4)
rep(c(1,3),each=3)

## ---------------------------------------------------------
x <- c("A","B","C")
x <- rep("A",5)
paste("X",1:5,sep="-")
substr("statistician",5,9)

## ---------------------------------------------------------
x <- c(-4,-3,1,3,5,8,0)
x[2]
x[c(2,5)]
x>0
x[x>0]

## ---------------------------------------------------------
x <- seq(-10,10,by=2)
y <- 1:length(x)
x+y
x*y
z <- x>0
x*z

## ----teacher=correct--------------------------------------
x <- c(1,3,8,9,11)
mean(x)
sum(x)
median(x)
var(x)

## ----teacher=correct--------------------------------------
rep(1:5,3)
rep(1:5,each=3)
rep(1:4,2:5)

## ----teacher=correct--------------------------------------
paste("A",0:10,")",sep="")

## ----teacher=correct--------------------------------------
index_q <- which(letters=="q")
paste(letters[1:index_q],1:index_q,sep="")

## ---------------------------------------------------------
m <- matrix(1:4,ncol=2)
m
m <- matrix(1:4,nrow=2)
m
m <- matrix(1:4,nrow=2,byrow=TRUE)
dim(m)

## ---------------------------------------------------------
m[2,1]

## ---------------------------------------------------------
m[1,] #première ligne
m[,2] #deuxième colonne

## ---------------------------------------------------------
det(m) #déterminant
solve(m) #inverse
t(m) #transposé
n <- matrix(5:8,nrow=2)
m+n
m*n #attention : produit de Hadamart
m%*%n #Produit matriciel
eigen(m) #Décomposition en valeurs propres

## ---------------------------------------------------------
mylist <- list(vector=rep(1:5),mat=matrix(1:8,nrow=2))
mylist
length(mylist)

## ---------------------------------------------------------
mylist[[1]]

## ---------------------------------------------------------
mylist$mat
mylist[["mat"]]

## ---------------------------------------------------------
name <- c("Paul","Mary","Steven","Charlotte","Peter")
sex <- factor(c("M","F","M","F","M"))
size <- c(180,165,168,170,175)
data <- data.frame(name,sex,size)
summary(data)

## ---------------------------------------------------------
data[2,3]
data[,2]
data$sex

## ---------------------------------------------------------
summary(data)
summary(1:10)

## ---------------------------------------------------------
x <- c(1,8,5,4)
sort(x)
order(x)

## ---------------------------------------------------------
V1 <- 1:10
V2 <- seq(-20,25,length=10)
df <- data.frame(V1,V2)
apply(df,1,mean)
apply(df,2,sum)

## ----teacher=correct--------------------------------------
mat <- matrix(c(1,0,3,4,5,5,0,4,5,6,3,4,0,1,3,2),ncol=4)
rownames(mat) <- paste("row-",1:4,sep="")
colnames(mat) <- paste("column ",1:4)
mat

## ----teacher=correct--------------------------------------
diag(mat)

## ----teacher=correct--------------------------------------
mat[1:2,]

## ----teacher=correct--------------------------------------
mat[,3:4]

## ----teacher=correct--------------------------------------
det(mat)
solve(mat)

## ---------------------------------------------------------
data(iris)
head(iris)

## ----teacher=correct--------------------------------------
mean(iris$Sepal.Width)
mean(iris$Petal.Length)
var(iris$Sepal.Width)
var(iris$Petal.Length)

## ----teacher=correct--------------------------------------
test <- iris$Species=="versicolor"
iris2 <- iris[test,]

## ----teacher=correct--------------------------------------
ord <- order(iris2$Sepal.Length,decreasing=TRUE)
iris3 <- iris2[ord,]

## ----teacher=correct--------------------------------------
mean(iris[iris$Species=="versicolor","Sepal.Length"])
mean(iris[iris$Species=="virginica","Sepal.Length"])
mean(iris[iris$Species=="setosa","Sepal.Length"])

## ----teacher=correct--------------------------------------
iris$sum.petal <- iris$Petal.Length+iris$Petal.Width
head(iris)

## ---------------------------------------------------------
library(lattice)
data("ethanol")

## ----teacher=correct--------------------------------------
summary(ethanol)
apply(ethanol,2,mean)

## ----teacher=correct--------------------------------------
quantile(ethanol$NOx,probs=c(0.25,0.5,0.75))
quantile(ethanol$C,probs=c(0.25,0.5,0.75))
quantile(ethanol$E,probs=c(0.25,0.5,0.75))
#ou mieux
apply(ethanol,2,quantile,probs=c(0.25,0.5,0.75))

## ----teacher=correct--------------------------------------
apply(ethanol,2,quantile,probs=seq(0.1,0.9,by=0.1))

## ---------------------------------------------------------
data("presidents")
df <- matrix(presidents,ncol=4,byrow=T)

## ----teacher=correct--------------------------------------
any(is.na(df[20,]))

## ----teacher=correct--------------------------------------
which(apply(is.na(df),1,any))

## ----teacher=correct--------------------------------------
ind_sup <- which(apply(is.na(df),1,any))
df1 <- df[-ind_sup,]
summary(df1)

## ----teacher=correct--------------------------------------
df2 <- na.omit(df)
all(df1==df2)

## ---------------------------------------------------------
path <- file.path("data/", "piscines.csv") #premier : répertoire, deuxième : fichier
piscines <- read.csv(path)
class(piscines)
summary(piscines)

## ---------------------------------------------------------
library(readr)
piscines <- read_csv("data/piscines.csv")
summary(piscines)
class(piscines)

## ----teacher=correct--------------------------------------
data1 <- read.csv("data/mydata.csv")
summary(data1)

## ----teacher=correct--------------------------------------
data1 <- read.csv("data/mydata.csv",sep=";",dec=c("."),row.names = 1)
summary(data1)

## ----teacher=correct--------------------------------------
data2 <- read.csv("data/mydata2.csv")
summary(data2)

## ----teacher=correct--------------------------------------
data2 <- read.csv("data/mydata2.csv",sep="",na.strings = ".")
summary(data2)

## ----teacher=correct--------------------------------------
data22 <- data2

## ----teacher=correct--------------------------------------
levels(data2$sex) <- c("woman","man")
summary(data2)

## ----teacher=correct--------------------------------------
library(tidyverse)
data22$sex <- recode_factor(data2$sex,"F"="woman","M"="man")
summary(data22)

## ---------------------------------------------------------
df1 <- tibble(name=c("Mary","Peter","John","July"),age=c(18,25,21,43))
df2 <- tibble(name=c("Zac","Julian"),age=c(23,48))
df3 <- tibble(size=c(154,178,182,134,142),name1=c("Peter","Mary","July","John","stef"))
df1
df2
df3

## ----teacher=correct--------------------------------------
df <- bind_rows(df1,df2)
mean(df$age)

## ----teacher=correct--------------------------------------
a1 <- full_join(df,df3,by=c("name"="name1"))
a1

## ----teacher=correct--------------------------------------
a2 <- inner_join(df,df3,by=c("name"="name1"))
a2

## ---------------------------------------------------------
piscines[piscines$Longitude>153,c("Longitude","Latitude")]

## ---------------------------------------------------------
library(tidyverse) #ou library(dplyr)
piscines %>% select(Longitude,Latitude) %>% filter(Longitude>153)

## ---- eval=FALSE, include=TRUE----------------------------
#  select(df, VAR1, VAR2, ...)

## ---------------------------------------------------------
coord <- select(piscines, Latitude, Longitude)
head(piscines, n=2)
head(coord, n=2)

## ---------------------------------------------------------
coord <- select(piscines, ends_with("tude"))
head(coord, n=2)

## ---- eval=FALSE, include=TRUE----------------------------
#  mutate(df, NEW.VAR = expression(VAR1, VAR2, ...))

## ---------------------------------------------------------
df <- mutate(piscines, phrase=paste("Swimming pool", Name, "is located at the address", Address))
select(df,phrase)

## ---------------------------------------------------------
mutate(piscines,
       phrase = paste("Swimming pool", Name, "is located at the address", Address),
       unused = Longitude + Latitude
)

## ---- eval=FALSE, include=TRUE----------------------------
#  filter(df, TEST)

## ---------------------------------------------------------
p1 <- filter(piscines, Longitude>153.02)
select(p1,Longitude)

## ---------------------------------------------------------
df <- filter(piscines, !grepl("Pool", Name))
select(df,Name)

## ---------------------------------------------------------
p2 <- filter(piscines, Longitude>153.02 | Latitude < -27.488)
select(p2, Longitude, Latitude)

## ---------------------------------------------------------
slice(piscines,5:8)

## ---- eval=FALSE, include=TRUE----------------------------
#  arrange(df, VAR) #tri croissant

## ---- eval=FALSE, include=TRUE----------------------------
#  arrange(df, desc(VAR)) #tri décroissant

## ---------------------------------------------------------
arrange(piscines, Longitude)

## ---------------------------------------------------------
arrange(piscines, desc(Longitude))

## ---------------------------------------------------------
summarise(piscines,
          mean_long = mean(Longitude),
          med_lat = median(Latitude),
          min_lat = min(Latitude),
          sum_long = sum(Longitude)
)

## ---------------------------------------------------------
summarise(piscines,n())
summarise(piscines,last(Longitude))

## ---------------------------------------------------------
summarise_at(piscines,3:4,mean)

## ---------------------------------------------------------
lat_mean <- piscines %>% summarise(mean(Latitude))
pisc1 <- piscines %>% mutate(lat_dis=factor(Latitude>as.numeric(lat_mean)))
levels(pisc1$lat_dis) <- c("Low","High")

## ---------------------------------------------------------
summarise(group_by(pisc1,lat_dis),mean_long=mean(Longitude))

## ----eval=FALSE-------------------------------------------
#  pisc1

## ----eval=FALSE-------------------------------------------
#  pisc1 %>% group_by(lat_dis)

## ---------------------------------------------------------
pisc1 %>% group_by(lat_dis) %>% summarise(mean_long=mean(Longitude))

## ---------------------------------------------------------
mean(1:10)
1:10 %>% mean()

## ---------------------------------------------------------
iris <- iris %>% as_tibble()

## ----teacher=correct--------------------------------------
iris %>% select(Petal.Width,Species)

## ----teacher=correct--------------------------------------
iris %>% filter(Species=="versicolor" | Species=="virginica")

## ----teacher=correct--------------------------------------
iris %>% filter(Species=="setosa") %>% summarise(n())

## ----teacher=correct--------------------------------------
iris %>% filter(Species=="versicolor") %>%
  summarise(Mean_PW=mean(Petal.Width))

## ----teacher=correct--------------------------------------
iris1 <- iris
iris1 %>% mutate(Sum_Petal=Petal.Width+Sepal.Width)

## ----teacher=correct--------------------------------------
iris %>% group_by(Species) %>%
  summarise(mean_PL=mean(Petal.Length),var_PL=var(Petal.Length)) %>%
  mutate(var_PL=round(var_PL,3))

## ---------------------------------------------------------
library(hflights)
hflights <- as_tibble(hflights)

## ---------------------------------------------------------
lut1 <- c("AA" = "American", "AS" = "Alaska", "B6" = "JetBlue", "CO" = "Continental",
         "DL" = "Delta", "OO" = "SkyWest", "UA" = "United", "US" = "US_Airways", 
         "WN" = "Southwest", "EV" = "Atlantic_Southeast", "F9" = "Frontier", 
         "FL" = "AirTran", "MQ" = "American_Eagle", "XE" = "ExpressJet", "YV" = "Mesa")

## ---------------------------------------------------------
lut2 <- c("A" = "carrier", "B" = "weather", "C" = "FFA", "D" = "security", "E" = "not cancelled")

## ---------------------------------------------------------
hflights1 <- hflights
hflights1$UniqueCarrier <- lut1[hflights1$UniqueCarrier]
hflights1$CancellationCode[hflights1$CancellationCode==""] <- "Z"
hflights1$CancellationCode <- lut2[hflights1$CancellationCode]

## ----teacher=correct--------------------------------------
ind <- match(c("Origin","Cancelled"),names(hflights1))
hflights1 %>% select(seq(ind[1],ind[2]))
#ou
hflights1 %>% select(Origin:Cancelled) 

## ----teacher=correct--------------------------------------
hflights1 %>% select(contains("Time"),contains("Delay"))

## ----teacher=correct--------------------------------------
hflights2 <- hflights1 %>% mutate(ActualGroundTime=ActualElapsedTime-AirTime)
head(hflights2)

## ----teacher=correct--------------------------------------
hflights3 <- hflights2 %>% 
  mutate(AverageSpeed=Distance/AirTime) %>%
  arrange(desc(AverageSpeed))

## ----teacher=correct--------------------------------------
filter(hflights3,Dest=="JFK")

## ----teacher=correct--------------------------------------
hflights3 %>% filter(Dest=="JFK") %>% summarise(numb_to_JFK=n())

## ----teacher=correct--------------------------------------
hflights1 %>% summarize(n_flights=n(),n_dest=n_distinct(Dest),n_carrier=n_distinct(UniqueCarrier))

## ----teacher=correct--------------------------------------
hflights1 %>% filter(UniqueCarrier=="American") %>% 
  summarize(n_fligths_Am=n(),n_can_Am=sum(Cancelled),
        mean_ArrDelay_am=mean(ArrDelay,na.rm=TRUE))

## ----teacher=correct--------------------------------------
hflights1 %>% group_by(UniqueCarrier) %>%
  summarise(n_flights=n(),mean_AirTime=mean(AirTime,na.rm=TRUE))

## ----teacher=correct--------------------------------------
hflights1 %>% 
  group_by(UniqueCarrier) %>%
  filter(!is.na(DepDelay) & DepDelay>0) %>%
  summarise(meanDepDelay = mean(DepDelay)) %>%
  arrange(meanDepDelay)

## ----teacher=correct--------------------------------------
FrenchOpen_men_2013 <- read_csv("data/FrenchOpen-men-2013.csv")
RG2013 <- FrenchOpen_men_2013
RG2013

## ----teacher=correct--------------------------------------
RG2013 %>% filter(Player1=="Roger Federer" | Player2=="Roger Federer") %>%
  select(Player1,Player2)

## ----teacher=correct--------------------------------------
RG2013 %>% filter(Round==6) %>% select(Player1,Player2)

## ----teacher=correct--------------------------------------
RG2013 %>% mutate(nb_points=TPW.1+TPW.2) %>% select(nb_points) %>% summarize_all(mean)

## ----teacher=correct--------------------------------------
RG2013 %>% mutate(nb_aces=ACE.1+ACE.2) %>% summarize(mean_aces=mean(nb_aces))

## ----teacher=correct--------------------------------------
RG2013 %>% group_by(Round) %>% mutate(nb_aces=ACE.1+ACE.2) %>%
  summarize(mean_aces=mean(nb_aces))

## ----teacher=correct--------------------------------------
RG2013 %>% mutate(nb_df=DBF.1+DBF.2) %>% 
  summarize(nb_dbfaults=sum(nb_df,na.rm=TRUE))

## ----teacher=correct--------------------------------------
WIMB2013 <- read_csv("data/Wimbledon-men-2013.csv")
WIMB2013

## ----teacher=correct--------------------------------------
RG_WIMB2013 <- bind_rows("RG"=RG2013,"WIMB"=WIMB2013,.id="Tournament")

## ----teacher=correct--------------------------------------
RG_WIMB2013  %>% filter(Player1=="Roger Federer" | 
                      Player2=="Roger Federer" |
                      Player1=="R.Federer" | 
                      Player2=="R.Federer") 

## ----teacher=correct--------------------------------------
RG_WIMB2013  %>% filter(grepl("Federer",Player2) | grepl("Federer",Player2))

## ----teacher=correct--------------------------------------
RG_WIMB2013 %>% group_by(Tournament,Round) %>% 
  mutate(nb_aces=ACE.1+ACE.2) %>% summarize(mean_ace=mean(nb_aces))

## ----teacher=correct--------------------------------------
RG_WIMB2013 %>% group_by(Tournament,Round) %>% 
  mutate(nb_aces=ACE.1+ACE.2) %>% 
  summarize(mean_ace=mean(nb_aces)) %>%
  pivot_wider(names_from = "Round",values_from = "mean_ace")

## ---------------------------------------------------------
df <- read_delim("data/tauxchomage.csv",delim=";") %>% select(-1)
df

## ---------------------------------------------------------
df1 <- df %>% pivot_longer(-NOM_DPT,names_to="Année",values_to="TCHOM") %>% 
  mutate(Année=fct_recode(Année,"2001"="TCHOMB1T01","2006"="TCHOMB1T06","2011"="TCHOMB1T11"))
df1

## ---------------------------------------------------------
ggplot(df1)+aes(x=Année,y=TCHOM)+geom_boxplot()

## ---------------------------------------------------------
df1 %>% pivot_wider(names_from="Année",values_from="TCHOM")

## ---------------------------------------------------------
df <- tibble(date=as.Date(c("01/03/2015","05/18/2017",
          "09/14/2018"),"%m/%d/%Y"),temp=c(18,21,15))
df
df1 <- df %>% separate(date,into = c("year","month","day"))
df1

## ---------------------------------------------------------
df1 %>% unite(date,year,month,day,sep="/")

## ---- eval=FALSE, echo=TRUE-------------------------------
#  demo(graphics)

## ---------------------------------------------------------
x <- seq(-2*pi,2*pi,by=0.05)
y <- sin(x)
plot(x,y) #points (par défaut)
plot(x,y,type="l") #représentation sous forme de ligne

## ---------------------------------------------------------
ozone <- read.table("data/ozone.txt")
summary(ozone)

## ---------------------------------------------------------
plot(ozone[,"T12"],ozone[,"maxO3"])

## ---------------------------------------------------------
plot(maxO3~T12,data=ozone)

## ---------------------------------------------------------
plot(ozone[,"T12"],ozone[,"maxO3"],xlab="T12",ylab="maxO3")

## ---------------------------------------------------------
hist(ozone$maxO3,main="Histogram")
barplot(table(ozone$vent)/nrow(ozone),col="blue")
boxplot(maxO3~vent,data=ozone)

## ----message=FALSE, warning=FALSE-------------------------
library(rAmCharts)
amHist(ozone$maxO3)
amPlot(ozone,col=c("T9","T12"))
amBoxplot(maxO3~vent,data=ozone)

## ----test-a,echo=cor,eval=cor-----------------------------
x <- seq(0,2*pi,length=1000)
plot(x,sin(x),type="l")
title("Représentation de la fonction sinus")

## ----echo=cor,eval=cor------------------------------------
x <- seq(-4,4,by=0.01)
plot(x,dnorm(x),type="l")
abline(v=0,lty=2)
lines(x,dt(x,5),col=2)
lines(x,dt(x,30),col=3)
legend("topleft",legend=c("normal","Student(5)","Student(30)"),
   col=1:3,lty=1)

## ----eval=cor,echo=cor------------------------------------
taches <- read.table("data/taches_solaires.csv",sep=";",header=TRUE,dec=",")

## ----eval=cor,echo=cor------------------------------------
library(tidyverse)
periode <- cut_interval(taches$annee,n=8)

## ----eval=cor,echo=TRUE-----------------------------------
couleurs <- c("yellow", "magenta", "orange", "cyan",
          "grey", "red", "green", "blue")

## ----eval=cor,echo=cor------------------------------------
levels(periode) <- couleurs

## ----eval=cor,echo=cor------------------------------------
coordx <- seq(along=taches[,1])

## ----eval=cor,echo=cor------------------------------------
plot(coordx,taches[,1],xlab="Temps",ylab="Nombre de taches",
 col=periode,type="p",pch="+")

## ----eval=cor,echo=cor------------------------------------
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(maxO3~T12,data=ozone)
hist(ozone$T12)
boxplot(ozone$maxO3)

## ----message=FALSE, warning=FALSE-------------------------
library(tidyverse)
set.seed(1234)
diamonds2 <- diamonds[sample(nrow(diamonds),5000),] 
summary(diamonds2)
help(diamonds)

## ---------------------------------------------------------
plot(price~carat,data=diamonds2)

## ---------------------------------------------------------
ggplot(diamonds2) #rien
ggplot(diamonds2)+aes(x=carat,y=price) #rien
ggplot(diamonds2)+aes(x=carat,y=price)+geom_point() #bon

## ----echo=cor,eval=cor------------------------------------
ggplot(diamonds2)+aes(x=carat)+geom_histogram()
ggplot(diamonds2)+aes(x=carat)+geom_histogram(bins=10)
ggplot(diamonds2)+aes(x=cut)+geom_bar()

## ----eval=FALSE-------------------------------------------
#  ggplot(diamonds2)+aes(x=carat,y=price)

## ----eval=FALSE-------------------------------------------
#  ggplot(diamonds2)+aes(x=carat,y=price,color=cut)

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price,color=cut)+geom_point()

## ----echo=cor,eval=cor------------------------------------
ggplot(diamonds2)+aes(x=cut)+geom_bar(fill="blue")

## ----echo=cor,eval=cor------------------------------------
ggplot(diamonds2)+aes(x=cut,fill=cut)+geom_bar()

## ----echo=cor,eval=cor------------------------------------
ggplot(diamonds2)+aes(x=cut)+geom_bar(fill=c("blue","red","green","yellow","black"))

## ---------------------------------------------------------
D <- data.frame(X=seq(-2*pi,2*pi,by=0.01))
ggplot(D)+aes(x=X,y=sin(X))+geom_line()

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=price)+geom_histogram(bins=40)

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=price,y=..density..)+geom_histogram(bins=40)

## ----eval=FALSE,echo=TRUE---------------------------------
#  ggplot(diamonds2)+aes(x=price,y=..density..)+stat_bin()

## ----echo=cor,eval=cor------------------------------------
X <- data.frame(X1=c("red","blue","green","black"),prob=c(0.3,0.2,0.4,0.1))
ggplot(X)+aes(x=X1,y=prob,fill=X1)+geom_bar(stat="identity")+
  labs(fill="Couleur")+xlab("")

## ----echo=cor,eval=cor------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price)+geom_smooth(method="loess")
ggplot(diamonds2)+aes(x=carat,y=price)+stat_smooth(method="loess")


## ----echo=cor,eval=cor------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price)+geom_smooth(method="loess",linetype="dotted")
ggplot(diamonds2)+aes(x=carat,y=price)+stat_smooth(method="loess",geom="point")

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price,color=cut)+geom_point()+
  scale_color_manual(values=c("Fair"="black","Good"="yellow",
                              "Very Good"="blue","Premium"="red","Ideal"="green"))

## ---------------------------------------------------------
p1 <- ggplot(diamonds2)+aes(x=cut)+geom_bar(aes(fill=cut))
p1

## ---------------------------------------------------------
p1+scale_fill_brewer(palette="Purples")

## ---------------------------------------------------------
p2 <- ggplot(diamonds2)+aes(x=carat,y=price)+geom_point(aes(color=depth))
p2

## ---------------------------------------------------------
p2+scale_color_gradient(low="red",high="yellow")

## ---------------------------------------------------------
p2+scale_x_continuous(breaks=seq(0.5,3,by=0.5))+
  scale_y_continuous(name="prix")+
  scale_color_gradient("Profondeur")

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price,group=cut)+
  geom_smooth(method="loess")

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price)+
  geom_smooth(method="loess")+facet_wrap(~cut)
ggplot(diamonds2)+aes(x=carat,y=price)+
  geom_smooth(method="loess")+facet_wrap(~cut,nrow=1)

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price)+geom_point()+
  geom_smooth(method="lm")+facet_grid(color~cut)
ggplot(diamonds2)+aes(x=carat,y=price)+geom_point()+
  geom_smooth(method="lm")+facet_wrap(color~cut)

## ----echo=TRUE,eval=FALSE---------------------------------
#  ggplot()+aes()+geom_()+scale_()

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price)+geom_point()
ggplot(diamonds2,aes(x=carat,y=price))+geom_point()
ggplot(diamonds2)+geom_point(aes(x=carat,y=price))

## ---------------------------------------------------------
X <- seq(-2*pi,2*pi,by=0.001)
Y1 <- cos(X)
Y2 <- sin(X)
donnees1 <- data.frame(X,Y1)
donnees2 <- data.frame(X,Y2)
ggplot(donnees1)+geom_line(aes(x=X,y=Y1))+
  geom_line(data=donnees2,aes(x=X,y=Y2),color="red")

## ---------------------------------------------------------
p <- ggplot(diamonds2)+aes(x=carat,y=price,color=cut)+geom_point()
p+theme_bw()
p+theme_classic()
p+theme_grey()
p+theme_bw()

## ----echo=cor,eval=cor------------------------------------
X <- seq(-2*pi,2*pi,by=0.001)
Y1 <- cos(X)
Y2 <- sin(X)
donnees1 <- data.frame(X,Y1)
donnees2 <- data.frame(X,Y2)
ggplot(donnees1)+geom_line(aes(x=X,y=Y1))+
  geom_line(data=donnees2,aes(x=X,y=Y2),color="red")

## ----echo=cor,eval=cor------------------------------------
donnees <- data.frame(X,Y1,Y2)
ggplot(donnees)+aes(x=X,y=Y1)+geom_line()+
  geom_line(aes(y=Y2),color="red")
#ou pour la légende
ggplot(donnees)+aes(x=X,y=Y1)+geom_line(aes(color="cos"))+
  geom_line(aes(y=Y2,color="sin"))+labs(color="Fonction")

## ----echo=cor,eval=cor------------------------------------
df <- data.frame(X,cos=Y1,sin=Y2)
df1 <- df %>% pivot_longer(cols=c(cos,sin),
                       names_to = "Fonction",
                       values_to = "value")
#ou
df1 <- df %>% pivot_longer(cols=-X,
                       names_to = "Fonction",
                       values_to = "value")
ggplot(df1)+aes(x=X,y=value,color=Fonction)+geom_line()

## ----echo=cor,eval=cor------------------------------------
ggplot(df1)+aes(x=X,y=value)+geom_line()+facet_wrap(~Fonction)

## ----echo=cor,eval=cor------------------------------------
library(gridExtra)
p1 <- ggplot(donnees1)+aes(x=X,y=Y1)+geom_line()
p2 <- ggplot(donnees2)+aes(x=X,y=Y2)+geom_line()
grid.arrange(p1,p2,nrow=1)

## ---------------------------------------------------------
data(mtcars)
summary(mtcars)

## ----echo=cor,eval=cor------------------------------------
ggplot(mtcars)+aes(x=mpg)+geom_histogram()
ggplot(mtcars)+aes(x=mpg)+geom_histogram(bins=10)

## ----echo=cor,eval=cor------------------------------------
ggplot(mtcars)+aes(x=mpg,y=..density..)+geom_histogram(bins=10)

## ----echo=cor,eval=cor------------------------------------
ggplot(mtcars)+aes(x=cyl)+geom_bar()

## ----echo=cor,eval=cor------------------------------------
ggplot(mtcars)+aes(x=disp,y=mpg,color=cyl)+geom_point()
ggplot(mtcars)+aes(x=disp,y=mpg,color=as.factor(cyl))+geom_point()+labs(color="cyl")

## ----echo=cor,eval=cor------------------------------------
ggplot(mtcars)+aes(x=disp,y=mpg,color=as.factor(cyl))+geom_point()+
  geom_smooth(method="lm")+labs(color="cyl")

## ----echo=cor,eval=cor------------------------------------
n <- 100
X <- runif(n)
eps <- rnorm(n,sd=0.2)
Y <- 3+X+eps
D <- data.frame(X,Y)

## ----echo=cor,eval=cor------------------------------------
model <- lm(Y~.,data=D)
co <- coef(model)
D$fit <- predict(model)
co <- coef(lm(Y~.,data=D))
ggplot(D)+aes(x=X,y=Y)+geom_point()+
  geom_abline(slope=co[2],intercept=co[1],color="blue")

## ----echo=cor,eval=cor------------------------------------
ggplot(D)+aes(x=X,y=Y)+geom_point()+geom_smooth(method="lm")

## ----echo=cor,eval=cor------------------------------------
ggplot(D)+aes(x=X,y=Y)+geom_point()+geom_smooth(method="lm")+
  geom_segment(aes(xend=X,yend=fit))

## ----echo=FALSE,eval=TRUE---------------------------------
ggplot(data=diamonds) + geom_boxplot(aes(x=cut,y=carat,fill=cut)) 
ggplot(data=diamonds) + geom_boxplot(aes(x=cut,y=carat,fill=cut))+coord_flip()
ggplot(data=diamonds) + geom_density(aes(x=carat,y=..density..)) +  facet_grid(cut~.)

## ----echo=cor,eval=FALSE----------------------------------
#  ggplot(data=diamonds) + geom_boxplot(aes(x=cut,y=carat,fill=cut))
#  ggplot(data=diamonds) + geom_boxplot(aes(x=cut,y=carat,fill=cut))+coord_flip()
#  ggplot(data=diamonds) + geom_density(aes(x=carat,y=..density..)) +  facet_grid(cut~.)

## ----echo=cor,eval=TRUE-----------------------------------
Q1 <- diamonds %>% group_by(cut) %>% 
  summarize(q1=quantile(carat,c(0.25)),q2=quantile(carat,c(0.5)),
        q3=quantile(carat,c(0.75)))
quantildf <- Q1%>% gather(key="alpha",value="quantiles",-cut)
ggplot(data=diamonds) + geom_density(aes(x=carat,y=..density..))+
  facet_grid(cut~.) +
  geom_vline(data=quantildf,aes(xintercept=quantiles),col=alpha("black",1/2))

## ----message=FALSE, warning=FALSE,echo=FALSE,eval=TRUE----
library(ggstance)
ggplot(data=diamonds) +
  geom_boxploth(data=diamonds,aes(y=-0.5,x=carat,fill=cut)) +
  geom_density(aes(x=carat,y=..density..)) +  facet_grid(cut~.) +
  geom_vline(data=quantildf,aes(xintercept=quantiles),col=alpha("black",1/2))+
  ylab("")

## ----message=FALSE, warning=FALSE,echo=cor,eval=FALSE-----
#  library(ggstance)
#  ggplot(data=diamonds) +
#    geom_boxploth(data=diamonds,aes(y=-0.5,x=carat,fill=cut)) +
#    geom_density(aes(x=carat,y=..density..)) +  facet_grid(cut~.) +
#    geom_vline(data=quantildf,aes(xintercept=quantiles),col=alpha("black",1/2))+
#    ylab("")

## ----message=FALSE, warning=FALSE-------------------------
library(leaflet)
leaflet() %>% addTiles()

## ---------------------------------------------------------
Paris <- c(2.351462,48.8567)
m2 <- leaflet() %>% setView(lng = Paris[1], lat = Paris[2], zoom = 12) %>% 
  addTiles()
m2 %>% addProviderTiles("Stamen.Toner")
m2 %>% addProviderTiles("Wikimedia")
m2 %>% addProviderTiles("Esri.NatGeoWorldMap")
m2 %>%
  addProviderTiles("Stamen.Watercolor") %>%
  addProviderTiles("Stamen.TonerHybrid")

## ---------------------------------------------------------
data(quakes)
leaflet(data = quakes[1:20,]) %>% addTiles() %>%
  addMarkers(~long, ~lat, popup = ~as.character(mag))

## ---------------------------------------------------------
content <- paste(sep = "<br/>",
  "<b><a href='http://www.samurainoodle.com'>Samurai Noodle</a></b>",
  "606 5th Ave. S",
  "Seattle, WA 98138"
)

leaflet() %>% addTiles() %>%
  addPopups(-122.327298, 47.597131, content,
    options = popupOptions(closeButton = FALSE)
  )

## ----message=FALSE, warning=FALSE-------------------------
if (!(require(jsonlite))) install.packages("jsonlite")
mygeocode <- function(adresses){
# adresses est un vecteur contenant toutes les adresses sous forme de chaine de caracteres
  nominatim_osm <- function(address = NULL){
    ## details: http://wiki.openstreetmap.org/wiki/Nominatim
    ## fonction nominatim_osm proposée par D.Kisler
    if(suppressWarnings(is.null(address)))  return(data.frame())
    tryCatch(
      d <- jsonlite::fromJSON(
        gsub('\\@addr\\@', gsub('\\s+', '\\%20', address),
             'http://nominatim.openstreetmap.org/search/@addr@?format=json&addressdetails=0&limit=1')
      ), error = function(c) return(data.frame())
    )
    if(length(d) == 0) return(data.frame())
    return(c(as.numeric(d$lon), as.numeric(d$lat)))
  }
  tableau <- t(sapply(adresses,nominatim_osm))
  colnames(tableau) <- c("lon","lat")
  return(tableau)
}


## ---------------------------------------------------------
mygeocode("Paris")

## ----teacher=correct--------------------------------------
R2 <- mygeocode("Universite Rennes 2 Villejean") %>% as_tibble()
info <- paste(sep = "<br/>",
  "<b><a href='https://www.univ-rennes2.fr'>Universite Rennes 2</a></b>",
  "Campus Villejean")


leaflet() %>% addTiles() %>%  
  addPopups(R2[1]$lon, R2[2]$lat, info,options = popupOptions(closeButton = FALSE))


## ----echo=FALSE,eval=FALSE--------------------------------
#  #Pour éviter les problèmes de changement
#  sta.Paris <- read_delim("https://opendata.paris.fr/explore/dataset/velib-disponibilite-en-temps-reel/download/?format=csv&timezone=Europe/Berlin&use_labels_for_header=true",delim=";")
#  write_csv(sta.Paris,path="data/sta.Paris.csv")

## ----echo=FALSE,eval=TRUE---------------------------------
sta.Paris <- read_csv("data/sta.Paris.csv")

## ---- echo=correct,eval=FALSE-----------------------------
#  lien <- "https://opendata.paris.fr/explore/dataset/velib-disponibilite-en-temps-reel/
#  download/?format=csv&timezone=Europe/Berlin&use_labels_for_header=true"
#  sta.Paris <- read_delim(lien,delim=";")

## ---- echo=correct,eval=TRUE------------------------------
sta.Paris1 <- sta.Paris %>% separate(`Coordonnées géographiques`,
                                 into=c("lat","lon"),sep=",") %>% 
  mutate(lat=as.numeric(lat),lon=as.numeric(lon))

## ----teacher=correct--------------------------------------
map.velib1 <- leaflet(data = sta.Paris1) %>% 
  addTiles() %>%
  addCircleMarkers(~ lon, ~ lat,radius=3,
stroke = FALSE, fillOpacity = 0.5,color="red"
  )

map.velib1

## ----teacher=correct--------------------------------------
map.velib2 <- leaflet(data = sta.Paris1) %>% 
  addTiles() %>% 
  addCircleMarkers(~ lon, ~ lat,radius=3,stroke = FALSE, 
               fillOpacity = 0.7,color="red", 
               popup = ~ sprintf("<b> Vélos dispos: %s</b>",
                                 as.character(`Nombre total vélos disponibles`)))

#ou sans sprintf

map.velib2 <- leaflet(data = sta.Paris1) %>% 
  addTiles() %>% 
  addCircleMarkers(~ lon, ~ lat,radius=3,stroke = FALSE, fillOpacity = 0.7,color="red", 
               popup = ~ paste("Vélos dispos :",
                               as.character(`Nombre total vélos disponibles`)))

map.velib2

## ----teacher=correct--------------------------------------
map.velib3 <- leaflet(data = sta.Paris1) %>% 
  addTiles() %>%
  addCircleMarkers(~ lon, ~ lat,radius=3,stroke = FALSE, 
               fillOpacity = 0.7,color="red", 
               popup = ~ paste(as.character(`Nom station`),", Vélos dispos :",
                               as.character(`Nombre total vélos disponibles`),
                               sep=""))

map.velib3

## ---------------------------------------------------------
ColorPal1 <- colorNumeric(scales::seq_gradient_pal(low = "#132B43", high = "#56B1F7",
                                               space = "Lab"), domain = c(0,1))
ColorPal2 <- colorNumeric(scales::seq_gradient_pal(low = "red", high = "black", 
                                               space = "Lab"), domain = c(0,1))

## ----teacher=correct--------------------------------------
map.velib4 <- leaflet(data = sta.Paris1) %>% 
  addTiles() %>%
  addCircleMarkers(~ lon, ~ lat,radius=3,stroke = FALSE, fillOpacity = 0.7,
               color=~ColorPal1(`Nombre total vélos disponibles`/
                                  `Capacité de la station`), 
               popup = ~ paste(as.character(`Nom station`),", Vélos dispos :",
                               as.character(`Nombre total vélos disponibles`),
                               sep=""))

map.velib4
map.velib5 <- leaflet(data = sta.Paris1) %>% 
  addTiles() %>%
  addCircleMarkers(~ lon, ~ lat,stroke = FALSE, fillOpacity = 0.7,
               color=~ColorPal2(`Nombre total vélos disponibles`/
                                  `Capacité de la station`),
               radius=~(`Nombre total vélos disponibles`/
                          `Capacité de la station`)*8,
               popup = ~ paste(as.character(`Nom station`),", Vélos dispos :",
                               as.character(`Nombre total vélos disponibles`),
                               sep=""))

map.velib5

## ----echo=correct,eval=TRUE-------------------------------
nom.station <- "Jussieu - Fossés Saint-Bernard"
local.station <- function(nom.station){
  df <- sta.Paris1 %>% filter(`Nom station`==nom.station)
  leaflet(data = sta.Paris1) %>% setView(lng=df$lon,lat=df$lat,zoom=15) %>%
addTiles() %>% 
addCircleMarkers(~ lon, ~ lat,stroke = FALSE, fillOpacity = 0.7,
                popup = ~ paste(as.character(`Nom station`),", Vélos dispos :",
                                as.character(`Nombre total vélos disponibles`),
                                sep="")) %>%
addMarkers(lng=df$lon,lat=df$lat,
           popup = ~ paste(nom.station,", Vélos dispos :",
                           as.character(df$`Nombre total vélos disponibles`),
                           sep=""),
           popupOptions = popupOptions(noHide = T))
}

## ---------------------------------------------------------
local.station("Jussieu - Fossés Saint-Bernard")
local.station("Gare Montparnasse - Arrivée")

## ----eval=FALSE-------------------------------------------
#  method(Y~X1+X3,data=df,...)

## ---------------------------------------------------------
n <- 1000
p <- 5
set.seed(1234)
X.mat <- matrix(rnorm(n*p),ncol=p)
eps <- rnorm(n,mean = 0,sd=0.5)
df <- data.frame(X.mat,eps)
df <- df %>% mutate(Y=X1+X2+X3+X4+X5+eps) %>% select(-eps)

## ----teacher=correct--------------------------------------
mod1 <- lm(Y~.,data=df)
coef(mod1)
summary(mod1)

## ---------------------------------------------------------
m <- 500
p <- 5
set.seed(12345)
X.mat <- matrix(rnorm(m*p),ncol=5)
eps <- rnorm(m,mean = 0,sd=0.5)
df.test <- data.frame(X.mat,eps)
df.test <- df.test %>% mutate(Y=X1+X2+X3+X4+X5+eps) %>% select(-eps)

## ----teacher=correct--------------------------------------
pred <- predict(mod1,newdata=df.test)
head(pred)

## ----teacher=correct--------------------------------------
pred.df <- data.frame(pred=pred,obs=df.test$Y)

## ----teacher=correct--------------------------------------
pred.df %>% summarize(MSE=mean((pred-obs)^2))

## ---------------------------------------------------------
n <- 1000
p <- 105
set.seed(1234)
X.mat <- matrix(rnorm(n*p),ncol=p)
eps <- rnorm(n,mean = 0,sd=0.5)
df <- data.frame(X.mat,eps)
df <- df %>% mutate(Y=X1+X2+X3+X4+X5+eps) %>% select(-eps)

## ----teacher=correct--------------------------------------
mod2 <- lm(Y~.,data=df)
summary(mod2)$coefficients %>% head()

## ----selection-step,teacher=correct,cache=TRUE------------
mod.step <- step(mod2,direction=c("backward"),k=log(n),trace=0)
summary(mod.step)

## ---------------------------------------------------------
m <- 300
p <- 105
set.seed(12345)
X.mat <- matrix(rnorm(m*p),ncol=p)
eps <- rnorm(m,mean = 0,sd=0.5)
df.test <- data.frame(X.mat,eps)
df.test <- df.test %>% mutate(Y=X1+X2+X3+X4+X5+eps) %>% select(-eps)

## ----teacher=correct--------------------------------------
p.full <- predict(mod2,newdata=df.test)
p.step <- predict(mod.step,newdata=df.test)
pred.df <- tibble(full=p.full,step=p.step,obs=df.test$Y)

## ----teacher=correct--------------------------------------
pred.df %>% summarize(MSE.full=mean((full-obs)^2),MSE.step=mean((step-obs)^2))
#ou
pred.df %>% summarize_at(1:2,~(mean((.-obs)^2)))

## ---------------------------------------------------------
library(kernlab)
data(spam)

## ----teacher=correct--------------------------------------
set.seed(4321)
perm <- sample(nrow(spam),3000)
dapp <- spam[perm,]
dtest <- spam[-perm,]

## ----teacher=correct--------------------------------------
m.logit <- glm(type~.,data=dapp,family="binomial")

## ----echo=correct,eval=FALSE------------------------------
#  m.step <- step(m.logit,direction="backward",trace=0)

## ----echo=FALSE,eval=FALSE--------------------------------
#  #pour aller plus vite
#  save(m.step,file="m.step.logit.RData")

## ----echo=FALSE,eval=correct------------------------------
load("m.step.logit.RData")

## ----teacher=correct--------------------------------------
library(rpart)
arbre <- rpart(type~.,data=dapp)

## ----teacher=correct--------------------------------------
library(rpart.plot)
rpart.plot(arbre)
library(visNetwork)
visTree(arbre)

## ----teacher=correct--------------------------------------
prev <- data.frame(
  logit=predict(m.logit,newdata=dtest,type="response") %>% round() %>% recode_factor(`0`="nonspam",`1`="spam"),
  step=predict(m.step,newdata=dtest,type="response") %>% round() %>% recode_factor(`0`="nonspam",`1`="spam"),
  arbre=predict(arbre,newdata=dtest,type="class"))

## ----teacher=correct--------------------------------------
prev1 <- prev %>% mutate(obs=dtest$type)

## ----teacher=correct--------------------------------------
prev1 %>% summarize_at(1:3,~(mean(obs!=.))) %>% round(3)

## ----teacher=correct--------------------------------------
score <- data.frame(
  logit=predict(m.logit,newdata=dtest,type="response"),
  step=predict(m.step,newdata=dtest,type="response"),
  arbre=predict(arbre,newdata=dtest,type="prob")[,2]) %>% 
  mutate(obs=dtest$type) %>% 
  pivot_longer(-obs,names_to = "Methode",values_to="score")

## ----teacher=correct--------------------------------------
library(plotROC)
ggplot(score)+aes(d=obs,m=score,color=Methode)+geom_roc()+theme_classic()

## ----teacher=correct--------------------------------------
score %>% group_by(Methode) %>% 
  summarize(AUC=as.numeric(pROC::auc(obs,score))) %>% 
  mutate(AUC=round(AUC,3)) %>%
  arrange(desc(AUC))

## ---------------------------------------------------------
dnorm(c(-1,0,1),mean=0,sd=1)

## ---------------------------------------------------------
qnorm(c(0.05,0.5,0.95),mean=0,sd=1)

## ---------------------------------------------------------
rnorm(10,mean=0,sd=1)

## ----teacher=correct--------------------------------------
dbinom(c(1,5,10,15),size=20,prob=0.6)

## ----teacher=correct--------------------------------------
library(tidyverse)
prob <- dbinom(0:20,size=20,prob=0.6)
df <- data.frame(x=0:20,prob=prob)
ggplot(df)+aes(x=x,y=prob)+geom_bar(stat="identity",width=0.15)+theme_classic()

## ----teacher=correct--------------------------------------
X <- rbinom(5000,size=20,prob=0.6)
df1 <- data.frame(X=X)
ggplot(df1)+aes(x=X,y=..prop..)+geom_bar(width=0.15)+theme_classic()+xlim(c(0,20))

## ----teacher=correct--------------------------------------
df <- data.frame(x=seq(-3,3,by=0.01)) %>% mutate(y=dnorm(x))
ggplot(df)+aes(x=x,y=y)+geom_line()

## ----teacher=correct--------------------------------------
df1 <- data.frame(X=rnorm(5000))
ggplot(df1)+aes(x=X,y=..density..)+geom_histogram()+theme_classic()+
  geom_line(data=df,aes(x=x,y=y),color="red",size=1)

## ---------------------------------------------------------
ech1 <- runif(20)
ech2 <- runif(20)
df <- data.frame(ech1,ech2)

## ---------------------------------------------------------
df %>% summarise_all(mean)

## ---------------------------------------------------------
set.seed(1234)
df <- matrix(runif(20*5000),nrow=20) %>% as.data.frame()

## ---------------------------------------------------------
moy <- df %>% summarize_all(mean)
head(t(moy))

## ---------------------------------------------------------
moy <- data.frame(M=t(moy))
ggplot(moy)+aes(x=M,y=..density..)+geom_histogram(bins=20)+theme_classic()

## ---------------------------------------------------------
x <- seq(0.25,0.75,by=0.001)
df <- data.frame(x=x,y=dnorm(x,0.5,1/(sqrt(12*20))))
ggplot(moy)+aes(x=M,y=..density..)+geom_histogram(bins=20)+
  geom_line(data=df,aes(x=x,y=y),color="red",size=2)+theme_classic()+xlab("x")


## ----teacher=correct--------------------------------------
df1 <- matrix(runif(20*5000),nrow=20) 
df2 <- matrix(runif(50*5000),nrow=50) 
df3 <- matrix(runif(100*5000),nrow=100) 
df4 <- matrix(runif(500*5000),nrow=500)
df <- data.frame(n20=apply(df1,2,mean),n50=apply(df2,2,mean),
                  n100=apply(df3,2,mean),n500=apply(df4,2,mean))
df1 <- df %>% gather(key="taille_ech",value=x)
ggplot(df1)+aes(x=x,y=..density..)+geom_histogram(bins=50)+facet_wrap(~fct_relevel(taille_ech,"n20","n50","n100","n500"))+theme_classic()

## ----teacher=correct--------------------------------------
df1 <- matrix(rnorm(20*5000,1,2),nrow=20) 
df2 <- matrix(rnorm(50*5000,1,2),nrow=50) 
df3 <- matrix(rnorm(100*5000,1,2),nrow=100) 
df4 <- matrix(rnorm(500*5000,1,2),nrow=500)
df <- data.frame(n20=apply(df1,2,mean),n50=apply(df2,2,mean),
                  n100=apply(df3,2,mean),n500=apply(df4,2,mean))
df1 <- df %>% gather(key="taille_ech",value=x)
ggplot(df1)+aes(x=x,y=..density..)+geom_histogram(bins=50)+facet_wrap(~fct_relevel(taille_ech,"n20","n50","n100","n500"),scales="fixed")+theme_classic()

## ----teacher=correct--------------------------------------
df1 <- matrix(rbinom(20*50000,1,0.6),nrow=20) 
df2 <- matrix(rbinom(50*50000,1,0.6),nrow=50) 
df3 <- matrix(rbinom(100*50000,1,0.6),nrow=100) 
df4 <- matrix(rbinom(500*50000,1,0.6),nrow=500)
df <- data.frame(n20=apply(df1,2,mean),n50=apply(df2,2,mean),
                  n100=apply(df3,2,mean),n500=apply(df4,2,mean))
df1 <- df %>% gather(key="taille_ech",value=x)
ggplot(df1)+aes(x=x,y=..density..)+geom_histogram(bins=30)+facet_wrap(~fct_relevel(taille_ech,"n20","n50","n100","n500"),scales="fixed")+theme_classic()

## ----teacher=correct--------------------------------------
ech <- rnorm(100,1,1)

## ----teacher=correct--------------------------------------
t.test(ech)$conf.int

## ----teacher=correct--------------------------------------
mu <- 1
n <- 100
B <- 5000
X <- matrix(rnorm(n*B,mean=mu),nrow=B)

## ----teacher=correct--------------------------------------
b1 <- apply(X,1, function(x) t.test(x)$conf.int[1:2])

## ----teacher=correct--------------------------------------
b2 <- as.data.frame(t(b1))
b2 %>% mutate(test=mu>V1 & mu<V2) %>% summarize(mean(test))

## ----teacher=correct--------------------------------------
c1 <- apply(X,1, function(x) t.test(x,conf.level=0.90)$conf.int[1:2])
c2 <- as.data.frame(t(c1))
c2 %>% mutate(test=mu>V1 & mu<V2) %>% summarize(mean(test))

## ----teacher=correct--------------------------------------
t.test(iris$Sepal.Length,conf.level=0.90)$conf.int

## ----teacher=correct--------------------------------------
sep_set <- iris %>% filter(Species=="setosa") %>% select(Sepal.Length) 
t.test(sep_set,conf.level=0.90)$conf.int
#ou
iris %>% filter(Species=="setosa") %>% select(Sepal.Length) %>% t.test(conf.level=0.9)

## ----teacher=correct--------------------------------------
sep_vervin <- iris %>% filter(Species=="versicolor" | Species =="virginica") %>% select(Sepal.Length) 
t.test(sep_vervin,conf.level=0.90)$conf.int

## ---------------------------------------------------------
set.seed(12345)
res <- rbinom(100,1,0.52)

## ----teacher=correct--------------------------------------
phat <- mean(res)
phat

## ----teacher=correct--------------------------------------
n <- length(res)
binf <- phat-qnorm(0.975)*sqrt(phat*(1-phat)/n)
bsup <- phat+qnorm(0.975)*sqrt(phat*(1-phat)/n)
c(binf,bsup)

## ----teacher=correct--------------------------------------
prop.test(sum(res),n)

## ---------------------------------------------------------
library(FactoMineR)
data(decathlon)

## ----teacher=correct--------------------------------------
perf.D <- decathlon %>% filter(Competition=="Decastar") %>% select(`100m`)
t.test(perf.D)

## ----teacher=correct--------------------------------------
perf.JO <- decathlon %>% filter(Competition=="OlympicG") %>% select(`100m`)
t.test(perf.JO)

## ----teacher=correct--------------------------------------
t.test(perf.D,perf.JO)

