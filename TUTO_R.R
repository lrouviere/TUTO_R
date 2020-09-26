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
c(T,F,T)

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

## ---------------------------------------------------------
data(iris)
head(iris)

## ---------------------------------------------------------
library(lattice)
data("ethanol")

## ---------------------------------------------------------
data("presidents")
df <- matrix(presidents,ncol=4,byrow=T)

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

## ---------------------------------------------------------
df1 <- tibble(name=c("Mary","Peter","John","July"),age=c(18,25,21,43))
df2 <- tibble(name=c("Zac","Julian"),age=c(23,48))
df3 <- tibble(size=c(154,178,182,134,142),name1=c("Peter","Mary","July","John","stef"))
df1
df2
df3

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

## ---------------------------------------------------------
library(rAmCharts)
amHist(ozone$maxO3)
amPlot(ozone,col=c("T9","T12"))
amBoxplot(maxO3~vent,data=ozone)

## ----echo=correct,eval=FALSE------------------------------
#  title("Représentation de la fonction sinus")

## ----echo=FALSE,eval=correct------------------------------
#  plot(x,sin(x),type="l")
#  title("Représentation de la fonction sinus")

## ----echo=correct,eval=FALSE------------------------------
#  abline(v=0,lty=2)

## ----echo=FALSE,eval=correct------------------------------
#  plot(x,dnorm(x),type="l")
#  abline(v=0,lty=2)

## ----echo=correct,eval=FALSE------------------------------
#  lines(x,dt(x,5),col=2)
#  lines(x,dt(x,30),col=3)

## ----echo=FALSE,eval=correct------------------------------
#  plot(x,dnorm(x),type="l")
#  abline(v=0,lty=2)
#  lines(x,dt(x,5),col=2)
#  lines(x,dt(x,30),col=3)

## ----echo=correct,eval=FALSE------------------------------
#  legend("topleft",legend=c("Normal","Student(5)","Student(30)"),
#     col=1:3,lty=1)

## ----echo=FALSE,eval=correct------------------------------
#  plot(x,dnorm(x),type="l")
#  abline(v=0,lty=2)
#  lines(x,dt(x,5),col=2)
#  lines(x,dt(x,30),col=3)
#  legend("topleft",legend=c("Normal","Student(5)","Student(30)"),
#     col=1:3,lty=1)

## ----eval=correct,echo=TRUE-------------------------------
#  couleurs <- c("yellow", "magenta", "orange", "cyan",
#            "grey", "red", "green", "blue")

## ---------------------------------------------------------
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

## ----eval=FALSE-------------------------------------------
#  ggplot(diamonds2)+aes(x=carat,y=price)

## ----eval=FALSE-------------------------------------------
#  ggplot(diamonds2)+aes(x=carat,y=price,color=cut)

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=carat,y=price,color=cut)+geom_point()

## ---------------------------------------------------------
D <- data.frame(X=seq(-2*pi,2*pi,by=0.01))
ggplot(D)+aes(x=X,y=sin(X))+geom_line()

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=price)+geom_histogram(bins=40)

## ---------------------------------------------------------
ggplot(diamonds2)+aes(x=price,y=..density..)+geom_histogram(bins=40)

## ----eval=FALSE,echo=TRUE---------------------------------
#  ggplot(diamonds2)+aes(x=price,y=..density..)+stat_bin()

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

## ---------------------------------------------------------
data(mtcars)
summary(mtcars)

## ----echo=FALSE,eval=TRUE---------------------------------
ggplot(data=diamonds) + geom_boxplot(aes(x=cut,y=carat,fill=cut)) 
ggplot(data=diamonds) + geom_boxplot(aes(x=cut,y=carat,fill=cut))+coord_flip()
ggplot(data=diamonds) + geom_density(aes(x=carat,y=..density..)) +  facet_grid(cut~.)

## ----echo=correct,eval=FALSE------------------------------
#  ggplot(data=diamonds) + geom_boxplot(aes(x=cut,y=carat,fill=cut))
#  ggplot(data=diamonds) + geom_boxplot(aes(x=cut,y=carat,fill=cut))+coord_flip()
#  ggplot(data=diamonds) + geom_density(aes(x=carat,y=..density..)) +  facet_grid(cut~.)

## ----echo=correct,eval=TRUE-------------------------------
Q1 <- diamonds %>% group_by(cut) %>% 
  summarize(q1=quantile(carat,c(0.25)),q2=quantile(carat,c(0.5)),
        q3=quantile(carat,c(0.75)))
quantildf <- Q1%>% gather(key="alpha",value="quantiles",-cut)
ggplot(data=diamonds) + geom_density(aes(x=carat,y=..density..))+
  facet_grid(cut~.) +
  geom_vline(data=quantildf,aes(xintercept=quantiles),col=alpha("black",1/2))

## ----echo=FALSE,eval=TRUE---------------------------------
library(ggstance)
ggplot(data=diamonds) +
  geom_boxploth(data=diamonds,aes(y=-0.5,x=carat,fill=cut)) +
  geom_density(aes(x=carat,y=..density..)) +  facet_grid(cut~.) +
  geom_vline(data=quantildf,aes(xintercept=quantiles),col=alpha("black",1/2))+
  ylab("")

## ----echo=correct,eval=FALSE------------------------------
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

## ---------------------------------------------------------
ColorPal1 <- colorNumeric(scales::seq_gradient_pal(low = "#132B43", high = "#56B1F7",
                                               space = "Lab"), domain = c(0,1))
ColorPal2 <- colorNumeric(scales::seq_gradient_pal(low = "red", high = "black", 
                                               space = "Lab"), domain = c(0,1))

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

## ---------------------------------------------------------
m <- 500
p <- 5
set.seed(12345)
X.mat <- matrix(rnorm(m*p),ncol=5)
eps <- rnorm(m,mean = 0,sd=0.5)
df.test <- data.frame(X.mat,eps)
df.test <- df.test %>% mutate(Y=X1+X2+X3+X4+X5+eps) %>% select(-eps)

## ---------------------------------------------------------
n <- 1000
p <- 105
set.seed(1234)
X.mat <- matrix(rnorm(n*p),ncol=p)
eps <- rnorm(n,mean = 0,sd=0.5)
df <- data.frame(X.mat,eps)
df <- df %>% mutate(Y=X1+X2+X3+X4+X5+eps) %>% select(-eps)

## ---------------------------------------------------------
m <- 300
p <- 105
set.seed(12345)
X.mat <- matrix(rnorm(m*p),ncol=p)
eps <- rnorm(m,mean = 0,sd=0.5)
df.test <- data.frame(X.mat,eps)
df.test <- df.test %>% mutate(Y=X1+X2+X3+X4+X5+eps) %>% select(-eps)

## ---------------------------------------------------------
library(kernlab)
data(spam)

## ----echo=correct,eval=FALSE------------------------------
#  m.step <- step(m.logit,direction="backward",trace=0)

## ----echo=FALSE,eval=FALSE--------------------------------
#  #pour aller plus vite
#  save(m.step,file="m.step.logit.RData")

## ----echo=FALSE,eval=correct------------------------------
#  load("m.step.logit.RData")

## ---------------------------------------------------------
dnorm(c(-1,0,1),mean=0,sd=1)

## ---------------------------------------------------------
qnorm(c(0.05,0.5,0.95),mean=0,sd=1)

## ---------------------------------------------------------
rnorm(10,mean=0,sd=1)

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


## ---------------------------------------------------------
set.seed(12345)
res <- rbinom(100,1,0.52)

## ---------------------------------------------------------
library(FactoMineR)
data(decathlon)

