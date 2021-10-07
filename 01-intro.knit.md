# Introduction {#intro}

**R** est un logiciel libre et gratuit principalement dédié aux analyses statistiques et aux représentations graphiques. Il est gratuit et librement distribué par le **CRAN** (Comprehensive R Archive Network) à l'adresse suivante : [https://www.r-project.org](https://www.r-project.org). 

L'installation varie d'un système d'exploitation à l'autre (Windows, Mac OS, Linux) mais elle est relativement simple, il suffit de suivre les instructions.


**RStudio** est une interface facilitant l'utilisation de **R**. Elle est également gratuite et librement distribuée à l'adresse [https://www.rstudio.com](https://www.rstudio.com). 

L'interface **RStudio** est divisée en 4 fenêtres :

* *Console* où on peut entrer et exécuter des commandes (taper 1+2)
* *Environnement, History* où on peut visualiser les objets construits (taper a <- 1+2 dans la console)
* *Files, Plots...* où on peut visualiser les répertoires et fichiers de l'espace de travail, les graphiques, intaller des packages...
* *R script* où on conserve les lignes de commandes ainsi que les commentaires sur le travail effectué. Il faut penser à sauvegarder régulièrement ce fichier.


## R Script

Il existe différentes façons de travailler sur RStudio. De façon classique, on peut

* ouvrir un **script**.
* entrer les commandes dans le script.
* regarder les sorties dans la console (en cliquant sur le bouton **run**).
* sauver le script.

## Packages

Un package est une ensemble de programmes et fonctions **R** qui complètent les fonctions existantes par défaut dans le logiciel. Un package est généralement dédié à une méthode ou un champ d'application spécifique. Il existe plus de 18 000 packages disponibles sur le **CRAN**  [https://cran.r-project.org](https://cran.r-project.org). On installe un package en

* utilisant le fonction **install.packages** dans la console.
ou
* ou cliquant sur le bouton **Packages**.

Une fois le package installé sur la machine, on l'installe avec la fonction **library** :


```r
install.packages(package.name)
library(packages.name)
```


\BeginKnitrBlock{exercise}\iffalse{-91-73-110-115-116-97-108-108-97-116-105-111-110-32-101-116-32-99-104-97-114-103-101-109-101-110-116-93-}\fi{}<div class="exercise"><span class="exercise" id="exr:exo-intro-package"><strong>(\#exr:exo-intro-package)  \iffalse (Installation et chargement) \fi{} </strong></span></div>\EndKnitrBlock{exercise}

1. Exécuter

    
    ```r
    iris %>% summarize(mean_Petal=mean(Petal.Length))
    ```
    Que se passe t-il ?














