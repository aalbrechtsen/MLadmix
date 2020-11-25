# MLadmix


```
git clone https://github.com/aalbrechtsen/MLadmix.git
cd MLadmix
g++ MLadmix.cpp -lz -lpthread  -O3 -o MLadmix
```



```


LIKE=/home/jonas/torch/EUR/eur.merged.loglike.gz
./MLadmix -l $LIKE -K 3 -Y 16 -P 50 -method 1 -headless 2
```
 

```
q<-t(read.table("eur.merged.loglike.gz.qopt"))
pop<-scan("/home/jonas/torch/EUR/eur.labels")
ord <- order(pop)
pop<-pop[ord]
q<-q[,ord]

barplot(q,space=0,col=c("darkred","darkblue","darkgreen"),border=0,xlab="admixture proportions")
    abline(v=cumsum(table(pop)))
    m <-table(pop)/2+c(0,cumsum(table(pop))[-length(table(pop))])
    mtext(names(table(pop)),1,1,adj=m/length(pop))

```
