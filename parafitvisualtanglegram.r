# This is a simple host association code with visualization of the tanglegram in the end

library(ape)
library(vegan)
library(tidyverse)
library(phytools)
library (patchwork)

links = "parafitlinks.csv"
tree_host ="host.nwk"
tree_parasite = "virus.nwk"
permutations = "permut"
N.perm = 1000


HP <- as.matrix(read.csv(links, header=TRUE))
NLinks = sum(HP)
TreeH <- read.tree(tree_host)
TreeP <- read.tree(tree_parasite)

host.D <- cophenetic(TreeH)
para.D <- cophenetic(TreeP)
host.D <- host.D[rownames(HP),rownames(HP)]
para.D <- para.D[colnames(HP),colnames(HP)]


PACo <- function (H.dist, P.dist, HP.bin){
  HP.bin <- which(HP.bin > 0, arr.in=TRUE)
  H.PCo <- pcoa(H.dist, correction="cailliez")$vectors
  P.PCo <- pcoa(P.dist, correction="cailliez")$vectors
  H.PCo <- H.PCo[HP.bin[,1],]
  P.PCo <- P.PCo[HP.bin[,2],]
  list (H.PCo = H.PCo, P.PCo = P.PCo)
}

PACo.fit <- PACo(host.D, para.D, HP)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo)

HostX <- HP.proc$X
ParY <- HP.proc$Yrot
plot(HostX, asp=1, pch=46)
points(ParY, pch=1)
arrows(ParY[,1], ParY[,2], HostX[,1], HostX[,2], length=0.12, angle=15,xpd=FALSE)
HostX <- unique(HP.proc$X)
ParY <- unique(HP.proc$Yrot)
#identify(ParY[,1], ParY[,2], rownames(ParY), offset=0.3, xpd=FALSE, cex=0.8)
#identify(HostX[,1], HostX[,2], rownames(HostX),offset=0.3, xpd=TRUE, cex= 0.8)


m2.obs <- HP.proc$ss

P.value = 0
set.seed(5)


for (n in c(1:N.perm))
{if (NLinks <= nrow(HP) | NLinks <= ncol(HP))
{flag2 <- TRUE
while (flag2 == TRUE)
{HP.perm <- t(apply(HP,1,sample))
if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE
else flag2 <- FALSE}
} else { HP.perm <- t(apply(HP,1,sample))}
  PACo.perm <- PACo(host.D, para.D, HP.perm)
  m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss
  write(m2.perm, file=permutations, sep="\t", append=TRUE)
  if (m2.perm <= m2.obs){
    cat("olala\n")
    P.value = P.value + 1
  }
}


P.value <- P.value/N.perm
cat("The observed m2 is ", m2.obs, "\np-value is", P.value, " based on ", N.perm,"permutations.")


HP.ones <- which(HP > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)
colnames (SQres.jackn) <- paste(rownames(HostX),rownames(ParY), sep="-")
t.critical = qt(0.975,NLinks-1)
for(i in c(1:NLinks))
{HP.ind <- HP
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo(host.D, para.D, HP.ind)
Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo)
res.Proc.ind <- c(residuals(Proc.ind))
res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
SQres.jackn [i, ] <- res.Proc.ind}
SQres.jackn <- SQres.jackn**2
SQres <- (residuals (HP.proc)**2)
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres))
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE)
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE)
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks)
pat.bar <- barplot(phi.mean, names.arg = " ", space = 0.25, col="white", xlab=
                     "Host-parasite link", ylab= "Squared residuals", ylim=c(0, max(phi.UCI)),
                   cex.lab=1.2)
text(pat.bar, par("usr")[3] - 0.001, srt = 330, adj = 0, labels =
       colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.6)
arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
abline(a=median(phi.mean), b=0, lty=2)


res <- parafit(host.D, para.D, HP, nperm=N.perm, test.links=TRUE, correction = "lingoes")
res
rest = as.data.frame(res$link.table)
rest


listofhost <- data.frame(Number = integer(), Name = character(), stringsAsFactors = FALSE)

for (i in 1:length(rownames(host.D))) {
  cat(i, rownames(host.D)[i], "\n")  
  listofhost <- rbind(listofhost, data.frame(Number = i, Name = rownames(host.D)[i], stringsAsFactors = FALSE))
}

listofpara <- data.frame(Number = integer(), Name = character(), stringsAsFactors = FALSE)

for (i in 1:length(rownames(para.D))) {
  cat(i, rownames(para.D)[i], "\n")  
  listofpara <- rbind(listofpara, data.frame(Number = i, Name = rownames(para.D)[i], stringsAsFactors = FALSE))
}

print(listofhost)
print(listofpara)

#perform the tree connections automatic 
## create co-phylogenetic object

TreeH2<-untangle(TreeH)
TreeP2<-untangle(TreeP)
#this is the corepa file association , in which you havae a 2 column list of associations!!!
links = read.delim("ProjectCat.txt", sep= " ", header = FALSE)
#show only the significant links!
linkssig =rest %>% filter(p.F2  < 0.05) 
linkssig

HP <- as.matrix(links)
assoc <- HP
assoc

# First, merge linkssig with listofhost to add the host names
linkssig <- merge(linkssig, listofhost, by.x = "Host", by.y = "Number", all.x = TRUE)
names(linkssig)[ncol(linkssig)] <- "HostName"  # Rename the added column to HostName

# Next, merge the result with linksofpara to add the parasite names
linkssig <- merge(linkssig, listofpara, by.x = "Parasite", by.y = "Number", all.x = TRUE)
names(linkssig)[ncol(linkssig)] <- "ParasiteName"  # Rename the added column to ParasiteName

# Now, linkssig has both HostName and ParasiteName columns
print(linkssig)

file <-file("linkssig.txt")

write_delim (linkssig, file = file, col_names = TRUE)



linkssig2 <- linkssig %>% select(ParasiteName, HostName)

hbv.cophylo2<-cophylo(TreeP2,TreeH2,
                      assoc=linkssig2)

# Plot the initial cophylo object with all associations
pdf("allassoc.pdf")
plot(hbv.cophylo, link.type = "curved", link.lwd = 1,
     link.lty = "solid", link.col = make.transparent("black", 0.80))
dev.off()
# Add the significant associations on top of the existing plot
pdf("assocsig.pdf")
plot(hbv.cophylo2, link.type = "curved", link.lwd = 1,
     link.lty = "solid", link.col = make.transparent("red", 0.90),
     add = TRUE)  # Make sure to use the add argument to overlay
dev.off()


