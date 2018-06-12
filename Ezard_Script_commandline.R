###################################################################
## script written by Shahan Derkarabetian, implementing method of Ezard et al. (2010).
## execute in command line: $ Rscript Ezard_Script_commandline.R <data_file> <taxon_file>
###################################################################

## you must have R and these packages installed.
library("pcaPP")
library("mclust")
library("mvoutlier")
library("vegan")
library("gplots")

## command line arguments 
args <- commandArgs(trailingOnly = TRUE)
file = args[1]
sppf = args[2]

table <- read.table(args[1], sep=",")
spps <- read.table(args[2], sep=",")
rows <- nrow(table)
cols <- ncol(table)
tt <- matrix(scan(args[1], sep=",", n=as.numeric(rows)*as.numeric(cols)), as.numeric(rows), as.numeric(cols), byrow=TRUE)

## performs operations outlined in Ezard et al. (2010). this part taken directly from their publication
m1   <- PCAgrid(tt, k=as.numeric(cols), scale=mad, center=median, method="qn")
nk   <- max(which(m1$sdev^2 > bstick(m1)))
pp   <- predict(m1)[,1:nk]
m2   <- Mclust(pp)
grps <- as.numeric(m2$classification)
out01<- sign2(as.data.frame(pp), qcrit=.975)
out  <- out01$wfinal01

## reports BIC for each model: T. Ezard (pers. comm.)
eigen <-bstick(m1)
mod1  = Mclust(tt[,1:4])
bic   <-mod1$BIC[1:9,]

## prints output to screen.
nk		## nk = # of components
m2		## m2 = Mclust model choice
out		## out = outlier report [0=outlier]
bic		## bic = BIC table

## writes all output to files. makes nice plots and csv files
dir.create("Ezard_output")
write.csv(pp, "Ezard_output/retained_components.csv")
write.csv(bic, "Ezard_output/BIC.csv")
write.csv(out, "Ezard_output/outliers.csv")
write.csv(grps, "Ezard_output/groups.csv")
write.csv(eigen, "Ezard_output/squared_eigenvalues.csv")
sink("Ezard_output/model_choice.txt")
m2
sink()

sp1 <- cbind(spps, grps, out, pp)
attach(sp1)
species=sp1[,1]
col.list<-rich.colors(length(unique(species)))
palette(col.list)

write.csv(sp1, "Ezard_output/full_results.csv")

pdf("Ezard_output/taxa.pdf")
plot(pp, col=species, xlab="Comp1", ylab="Comp2", main="Taxa", cex=.75, pch=19)
dev.off()

pdf("Ezard_output/samples.pdf")
sample <- sp1[,2]
plot(pp, xlab="Comp1", ylab="Comp2", main="Sample")
text(pp, labels=sample)
dev.off()

pdf("Ezard_output/clusters.pdf")
plot(pp, pch=grps, xlab="Comp1", ylab="Comp2", main="Clusters") 
dev.off()

pdf("Ezard_output/taxa_key.pdf")
plot(0,0, type="n", col.axis="white", col.lab="white", bty="n", xaxt="n", yaxt="n")
legend(-0.25, 0.5, unique(species), col=unique(species), cex=.75, pch=19)
dev.off()

