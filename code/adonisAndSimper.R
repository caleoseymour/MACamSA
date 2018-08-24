library("phyloseq")
library("grid")
library("RColorBrewer")
library("vegan")

#otumat=read.table('otu_p_CM.txt', head=TRUE, sep="\t", row.names=1)
#OTU.a=otu_table(otumat, taxa_are_rows=TRUE)
#temp_tax=read.table('tax_p_CM.txt', head=TRUE, sep="\t", row.names=1)
#TAX.a=tax_table(as.matrix(temp_tax))
#mydat=phyloseq(OTU.a, TAX.a)
#sdat=subset_taxa(mydat, Domain=="Bacteria")
f_physeq = prune_samples(rownames(sample_data(physeq))[sample_data(physeq)$trmnt_type_s == "CamSa"],physeq)
f_physeq = prune_taxa(taxa_sums(f_physeq) > 0, f_physeq)
PT = f_physeq
#PT=physeq
#keep = c("MAAnti1","MAAnti2","MAAnti3","MANoAnti1","MANoAnti2","MANoAnti3","MMAntiD11","MMAntiD12","MMAntiD13","MMAntiD21","MMAntiD22","MMAntiD23","MMAntiD31","MMAntiD32","MMAntiD33","MMNoAntiCK1","MMNoAntiCK2","MMNoAntiCK3","MMclindamycin1","MMclindamycin2","MMclindamycin3","MMwater1","MMwater2","MMwater3")
#ham = c("MAAnti1","MAAnti2","MAAnti3","MANoAnti1","MANoAnti2","MANoAnti3")
#mice = c("MMAntiD11","MMAntiD12","MMAntiD13","MMAntiD21","MMAntiD22","MMAntiD23","MMAntiD31","MMAntiD32","MMAntiD33","MMNoAntiCK1","MMNoAntiCK2","MMNoAntiCK3","MMclindamycin1","MMclindamycin2","MMclindamycin3","MMwater1","MMwater2","MMwater3")

##kPT = prune_samples(sample_names(PT) %in% keep, PT)
##hPT = prune_samples(sample_names(PT) %in% ham, PT)
##mPT = prune_samples(sample_names(PT) %in% mice, PT)

##keepPT = prune_taxa(taxa_sums(kPT) > 0, kPT)
##hamPT = prune_taxa(taxa_sums(hPT) > 0, hPT)
##micePT = prune_taxa(taxa_sums(mPT) > 0, mPT)

#groups = substr(sample_names(keepPT),1,nchar(sample_names(keepPT))-1)
keepPT = rarefy_even_depth(PT,rngseed=rseed)

groups = na.omit(sample_data(keepPT)$trmnt_day_s)
#groups = c(replicate(3,'Anti'),replicate(5,'Day0'),replicate(4,'CamSa'),replicate(3,'Day0'))
u.unifrac.dm = as.matrix(read.table('../exported-metrics/unweighted_unifrac_distance_matrix.tsv',
                                    sep='\t', row.names = 1, header=TRUE))
u.unifrac.dm.f = u.unifrac.dm[rownames(u.unifrac.dm) %in% sample_names(keepPT),
                            colnames(u.unifrac.dm) %in% sample_names(keepPT)]
                            
w.unifrac.dm = as.matrix(read.table('../exported-metrics/weighted_unifrac_distance_matrix.tsv',
                                    sep='\t', row.names = 1, header=TRUE))
w.unifrac.dm.f = w.unifrac.dm[rownames(w.unifrac.dm) %in% sample_names(keepPT),
                            colnames(w.unifrac.dm) %in% sample_names(keepPT)]

ugroups = unique(groups)

all.comparisons = combn(ugroups,2)
ano.list = list()
ado.list = list()
simp.list = list()
for(i in 1:ncol(all.comparisons))
{
    comp = all.comparisons[,i]
    tempPT = prune_samples(groups %in% comp, keepPT)
    temp = vegdist(t(attributes(tempPT)[[1]]), distance='bray')
    gs = data.frame(group = groups[groups %in% comp])
    ano = adonis(temp ~ group, gs)
    ado = anosim(temp, gs$group)
    simp = simper(t(attributes(tempPT)[[1]]), gs$group)
    
    ado.list[[i]] = ado
    ano.list[[i]] = ano
    simp.list[[i]] = simp
    
    names(ado.list)[i] = paste(comp,collapse='-')
    names(ano.list)[i] = paste(comp,collapse='-')
    names(simp.list)[i] = paste(comp,collapse='-')
}
otu_mat = as.matrix(t(attributes(keepPT)[[1]]))
bc_dist = vegdist(otu_mat, distance='bray')
gs = data.frame(group = groups)
ano.adonis=adonis(bc_dist ~ gs$group)
simp = simper(otu_mat, gs$group)
ano.anosim=anosim(bc_dist , gs$group)

#simp.df = cbind(data.frame(val = c(simp[[1]]$cusum[1],diff(simp[[1]]$cusum)),
#cumsum = simp[[1]]$cusum, row.names=simp[[1]]$species),tax_table(keepPT))
insert = data.frame(contrib=c(simp[[1]]$cusum[1],diff(simp[[1]]$cusum)))
rownames(insert) = rownames(summary(simp)[[1]])
rows = match(rownames(summary(simp)[[1]]),taxa_names(keepPT))
simp.df = cbind(summary(simp)[[1]],insert,tax_table(keepPT)[rows,])
sink('adonisAndanosim.txt')
cat('\n--------------------------------ADONIS RESULTS--------------------------------')
print(ano.adonis)
cat('\n--------------------------------ANOSIM RESULTS--------------------------------')
print(ano.anosim)
sink()
write.csv(simp.df,'simper.csv')
           
#pal = palette()[2:(length(unique(groups))+1)]
pal = c('#FEBF00','#006FBF')
cols = groups
for (i in 1:length(unique(groups)))
{
    cols[cols == unique(groups)[i]] = pal[i]
}

# setEPS()
# postscript("whatever.eps")
# plot(rnorm(100), main="Hey Some Data")
# dev.off()
#ord.data = ordinate(keepPT, method='NMDS', distance='bray',nsim=5000)
ord.data = ordinate(keepPT, method='NMDS', distance='bray',nsim=5000)
svglite('nmds.svg',pointsize=24,width=5.5*(9.49/7.12),height=5.5*(9.49/7.12))
par(mar = c(3.1,3.1,0.5,0.5), mgp = c(2,0.5,0), family = 'sans')
plot(ord.data,type='n',xlim=c(-0.4,0.4),ylim=c(-0.4,0.4),ps=24,xaxt='n',yaxt='n',bty='n',xlab='Axis 1', ylab = 'Axis 2')
#plot.window(, asp = 1)
labs = gsub('MA','',rownames(ord.data$points))
labs = gsub('0','0-',labs)
axis(side=1,tick=FALSE,labels=c('-0.4','-0.2','0','0.2','0.4'), at = c(-0.4,-0.2,0,0.2,0.4))
axis(side=1,tick=FALSE,labels=NULL, at = c(-0.3,-0.1,0.1,0.3))
axis(side=2,tick=FALSE,labels=c('-0.4','-0.2','0','0.2','0.4'), at = c(-0.4,-0.2,0,0.2,0.4), las = 2)
axis(side=2,tick=FALSE,labels=NULL, at = c(-0.3,-0.1,0.1,0.3), las = 2)
box(lwd=2)
NMDS = data.frame(x=ord.data$points[,1],y=ord.data$points[,2],lab=labs,group=as.character(paste0('Day ',groups)))
ord = ordiellipse(ord.data, NMDS$group, display = "sites", lty = 1, lwd = 2,
                   kind = "se", conf = 0.95, show.groups = 'Day 0')
ord = ordiellipse(ord.data, NMDS$group, display = "sites", lty = 1, lwd = 2, 
                   kind = "se", conf = 0.95, show.groups = 'Day 30')
points(x=ord.data$points[,1], y=ord.data$points[,2], bg=cols, pch=22, cex=1)
text(labels=NMDS$lab,x=NMDS$x, y=NMDS$y, cex=1, adj=c(0,0), ps = 24)
text(x=-0.4,y=0.4, labels=paste0('Stress: ', round(ord.data$stress,3)), cex=1, adj=c(0,0.11), ps = 24)
text(x=-0.4,y=0.4, labels=paste0('Adonis Pr(<F): ', round(ano.adonis$aov.tab$`Pr(>F)`[1],3),'**'), cex=1, adj=c(0,2.1), ps = 24)
dev.off()
# veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
# {
    # theta <- (0:npoints) * 2 * pi/npoints
    # Circle <- cbind(cos(theta), sin(theta))
    # t(center + scale * t(Circle %*% chol(cov)))
# }
# df_ell <- data.frame()
# for(g in levels(NMDS$group)){
    # df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                    # veganCovEllipse(cov.wt(cbind(x,y),wt=rep(1/length(x),length(x)))$cov,center=c(mean(x),mean(y)))))
                    # ,group=g))
# }
# df_ell <- data.frame()
# for(g in levels(NMDS$group)){
    # df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                  # veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                # ,group=g))
# }
# df_0 = df_ell[df_ell$group == 'day0',]
# df_30 = df_ell[df_ell$group == 'day30',]
# lines(x = df_0$x, y = df_0$y, col = "blue")
# lines(x = df_30$x, y = df_30$y, col = "red")

#labcoords = mvlabs(d.f,cols='red')

#text(labels=d.f$lab,x=d.f$x, y=d.f$y, cex=0.8, adj=c(1.1,0.5))

