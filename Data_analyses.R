============================================================================
PERMANOVA -  PERMDISP - ANOVA
https://rstudio-pubs-static.s3.amazonaws.com/288625_b0cf85c5be974eb795cef778535289d8.html
============================================================================

# PERMANOVA by sediment habitat 
per<-read.csv("permdisp_habitat.csv", sep=",", header=TRUE, row.names=1)
str(per)

permanova<-adonis(per~habitat, method="bray", perm=999, by=NULL)
permanova

dis <- vegdist(per)

mod <- betadisper(dis, habitat, type=c("median"))
anova(mod)

plot(mod, label = F, hull = F, ellipse = T)
boxplot(mod)

TukeyHSD(mod)

# PERMANOVA by latitudinal gradient
per<-read.csv("permdisp_lat.csv", sep=",", header=TRUE, row.names=1)
str(per)

permanova<-adonis(per~lat, method="bray", perm=999, by=NULL)
permanova

dis <- vegdist(per)

mod <- betadisper(dis, lat, type=c("median"))
anova(mod)

plot(mod, label = F, hull = F, ellipse = T)
boxplot(mod)

TukeyHSD(mod)

#ANOVA eDNA experiment
per<-read.csv("ANOVA_eDNA.csv", sep=",", header=TRUE, row.names=1)
str(per)

permanova<-adonis(per~time, method="bray", perm=999, by=NULL)
permanova

dis <- vegdist(per)

mod <- betadisper(dis, time, type=c("median")) 
anova(mod)

plot(mod, label = F, hull = F, ellipse = T)
boxplot(mod)

TukeyHSD(mod)


  
============================================================================
nMDS 
Tutorial: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/comment-page-1/
https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
============================================================================
  

library(vegan)

mds<-read.csv("mds.csv", sep=",", header=TRUE, stringsAsFactors=FALSE) ### TO GET LEGEND - Con summary de las muestras by region 
mds<-as.matrix(mds[,2:18]) #The first column is nominal
mds<-na.omit(mds) ; 

treat=c(rep("North_Mangrove",1), rep("Center_Mangrove",8), rep("South_Mangrove",8), rep("North_Seagrass",6), rep("Center_Seagrass",14), rep("South_Seagrass",15))
symbol<-c(rep("+",17), rep("•",35)) 

png("~/MNDS.pdf", width=5*1200, height=5*1200, res=1000,  bg="transparent")

example_NMDS=metaMDS(mds, k=3)  
example_NMDS=metaMDS(mds,k=3,trymax=1000) 

ordiplot(example_NMDS,type="n", xlim = c(-1.0, 0.01), ylim = c(-1.5,3.5)))
orditorp(example_NMDS, display="sites", air=5,cex=0.7, label=F, pch=symbol)

dev.off()
  

============================================================================
Stable isotope analyses
Tutorial: https://cran.r-project.org/web/packages/simmr/vignettes/simmr.html
============================================================================

library(siar)

#Matrix 
mix = matrix(c(0.677335527	,0.254120615	, -0.05915109	,0.185506196	,2.103220114	, 0.833874026	,1.060078712	, 0.680948743	,1.826625606	,3.231257033	,-2.635782289	,-1.286943538	,0.425230717	,1.443986162	,
               2.487005521	,0.414074113	,1.466323123	,0.445839497	,1.066581123	,0.782636411	,0.545450232	,2.293174113	,
               0.42004646	,0.6943184	,2.712513687	,2.339588399	,2.783599895	,-22.88690968	,-20.78729667	,-15.31144799	,-22.98093588	,-18.4015875	,-16.68787943	,-22.34520333	,-19.83985568	, -21.83007295	,-28.10058373	,-14.06074617	,-14.40495725	,-14.93178715	,-12.20371268	,
               -14.73002879	, -11.78834974	,-22.8399449	,-6.90650433	,-11.10461174	,-25.46194251	, -4.920628591	,-33.54703808	,-2.596756046	,-12.8711635	,-14.58568559	,-13.84877625	,
               -16.20385296	), ncol=2, nrow=27)
colnames(mix) = c('δ15N','δ13C')
s_names = c("Embryophyta","Mangrove", "Seagrass", "Macroalgae")
s_means = matrix(c(4.190007778	,1.719213125	, 0.196152839	, 1.744006061	,-24.2107	,-26.5858975	,-7.777040516	,
                   -13.73921176	), ncol=2, nrow=4)
s_sds = matrix(c(2.907674494	,2.297554685	, 2.939509648	,1.588795406	,6.249038146	,
                 1.05478008	,1.581554854	,4.21474864	), ncol=2, nrow=4)
grp = as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 1,1,2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, 2, 2))


simmr_groups = simmr_load(mixtures=mix,
                          source_names=s_names,
                          source_means=s_means,
                          source_sds=s_sds,
                          group=grp)


plot(simmr_groups,group=1:2,xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     title='Isospace plot of mangrove and seagrass sediment',mix_name='Sediment')

#Running the model
simmr_groups_out = simmr_mcmc(simmr_groups)
summary(simmr_groups_out,type='diagnostics')

#Exploring the results
summary(simmr_groups_out,type='quantiles', group=c(1,2))
summary(simmr_groups_out,type=c('quantiles','statistics'),group=c(1,2))

plot(simmr_groups_out,type='boxplot',group=2,title='Seagrass')
plot(simmr_groups_out,type='boxplot',group=1,title='Mangrove')
plot(simmr_groups_out,type=c('density','matrix'),group=1,title='Mangrove')
plot(simmr_groups_out,type=c('density','matrix'),group=2,title='Seagrass')

compare_groups(simmr_groups_out,source='Embryophyta',groups=1:2)
compare_groups(simmr_groups_out,source='Macroalgae',groups=1:2)
compare_groups(simmr_groups_out,source='Seagrass',groups=1:2)
compare_groups(simmr_groups_out,source='Mangrove',groups=1:2)
