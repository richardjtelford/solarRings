#load packages
library("nest")
library("R.matlab")
library("zoo")
library("tidyverse")
library("assertthat")

##read data

#sunspots
# # #http://www.sidc.be/silso/datafiles
# Filename: SN_y_tot_V2.0.txt
# Format: plain ASCII text
# 
# Contents:
# Column 1: Gregorian calendar year (mid-year date)
# Column 2: Yearly mean total sunspot number.
# Column 3: Yearly mean standard deviation of the input sunspot numbers from individual stations.
# Column 4: Number of observations used to compute the yearly mean total sunspot number.
# Column 5: Definitive/provisional marker. A blank indicates that the value is definitive. A '*' symbol indicates that the yearly average still contains provisional daily values and is subject to a possible revision.

ss <- read.table("SN_y_tot_V2.0.txt", header = FALSE)
names(ss) <- c("year", "sunspots", "sunspotSD", "nobs")
ss <- ss %>% filter(between(year, 1750, 1951))


#ring widths
rw <- readMat("NorthernHemisphereRW.mat")

#reformat rw
rwl<-lapply(1:length(rw[[1]][1,1,1]$meta.data[1,1,1]$site.name), function(i){
       out<-list()
       out$site.name <- rw[[1]][1, 1, 1]$meta.data[1, 1, 1]$site.name[[i]][[1]][1, 1]
       out$last.name <- rw[[1]][1, 1, 1]$meta.data[2, 1, 1]$last.name[[i]][[1]][1, 1]
       out$species.code <- rw[[1]][1, 1, 1]$meta.data[3, 1, 1]$species.code[[i]][[1]][1, 1]
       out$lat <- rw[[1]][1, 1, 1]$meta.data[4, 1, 1]$lat[i, 1]
       out$long <- rw[[1]][1, 1, 1]$meta.data[5, 1, 1]$lon[i, 1]
       out$elev <- rw[[1]][1, 1, 1]$meta.data[6, 1, 1]$elev[i, 1]
       tim <- rw[[1]][2, 1, 1]$time[, 1]
       crns <- rw[[1]][3, 1, 1]$crns[, i]
       out$crns <- crns[tim > 1749 & tim < 1951]
       out
})
rm(rw)#clean up

meta <- purrr::map_df(rwl, function(x){x[1:6]})

#check data
assert_that(all(sapply(rwl, function(x)sum(!is.finite(x$crns))) == 0))

#
#map
mp <- map_data("world", ylim = c(-10, 90))
map1 <- ggplot(meta, aes(x = long, y = lat, colour = I("red"))) +
  geom_map(data = mp, map = mp, aes(map_id = region), colour  = "grey50", fill = "grey50") +
  geom_point() +
  coord_quickmap(ylim = c(0, 82)) +
  labs(x = "Longitude E°", y = "Latitude N°")
map1

#frequency domain
ssf <- nestspec.sig(data=zoo(select(ss, sunspots)), methtest = "medsmooth", pval = 0.9, methspec = "multitaper", detr = 0, W = 15)
autoplot.nestspec(ssf, speconly = TRUE)
 

z <- zoo(rwl[[1]]$crns, order.by=1750:1950)
n <- nestspec.sig(data=z, methtest = "medsmooth", pval = 0.9, methspec = "multitaper", detr = 0, W = 15)
n$spec > n$powq
  
autoplot.nestspec(n, speconly = TRUE)

sig.9 <- sapply(rwl, function(x) {
  z <- zoo(x$crns, order.by = 1750:1950)
  n <- nestspec.sig(
      data = z,
      methtest = "medsmooth",
      pval = 0.9,
      methspec = "multitaper",
      detr = 0,
      W = 15
    )
  n$spec > n$powq
})

#remove first few frequencies - poor estimate of median
trim <- 4


summary.plot <- data_frame(frequency = n$freq, percent = rowMeans(sig.9) * 100) %>% 
  slice(-(1:trim)) %>% 
  ggplot(aes(x = frequency, y = percent)) + 
  geom_hline(yintercept = 10, colour = "grey", linetype = "dashed") +
  geom_path() +
  labs(x = expression(Frequency~years^-1), y = "Percent")
summary.plot  

#plot
#upper - spectrum ss
#lower - proportion significant
xlim <- xlim(0, 0.5)
no_x_axis <- theme( axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

cowplot::plot_grid(
  autoplot.nestspec(ssf, speconly = TRUE, addtitle = FALSE) + 
    xlim + no_x_axis +
    theme(legend.position = c(0.85, 0.6), legend.title = element_blank()),
  autoplot.nestspec(n, speconly = TRUE, addtitle = FALSE) + 
    xlim +  no_x_axis +
    theme(legend.position = "none"),
  summary.plot + xlim, 
  ncol = 1, 
  align = "v", 
  rel_heights = c(1, 1, 1.4)
)



# x11(height=9, point=16*3/2)
# par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(2.7,2.7,1,2.7), mgp=c(1.2,.1,0), tcl=0.2)
# with(ssf,matplot(freq[-(1:trim)],cbind(background, powq, spec)[-(1:trim),], col=c(1,2,4), type="l", lty=c(2,3,1), log="", bty="n", xaxt="n", lwd=2))
# legend("topright", legend=c("Spectral power", "Median smoother", "90% significance level"), lty=c(1,2,3), col=c(4,1,2), cex=.8, bty="n", lwd=2)
# title(main="a)", line=-1, adj=0.05)
# with(n,matplot(freq[-(1:trim)],cbind(background, powq, spec)[-(1:trim),], col=c(1,2,4), type="l", lty=c(2,3,1), log="", bty="n", xaxt="n", yaxt="n",lwd=2))
# title(main="b)", line=-1, adj=0.05)
# axis(4, at=seq(0,.3,.1))
# plot(n$freq[-(1:trim)], rowMeans(sig.9)[-(1:trim)]*100, type="l", bty="u", ylim=c(0,21),lwd=2, yaxt="n", col=4)
# axis(2, at=seq(0,30,10))
# title(main="c)", line=-1, adj=0.05)
# abline(h=10, col=1, lty=3,lwd=2)
# box("inner")
# title(xlab=expression("Frequency years"^"-1"), outer=TRUE)
# title(ylab="Percent records significant", adj=0.0, outer=TRUE)
# title(ylab="Spectral power", adj=0.9, outer=TRUE)
# mtext(side=4, text="Spectral power", outer=TRUE, line=1.2, cex=.66)

#map significant
mx <- which.max(ssf$spec)
1/ssf$freq[mx]
1/ssf$freq

sigr<-apply(sig.9[(mx-2):(mx+2),], 2, any)

map1 + aes(colour = sigr) +
  scale_colour_discrete(limits = c("TRUE", "FALSE"))
  labs(colour = "Significant\nperiodicity")


write.table(cbind(meta, sig=sigr), file="sig.trees.txt", quote=FALSE, sep="\t", row.names=FALSE)

#sensitivity - signal to noise ratio
#standardise ss and treering data
#add/subtract ss * k to tr
# k in 1/c(1, 10, 100) etc
#spectral analysis
#proportion records significant

snf <- c(1, 1.5, 2, 3, 5, 10, 50, 100)
signalNoise <- c(0, 1/c(-snf, snf)) 
scaled_ss <- ss %>% select(sunspots) %>% mutate(sunspots = scale(sunspots))

signal_noise_analysis <- sapply(signalNoise, function(k){
  scaled_ss <- scaled_ss * k
  sig.9 <- sapply(rwl, function(x) {
    z <- zoo(scale(x$crns) + scaled_ss, order.by = 1750:1950)
    n <- nestspec.sig(
      data = z,
      methtest = "medsmooth",
      pval = 0.9,
      methspec = "multitaper",
      detr = 0,
      W = 15
    )
    sig <- n$spec > n$powq
    sig_11 <- sig[mx] 
    sig_11
  })
  mean(sig.9)
})


sensivity_plot <- data_frame(sign = factor(sign(signalNoise)), 
                             signalNoise = abs(signalNoise), 
                             power = signal_noise_analysis) %>% 
  ggplot(aes(x = signalNoise, y = power, colour = sign, group = sign)) + 
    geom_line() +
    geom_point() 
sensivity_plot

sensivity_plot + aes(weight = length(rwl)) +
     geom_smooth(method = "gam", formula = y ~ x + I(x^2), method.args = list(family = binomial))



plot(0:750, dbinom(0:750, size = 750, prob = 0.9))


cor.mat <- sapply(rwl, function(x)x$crns) %>% cor() %>% as.dist() %>% hist()


#time domain
#correlations with lags

cor.9 <- sapply(rwl, function(x) {
  cor(x$crns, ss)
})
