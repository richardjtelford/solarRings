## ---- loadPackages
#load packages
library("nest")
library("R.matlab")
library("zoo")
library("tidyverse")
library("assertthat")
library("ggfortify")
library("nlme")
library("lmtest")
source("autoplot.nestpec.r")

## ---- importData
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

ss <- read.table("data/SN_y_tot_V2.0.txt", header = FALSE)
names(ss) <- c("year", "sunspots", "sunspotSD", "nobs")
ggplot(ss, aes(x = year, y = sunspots)) + geom_line()


#ring widths
rw <- readMat("data/NorthernHemisphereRW.mat")

#find period with all series
series_depth <- data_frame(time = rw[[1]][2, 1, 1]$time[, 1], depth = rowSums(is.finite(rw[[1]][3, 1, 1]$crns))) 


ggplot(series_depth, aes(x = time, y = depth)) + geom_line() + xlim(1700, 2010)

time_limits <- series_depth %>% 
  filter(depth == ncol(rw[[1]][3, 1, 1]$crns)) %>% 
  summarise(max = max(time), min = min(time))

#extract solar data with all ring with series
ss <- ss %>% filter(between(year, time_limits$min, time_limits$max + 1))

#reformat rw
rwl <- lapply(1:length(rw[[1]][1,1,1]$meta.data[1,1,1]$site.name), function(i){
       out<-list()
       out$site.name <- rw[[1]][1, 1, 1]$meta.data[1, 1, 1]$site.name[[i]][[1]][1, 1]
       out$last.name <- rw[[1]][1, 1, 1]$meta.data[2, 1, 1]$last.name[[i]][[1]][1, 1]
       out$species.code <- rw[[1]][1, 1, 1]$meta.data[3, 1, 1]$species.code[[i]][[1]][1, 1]
       out$lat <- rw[[1]][1, 1, 1]$meta.data[4, 1, 1]$lat[i, 1]
       out$long <- rw[[1]][1, 1, 1]$meta.data[5, 1, 1]$lon[i, 1]
       out$elev <- rw[[1]][1, 1, 1]$meta.data[6, 1, 1]$elev[i, 1]
       tim <- rw[[1]][2, 1, 1]$time[, 1]#oldest first
       crns <- rw[[1]][3, 1, 1]$crns[, i]
       out$crns <- crns[tim >= time_limits$min & tim <= time_limits$max]
       out
})
rm(rw)#clean up

meta <- purrr::map_df(rwl, function(x){x[1:6]})

#check data
assert_that(all(sapply(rwl, function(x)sum(!is.finite(x$crns))) == 0))

## ---- map1
#map
mp <- map_data("world", ylim = c(-10, 90))
map1 <- ggplot(meta, aes(x = long, y = lat, colour = I("red"))) +
  geom_map(data = mp, map = mp, aes(map_id = region), colour  = "grey50", fill = "grey50") +
  geom_point() +
  coord_quickmap(ylim = c(0, 82)) +
  labs(x = "Longitude E°", y = "Latitude N°")
map1

## ---- frequency
#frequency domain
ssf <- nestspec.sig(data=zoo(select(ss, sunspots)), methtest = "medsmooth", pval = 0.9, methspec = "multitaper", detr = 0, W = 15)
#autoplot.nestspec(ssf, speconly = TRUE)
 
example <- 1
z <- zoo(rwl[[example]]$crns)
n <- nestspec.sig(data=z, methtest = "medsmooth", pval = 0.9, methspec = "multitaper", detr = 0, W = 15)
n$spec > n$powq
  
#autoplot.nestspec(n, speconly = TRUE)

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


## ---- fig2
#remove first few frequencies - poor estimate of median
trim <- 4


summary.plot <- data_frame(frequency = n$freq, percent = rowMeans(sig.9)) %>% 
  slice(-(1:trim)) %>% 
  ggplot(aes(x = frequency, y = percent)) + 
  geom_hline(yintercept = 0.1, colour = "grey", linetype = "dashed") +
  geom_path() +
  labs(x = expression(Frequency~years^-1), y = "Periodicity (p < 0.1)") +
  scale_y_continuous(labels = scales::percent)
summary.plot  

#upper - spectrum ss
#lower - proportion significant
xlim <- scale_x_continuous(limits = c(0, 0.5), expand = c(0.02, 0))
no_x_axis <- theme( axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig2 <- cowplot::plot_grid(
  autoplot.nestspec(ssf, speconly = TRUE, addtitle = FALSE) + 
    xlim + no_x_axis +
    theme(legend.position = c(0.85, 0.6), legend.title = element_blank()),
  autoplot.nestspec(n, speconly = TRUE, addtitle = FALSE) + 
    xlim +  no_x_axis +
    theme(legend.position = "none"),
  summary.plot + xlim, 
  ncol = 1, 
  align = "v", 
  rel_heights = c(1, 1, 1.4), 
  labels = "auto"
)
fig2

## ---- mapSignificant
#mapsignificant
mx <- which.max(ssf$spec)
1/ssf$freq[mx]
1/ssf$freq

sigr <- apply(sig.9[(mx-2):(mx+2),], 2, any)

map1 + aes(colour = sigr) +
  scale_colour_discrete(limits = c("TRUE", "FALSE"))
  labs(colour = "Significant\nperiodicity")

#significant species  
meta %>% 
  group_by(species.code) %>%  
  mutate(n = n(), species.code2 = if_else(n < 5, "rare", species.code)) %>% 
  ungroup() %>% 
  mutate( sig = as.character(sigr))  %>% 
  count(species.code2, sig) %>% 
  spread(key = sig, value = nn) %>% 
  print(n = 100)
  

map1 %+% (meta %>% mutate(sig = sigr) %>% group_by(species.code) %>% filter(n() > 5)) + 
  aes(colour = sig) +
  scale_colour_discrete(limits = c("TRUE", "FALSE")) +
  labs(colour = "Significant\nperiodicity") +
  facet_wrap(~species.code)

## ---- spectralSensitivity

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

## ---- timeDomain
#correlations #with lags 

cor.9 <- map_df(rwl, function(x) {
  correlation <- cor.test(x$crns, ss$sunspots)[c("estimate", "p.value")]
  correlation
})
cor.9 <- cor.9 %>% mutate(sig = p.value < 0.1)
mean(cor.9$sig)

ggplot(cor.9, aes(x = estimate, fill = sig)) + geom_histogram()
nbins <- 30
ggplot(cor.9, aes(x = p.value)) + geom_histogram(breaks = seq(0, 1, length = nbins)) + geom_hline(yintercept = length(rwl)/nbins)

#regression - makes more sense than correlation as we know the direction of causality
reg_res <- map_df(rwl, function(x){
  dat <- data_frame(rw = x$crns, s = ss$sunspots)
  lm_mod <- lm(rw ~ s, data = dat)
  dw_res <- dwtest(lm_mod)
  list(coef = coef(lm_mod)[2], lm_pvalue = anova(lm_mod)$`Pr(>F)`[1], dw_stat = dw_res$statistic, dw_pvalue = dw_res$p.value)
})

ggplot(reg_res, aes(x = coef, fill = lm_pvalue < 0.05)) + 
  geom_histogram()


ggplot(reg_res, aes(x = lm_pvalue)) + 
  geom_histogram(breaks = seq(0, 1, length = nbins)) + 
  geom_hline(yintercept = length(rwl)/nbins)



#but data are autocorrelted
ggplot(reg_res, aes(x = dw_stat, fill = dw_pvalue < 0.05)) + 
  geom_histogram()

#acf/pacf for ring
plot_p_acf <- function(i, s = ss$sunspots){
  mod <- lm(rwl[[i]]$crns ~ s)
  resids <- resid(mod)
  ACF <- autoplot(acf(resids, plot = FALSE))
  PACF <- autoplot(pacf(resids, plot = FALSE))
  gridExtra::grid.arrange(ACF, PACF)
}
plot_p_acf(1)
plot_p_acf(10)


#know direction of supposed causality can use regression

gls_res <- map_df(rwl, function(x) {print(".")
  dat <- data_frame(rw = x$crns, s = ss$sunspots)
  mods <- list()
  mods$mod1 <- gls(rw ~ s, dat = dat)
  mods$mod2 <- gls(rw ~ s, dat = dat, correlation = corARMA(p = 1))
  mods$mod3 <- gls(rw ~ s, dat = dat, correlation = corARMA(q = 1))
  try(mods$mod4 <- gls(rw ~ s, dat = dat, correlation = corARMA(p = 1, q = 1)))
  aics <- sapply(mods, AIC)
  lowest <- which.min(aics)
  if(length(aics) == 3) aics <- c(aics, NA)
  c(list(
    coef = coef(mods[[lowest]])[2], 
    pvalue = anova(mods[[lowest]])$`p-value`[2],
    lowest = lowest
    ), as.list(aics - min(aics, na.rm = TRUE)))
})

gls_res %>% count(lowest)

gls_res <- gls_res %>% mutate(sig = pvalue < 0.1)
mean(gls_res$sig)
gls_res %>% ggplot(aes(x = pvalue)) + 
  geom_histogram(breaks = seq(0, 1, length = nbins)) + 
  geom_hline(yintercept = length(rwl)/nbins)

gls_res %>% ggplot(aes(x = mod2, y = mod4)) + geom_point()






