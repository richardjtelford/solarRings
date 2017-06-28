fortify.nestspec <- function(nestobj){
  df <- data_frame(spec = nestobj$spec,
             freq = nestobj$freq,
             background = nestobj$background,
             periodogram_p = nestobj$periodogram$p,
         #    periodogram_f = nestobj$periodogram$f,
             powq = nestobj$powq
             )
  df %>% mutate(indsig = row_number() %in% nestobj$indsig)
}


autoplot.nestspec <- function (nestobj, loglog = "", fax = TRUE, label.sig = FALSE, 
                               speconly = FALSE, addtitle = TRUE, addlegend = TRUE, ...) 
{
  nestobj2 <- fortify.nestspec(nestobj)
  if (fax) {
    xaxval = nestobj$freq
    xlabval = "Frequency"
  }
  else {
    nestobj$freq = 1/nestobj$freq
    xlabval = "Period"
  }
  if (speconly) {
    nestobj2 <- select(nestobj2, -periodogram_p)
  }  
  nestobj3 <- gather(nestobj2, key = key, value = value, -freq, -indsig)
  
  g <- ggplot(nestobj3, aes(x = freq, y = value, colour = key, linetype = key)) +
         geom_line() +
         geom_point(data = filter(nestobj2, indsig), aes(x = freq, y = spec), show.legend = FALSE, inherit.aes = FALSE) +
         labs(x = xlabval, y = "Power", colour = "", linetype = "")
     
  if(grepl("x", loglog)){
    g <- g + scale_x_log10()
  }
  if(grepl("y", loglog)){
    g <- g + scale_y_log10()
  }
  
  if(addtitle){
    g <- g + ggtitle(sprintf("Spectrum from %s with %g%% CI from %s, W=%g", 
                     nestobj$methspec, nestobj$quantile * 100, nestobj$methtest, 
                     nestobj$nopts))
  }
         
  lims <- c( "periodogram_p", "spec", "background", "powq")
  labs <- c("periodogram", 
            sprintf("spectrum, M=%g", nestobj$M), 
            "background", 
            sprintf("%g%% CI", nestobj$quantile * 100))
  cols <- c("grey", "blue", "black", "red")
  lty <- c(rep("solid", 3), "dashed")
  if(speconly){
    lims <- lims[-1]
    labs <- labs[-1]
    cols <- cols[-1]
    lty <- lty[-1]
  }       

    g <- g + scale_colour_manual(limits = lims, labels = labs, values = cols) +
             scale_linetype_manual(limits = lims, labels = labs, values = lty)
  g
}
