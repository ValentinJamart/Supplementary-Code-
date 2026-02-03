#Title: Astronomical calibration of the middle Cambrian in Baltica: global carbon cycle synchronization and climate dynamics
#Authors: Valentin JAMART, Damien PAS, Linda A. HINNOV, Jorge E. SPANGENBERG, Thierry ADATTE, Arne T. NIELSEN, Niels H. SCHOVSBO, Nicolas THIBAULT, Michiel ARTS, Allison C. DALEY
#Last update: July 23, 2025 - Valentin Jamart

--------------------------------------------------------
   
#Work environment, packages and data loading ----


#Setting of the work environment
#setwd("C:/Users/.....")

#Loading of the packages that will be used in this study
library(WaverideR)
library(astrochron)
library(matrixStats)

#Loading of the dataset file
Alb <- read.csv("Z score.csv", sep=";")

--------------------------------------------------------
  
#Step 1: Data selection, resampling and detrending ----

#* 1.1. Selecting Titanium data ----

# Selection of the data
alb_Ti <- cbind(Alb$DepthAdj,Alb$Ti)
alb_Ti <- na.omit(alb_Ti)
alb_Ti[!is.finite(alb_Ti)] <- NA
alb_Ti <- na.omit(alb_Ti)

#Isolation of the studied time interval
alb_Ti <- iso(dat=alb_Ti, xmin=73.452, xmax=100)
alb_Ti <- linterp(alb_Ti, genplot=T)

#Saving of the 1 mm (non resampled) Ti series
alb_Ti_1mm<-linterp(alb_Ti, genplot=T)

#Resampling of the Ti series every 5 mm
alb_Ti <-  linterp(alb_Ti, dt=0.005, genplot=T)

#Saving of the non detrended Ti series resampled every 5 mm
alb_Ti_ndet <-  linterp(alb_Ti, genplot=T)
graphics.off()

#Plotting of the non detrended 5mm Ti series
plot(alb_Ti_ndet, type = "l", 
     xlab = "Depth Adjusted (m)", 
     ylab = "Ti Z score non-detrended"
     )

#The Ti data have been selected, from the composite core, 
#extending from 73.452 m (= anchoring depth with Zhao et al (2022b) dataset) 
#to 93.609 m (= end of the Alum Shale Formation). Then the data are interpolated 
#to ensure continuity of the Ti record every 5 mm. These data are not detrended 
#and will be used in the Continuous Wavelet analysis (CWT).




#* 1.2. Detrending of the Ti data----

#Detrending (20% LOWESS) of the Ti series resampled every 5 mm
alb_Ti <- noLow(alb_Ti, smooth = 0.2, output = 1, genplot = T) 
alb_Ti <- linterp(alb_Ti, genplot=T)
graphics.off()

plot(alb_Ti, type = "l", 
     xlab = "Depth Adjusted (m)", 
     ylab = "Ti Z score detrended"
     )

#The 5mm spaced Ti series is detrended using a 20% LOWESS regression. 
#These data will be used for step 2 to 10 of the protocol. The step 1 of the 
#protocol was also performed for the other detrital elements (Al, Si, K and Zr)

--------------------------------------------------------
  
#Step 2: MTM in depth domain ----

#* 2.1. MTM on the non-detrened dataset ----

alb_Ti_MTM_ndet<-mtm(alb_Ti_ndet, xmax=15, siglevel = 0.9,ar1 = T, tbw= 2, output = 1)
alb_Ti_MTM_ndet <- iso(dat=alb_Ti_MTM_ndet,xmin=0,xmax=15, genplot=F)
graphics.off()

plot(alb_Ti_MTM_ndet$Frequency, 
     alb_Ti_MTM_ndet$Power, 
     type = 'l', 
     xlab = "Frequency (cycles/m)", 
     ylab = "Variance"
     )

lines(alb_Ti_MTM_ndet$Frequency, 
      alb_Ti_MTM_ndet$AR1_fit, 
      type = 'l', 
      col = 'black', 
      lwd = 2
      )

lines(alb_Ti_MTM_ndet$Frequency, 
      alb_Ti_MTM_ndet$AR1_90_power, 
      type = 'l', 
      col = 'red', 
      lwd = 2
      )

lines(alb_Ti_MTM_ndet$Frequency, 
      alb_Ti_MTM_ndet$AR1_95_power, 
      type = 'l', 
      col = 'blue', 
      lwd = 2
      )

lines(alb_Ti_MTM_ndet$Frequency, 
      alb_Ti_MTM_ndet$AR1_99_power, 
      type = 'l', 
      col = 'green', 
      lwd = 2
      )

#The 2 pi-MTM on the non-detrended data is performed to show the entire 
#range of cyclicities in the signal. This MTM is not used afterward.




#* 2.2. MTM on the detrended dataset ----

alb_Ti_MTM <- mtm(alb_Ti, xmax=15, siglevel = 0.9, ar1 = T, tbw= 2, output = 1)
alb_Ti_MTM <- iso(dat=alb_Ti_MTM, xmin=0, xmax=15, genplot=F)
graphics.off()

plot(alb_Ti_MTM$Frequency, 
     alb_Ti_MTM$Power, 
     type = 'l', 
     xlab = "Frequency (cycles/m)", 
     ylab = "Variance",
     )

lines(alb_Ti_MTM$Frequency, 
      alb_Ti_MTM$AR1_fit, 
      type = 'l', 
      col = 'black', 
      lwd = 2
      )

lines(alb_Ti_MTM$Frequency, 
      alb_Ti_MTM$AR1_90_power, 
      type = 'l', 
      col = 'red', 
      lwd = 2
      )

lines(alb_Ti_MTM$Frequency, 
      alb_Ti_MTM$AR1_95_power, 
      type = 'l', 
      col = 'blue', 
      lwd = 2
      )

lines(alb_Ti_MTM$Frequency, 
      alb_Ti_MTM$AR1_99_power, 
      type = 'l', 
      col = 'green', 
      lwd = 2
      )

#The 2 pi-MTM on the detrended data is performed to show the range of 
#cyclicities in the detrended signal. Allowing to visually identify > 90 % CL 
#periodicities.




#* 2.3. Identification and extraction of the >90% CL periodicities ----

alb_Ti_MTM<-mtm(alb_Ti, xmax=15, siglevel = 0.9, ar1 = T, tbw= 2, output = 2)
alb_Ti_MTM <- iso(dat=alb_Ti_MTM, xmin=0, xmax=15, genplot=F)
graphics.off()

write.csv(alb_Ti_MTM,"alb_Ti_sigfreq_MTM_rsp0.005_det0.20_Zscore_ACL.csv")

--------------------------------------------------------
  
#Step 3: Ratio of frequencies ----

#* 3.1. Ratio of the significant frequencies extracted from the MTM (in Excel) ----

#The comparison between the ratios of the observed >90% CL and the theoretical 
#ratio of the Milankovitch cycles is conducted in a separate Excel file. 
#See Supplementary Materials 6 for more details. Based on Zhao et al (2022b) paper, 
#the 405 kyr cycle should be recorded between 1.75 and 2.25 m but nothing is 
#visible in the Ti data (but is is observable within the other detrital elements). 
#However, all the other frequency bands observed by Zhao et al (2022b), and 
#corresponding to Milankovitch periodicities are identify in the Ti data.
#In addition, a significant peak is located around 0.85 m where the 173-kyr long 
#obliquity (Inclination metronome) should be identified and this explains why we 
#focus our study on the 173 kyr metronome. The duration of target cycles are 
#recalculated from the literature (Laskar et al (2004); Laskar (2020); 
#Wu et al (2024)).

--------------------------------------------------------
  
#Step 4: EHA and CWT ----

#* 4.1. Evolutive Harmonic Analysis (EHA) ----

eha(alb_Ti, 
    tbw=2, 
    fmin=0.01, 
    fmax=15, 
    step=0.005, 
    win=6, 
    demean=T, 
    detrend=T, 
    siglevel=0.90, 
    sigID=F, 
    ydir=1, 
    output=0, 
    pl=1, 
    palette=1, 
    centerZero=T, 
    ncolors=100, 
    genplot=4, 
    verbose=T
    )

#The EHA is performed on the Ti detrended series and confirms the presence of 
#recurrent periodicities along the ACL as suggested by the > 90 % CL frequencies 
#observed in the MTM



#* 4.2. WaverideR (non-detrended data) ----

#** 4.2.1. Continuous Wavelet Transform (CWT) ----

alb_Ti_wt <- analyze_wavelet(alb_Ti_ndet, dj = 1/100, 
                             lowerPeriod = 0.05, 
                             upperPeriod = 25, 
                             verbose = FALSE, 
                             omega_nr = 10
                             )

plot_wavelet(wavelet = alb_Ti_wt,  
             lowerPeriod = NULL, 
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE, 
             periodlab = "Period (metres)", 
             x_lab = "Depth (metres)", 
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = NULL, 
             add_points = NULL, 
             add_abline_h = NULL, 
             add_abline_v = NULL, 
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE, 
             demean_mtm = TRUE, 
             detrend_mtm = TRUE, 
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )

#The Continuous Wavelet Transfrom (CWT) is performed on the Ti non-detrended 
#series and exhibits similar results as the EHA analysis with a lot of power 
#between 0.7-1m




#** 4.2.2. Tracking (in depth domain) of the 173 kyr cycle ----

#The powerful 0.7-1m (supposed to correspond to the 173 kyr cycle) is tracked 
#along the ACL

#Tracking the wavelet ("if you want to track the wavelet, remove the "#" 
#in front of the 3 next command lines)

#alb_Ti_track <- track_period_wavelet(astro_cycle = 173, 
#                                    wavelet = alb_Ti_wt, 
#                                     n.levels = 100, 
#                                     periodlab = "Period (metres)", 
#                                     x_lab = "depth (metres)", 
#                                     palette_name = "rainbow", 
#                                     color_brewer = "grDevices", 
#                                     plot_horizontal = TRUE
#                                     )

#alb_Ti_track_comp <- completed_series(wavelet = alb_Ti_wt, 
#                                      tracked_curve = alb_Ti_track, 
#                                      period_up = 1, 
#                                      period_down = 0.75, 
#                                      extrapolate = TRUE, 
#                                      genplot = FALSE, 
#                                      keep_editable = FALSE
#                                      )

#alb_Ti_track_comp <- loess_auto(alb_Ti_track_comp)

#write.csv(alb_Ti_track_comp,"alb_Ti_track.csv")


# Loading the tracked curve
alb_Ti_track_comp <- read.csv("alb_Ti_track.csv")
alb_Ti_track_comp <- alb_Ti_track_comp[,c(2,3)]


plot_wavelet(wavelet = alb_Ti_wt, 
             lowerPeriod = NULL, 
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE, 
             periodlab = "Period (metres)", 
             x_lab = "Depth (metres)", 
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = cbind(alb_Ti_track_comp[,1],alb_Ti_track_comp[,2]),
             add_points = NULL, 
             add_abline_h = NULL, 
             add_abline_v = NULL, 
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE, 
             demean_mtm = TRUE, 
             detrend_mtm = TRUE, 
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )



#** 4.2.3. Depth to time domain convertion using the tracked cycle ----

alb_Ti_time <- curve2tune(data = alb_Ti_ndet, 
                          tracked_cycle_curve = alb_Ti_track_comp, 
                          tracked_cycle_period = 173, 
                          genplot = FALSE, 
                          keep_editable = FALSE
                          )

alb_Ti_time <- linterp(alb_Ti_time)  

alb_Ti_time_wt <- analyze_wavelet(alb_Ti_time, 
                                  dj = 1/100, 
                                  lowerPeriod = 1, 
                                  upperPeriod = 3000, 
                                  verbose = FALSE, 
                                  omega_nr = 10
                                  )

plot_wavelet(wavelet = alb_Ti_time_wt,
             lowerPeriod = NULL,
             upperPeriod = NULL,
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices",
             useRaster = TRUE,
             periodlab = "Period (kyr)",
             x_lab = "Time (kyr)", 
             keep_editable = FALSE,
             dev_new = F,
             add_lines = NULL,
             add_points = NULL,
             add_abline_h = c(20,32,100,110,125,173,405),
             add_abline_v = NULL,
             add_MTM_peaks = FALSE,
             add_data = TRUE,
             add_avg = TRUE,
             add_MTM = FALSE,
             demean_mtm = TRUE,
             detrend_mtm = TRUE,
             padfac_mtm = 5,
             tbw_mtm = 3,
             plot_horizontal = TRUE
             )

#The CWT performed in time domain shows a lot of power and/or prominent peaks 
#in the Milankovitch bands (precession: ~20 kyr, obliquity: ~32 kyr, 
#short eccentricity: ~100-125 kyr and 173 kyr) in line with the expectations 
#from the literature.

--------------------------------------------------------
  
#Step 5: Average Spectral Misfit (ASM) ----

freq <- c(1.151,1.458,1.548,1.647,1.746,1.964,2.093,5.228,5.704,5.853,
          6.240,7.103,7.202,10.317,11.101,12.55,12.659) 

target <- c(1/405,1/173,1/131.3,1/123.9,1/99.2,1/94.9,1/38.8,1/32.1,
            1/31.6,1/31.3,1/30.8, 1/24.4,1/20.2,1/19.8,1/19.3,1/16.8,1/16.7)

rayleigh <- 0.1

nyquist <- 15    

asm(freq, 
    target, 
    fper=NULL, 
    rayleigh, 
    nyquist, 
    sedmin=0.05, 
    sedmax=5, 
    numsed=500, 
    linLog=1, 
    iter=100000, 
    output=F, 
    genplot=T
    )

#The ASM analysis is performed confronting the > 90 % CL frequencies from the 
#MTM analysis with the theroretical duration of Milankovitch cycles.

#The ASM results in a sedimentation rate of 0.493 cm/ka, consistent with 
#the 0.4-0.6 cm/kyr expectation for Blatica during the Cambrian.

--------------------------------------------------------
  
#Step 6: Milankovitch bands filters in depth domain ----

#405 kyr cycle (1.75-2.25 m)
alb_Ti_405 <- taner(alb_Ti, xmax=5, fhigh=1/1.75, flow=1/2.25, demean = T)
#173 kyr cycle (0.76-1 m)
alb_Ti_173 <- taner(alb_Ti,xmax=5, fhigh=1/0.76, flow=1, demean = T)	
#short eccentricity band (0.42-0.72 m)
alb_Ti_100 <- taner(alb_Ti, xmax=5, fhigh=1/0.42, flow=1/0.72, demean = T)
#obliquity band (0.12-0.2 m)
alb_Ti_31 <- taner(alb_Ti, xmax=10, fhigh=1/0.12, flow=1/0.2, demean=T)	
#precession band (0.08-0.1 m)
alb_Ti_20 <- taner(alb_Ti, xmax=15, fhigh=1/0.08, flow=1/0.1, demean=T)	

#The filtering of the Milankovitch bands in depth domain allows to see the 
#imprint of the astronomical cycles on the Ti sedimentary record.

--------------------------------------------------------
  
#Step 7: Tuning of the Ti series ----

#* 7.1. Filtering of the 0.76-0.97 m band ----

Ti_173R <- taner(alb_Ti, xmax=2, fhigh=1/0.97, flow=1/0.76, demean=T)



#* 7.2. Identification of the trough of the filtered signal ----

Ti_173_min <- trough(Ti_173R, level = 200, genplot = T)

Ti_173_min <- Ti_173_min[,c(2,3)]

#write.csv(Ti_173_min,"Ti_173_rsp0.005_minTan_0.76-0.97m.csv")



#* 7.3. Trough = multiple of 173 kyr cycle ----

#Once the ".csv" file is created, go into it and for each location, 
#add a multiple of 173 for the trough value as follows:
  
#Location; Trough_Value

#74.127; 173
#75.002; 346
#75.852; 519
#76.697; 692
#…; …
#93.527; 4325



#* 7.4. Tuning of the 5mm spaced Ti detrended dataset ----

Ti_173_min <- read.csv("Ti_173_rsp0.005_minTan_0.76-0.97m.csv", sep=";")

Ti_173_tuned <- tune(alb_Ti, Ti_173_min, extrapolate=T, genplot=T, check=T, verbose=T)



#* 7.5. Evenly spaced dataset in time domain ----

#Interpolation of the tuned dataset every 1.067338 kyr (mean sampling interval) 
#to have an evenly space tune dataset

alb_Ti_time <- linterp(Ti_173_tuned, dt=1.067338, genplot=T)

--------------------------------------------------------
  
#Step 8: MTM and EHA in time domain ----

#* 8.1. MTM in time domain ----

alb_Ti_time_MTM<-mtm(alb_Ti_time, xmax=0.07, siglevel = 0.9, ar1 = T, tbw= 2, output = 1)
alb_Ti_time_MTM <- iso(dat=alb_Ti_time_MTM,xmin=0,xmax=0.07, genplot=F)
graphics.off()

plot(alb_Ti_time_MTM$Frequency, 
     alb_Ti_time_MTM$Power, 
     type = 'l', 
     xlab = "Frequency (cycles/kyr)", 
     ylab = "Variance"
     )

lines(alb_Ti_time_MTM$Frequency, 
      alb_Ti_time_MTM$AR1_fit, 
      type = 'l', 
      col = 'black', 
      lwd = 2
      )

lines(alb_Ti_time_MTM$Frequency, 
      alb_Ti_time_MTM$AR1_90_power, 
      type = 'l', 
      col = 'red', 
      lwd = 2
      )

lines(alb_Ti_time_MTM$Frequency, 
      alb_Ti_time_MTM$AR1_95_power, 
      type = 'l', 
      col = 'blue', 
      lwd = 2
      )

lines(alb_Ti_time_MTM$Frequency, 
      alb_Ti_time_MTM$AR1_99_power, 
      type = 'l', 
      col = 'green', 
      lwd = 2
      )

#The 2 pi-MTM on the tuned dataset shows the range of cyclicities in the signal 
#in time domain. Allowing to visually identify > 90 % CL periodicities.



#* 8.2. EHA in time domain ----

alb_Ti_time_EHA<- eha(alb_Ti_time, 
                      tbw=2, 
                      fmin=0.001, 
                      fmax=0.07, 
                      step=1, 
                      win=1100, 
                      demean=T, 
                      detrend=T, 
                      siglevel=0.90, 
                      sigID=F, 
                      ydir=1, 
                      output=1, 
                      pl=1, 
                      palette=1, 
                      centerZero=T, 
                      ncolors=200, 
                      genplot=4, 
                      verbose=T
                      )

#The EHA is performed on the Ti tuned series and confirms the presence of 
#recurrent periodicities in the Milankovitch bands along the ACL as suggested 
#by the > 90 % CL frequencies observed in the 2 pi-MTM

--------------------------------------------------------
  
#Step 9: Hilbert transform analysis between obliquity and the 173 kyr cycle ----

#Filter out the obliquity cycle
alb_Ti_time_obl <- taner(alb_Ti_time, xmax=1/20, fhigh=1/27, flow=1/40)

#Do the Hilbert transform to extract the amplitude
alb_Ti_time_obl_hilbert <- hilbert(alb_Ti_time_obl)

#Filter out the 173 kyr cycle from the amplitude modulation of the obliquity cycle
alb_Ti_time_obl_hilbert_173 <-taner(alb_Ti_time_obl_hilbert, 
                                    xmax=1/100, 
                                    fhigh=1/195, 
                                    flow=1/155, 
                                    detrend=TRUE
                                    )

#Filter out the 173 kyr obliquity cycle directly from the proxy record 
alb_Ti_time_173 <- taner(alb_Ti_time, 
                         xmax=1/50, 
                         fhigh=1/195, 
                         flow=1/155, 
                         demean=T
                         )
graphics.off()

#Plot the 173 kyr obliquity cycle extracted from the amplitude modulation of the 
#obliquity cycle vs the 173 kyr obliquity cycle directly extracted from the proxy
#record. 

plot(alb_Ti_time_173[,1],
     alb_Ti_time_173[,2]-mean(alb_Ti_time_173[,2]),
     type="l", 
     xlab ="Time (kyr)", 
     ylab = "Power"
     )

lines(alb_Ti_time_obl_hilbert_173[,1],
      alb_Ti_time_obl_hilbert_173[,2]-mean(alb_Ti_time_obl_hilbert_173[,2]), 
      col="red",
      lwd=2
      ) 

--------------------------------------------------------
  
#Step 10: Hilbert transform analysis between the 173 kyr and the 1.2 Myr cycles ----

#* 10.1. Tuning of the 1mm non-detrended Ti series ----

#Same tuning procedure as for step 7 but on the non resampled and non detrended 
#Ti dataset

#tuning of the non detrended 1mm Ti series
Ti_1mm_173R <- taner(alb_Ti_1mm, xmax=2, fhigh=1/0.97, flow=1/0.76, demean=T)
Ti_1mm_173_min <- trough(Ti_1mm_173R, level = 200, genplot = T)
Ti_1mm_173_min <- Ti_1mm_173_min[,c(2,3)]

#write.csv(Ti_1mm_173_min,"Ti_1mm_173_minTan_0.76-0.97m.csv")

#Before continuing, go to the .csv file and for each location, add a multiple of
#173 for the trough value as explain in step 7.3.  

Ti_1mm_173_min <- read.csv("Ti_1mm_173_minTan_0.76-0.97m.csv", sep=";")

Ti_1mm_173_tuned <- tune(alb_Ti_1mm, 
                         Ti_1mm_173_min, 
                         extrapolate=T, 
                         genplot=T, 
                         check=T, 
                         verbose=T
                         )

#Interpolation every 0.2131761 kyr (= mean sampling interval)
alb_Ti_1mm_time <- linterp(Ti_1mm_173_tuned, dt=0.2131761, genplot=T)



#* 10.2. Hilbert analysis ----

#Filter out the 173 kyr long obliquity cycle
alb_Ti_1mm_time_173 <- taner(alb_Ti_1mm_time, xmax=1/100, fhigh=1/155, flow=1/195)

#Do the Hilbert transform to extract the amplitude
alb_Ti_1mm_time_173_hilbert <- hilbert(alb_Ti_1mm_time_173)

#Filter out the 1.2 Myr cycle from the amplitude modulation of the 173 kyr cycle
alb_Ti_1mm_time_173_hilbert_1300 <- taner(alb_Ti_1mm_time_173_hilbert, 
                                          xmax=1/500, 
                                          fhigh=1/1300, 
                                          flow=1/1150, 
                                          detrend=T
                                          )

#Filter out the 1.2 Myr obliquity cycle directly from the proxy record 
alb_Ti_1mm_time_1300 <- taner(alb_Ti_1mm_time, 
                              xmax=1/500, 
                              fhigh=1/1300, 
                              flow=1/1150, 
                              demean=T
                              )
graphics.off()

#Plot the 1.2 Myr obliquity cycle extracted from the amplitude modulation of the
#173 kyr long obliquity cycle vs the 1.2 Myr obliquity cycle directly extracted
#from the proxy record. 

plot(alb_Ti_1mm_time_1300[,1],
     alb_Ti_1mm_time_1300[,2]-mean(alb_Ti_1mm_time_1300[,2]), 
     type="l", 
     xlab = "Time (kyr)", 
     ylab = "Power"
     )

lines(alb_Ti_1mm_time_173_hilbert_1300[,1],
      alb_Ti_1mm_time_173_hilbert_1300[,2]-mean(alb_Ti_1mm_time_173_hilbert_1300[,2]),
      col="red",
      lwd=2
      )

--------------------------------------------------------
  
#Step 11: Age-depth model and sedimentation rate based on 4 proxy comparison ----

#* 11.1. Age-Depth model (bandpass filtering) ----

#** 11.1.1. Minimal_tuning function ----

minimal_tuning <- function (data = NULL, pts = 5, cycle = 173, tune_opt = "max", 
                            output = 0, genplot = FALSE, keep_editable = FALSE) 
{
  astro_mindetect <- as.data.frame(data)
  astro_mindetect$min <- 0
  for (i in pts:(nrow(data) - pts)) {
    if ((data[i, 2] - data[(i + pts), 2] < 0) & (data[i, 
                                                      2] - data[(i - (pts - 1)), 2] < 0)) {
      astro_mindetect[i, 3] <- 1
    }
  }
  astro_mindetect_error_corr <- astro_mindetect
  astro_mindetect_error_corr <- astro_mindetect_error_corr[astro_mindetect_error_corr$min == 
                                                             1, ]
  astro_maxdetect <- as.data.frame(data)
  astro_maxdetect$max <- 0
  for (i in pts:(nrow(data) - pts)) {
    if ((data[i, 2] - data[(i + pts), 2] > 0) & (data[i, 
                                                      2] - data[(i - (pts - 1)), 2] > 0)) {
      astro_maxdetect[i, 3] <- 1
    }
  }
  astro_maxdetect_error_corr <- astro_maxdetect
  astro_maxdetect_error_corr <- astro_maxdetect_error_corr[astro_maxdetect_error_corr$max == 
                                                             1, ]
  max <- astro_maxdetect_error_corr
  colnames(max) <- c("A", "B", "C")
  min <- astro_mindetect_error_corr
  colnames(min) <- c("A", "B", "C")
  min[, 3] <- -1
  peaks <- rbind(max, min)
  peaks <- peaks[order(peaks[, 1]), ]
  i <- 1
  res_rownr <- nrow(peaks)
  while (i < res_rownr) {
    if ((i < res_rownr) & (peaks[i, 3] == peaks[(i + 1), 
                                                3])) {
      if ((i < res_rownr) & (peaks[i, 3] == 1 & peaks[(i + 
                                                       1), 3] == 1) & (peaks[i, 2] > peaks[(i + 1), 
                                                                                           2])) {
        peaks[(i + 1), ] <- NA
        peaks <- na.omit(peaks)
        res_rownr <- res_rownr - 1
      }
      if ((i < res_rownr) & (peaks[i, 3] == 1 & peaks[(i + 
                                                       1), 3] == 1) & (peaks[i, 2] < peaks[(i + 1), 
                                                                                           2])) {
        peaks[i, ] <- NA
        peaks <- na.omit(peaks)
        res_rownr <- res_rownr - 1
      }
      if ((i < res_rownr) & (peaks[i, 3] == -1 & peaks[(i + 
                                                        1), 3] == -1) & (peaks[i, 2] < peaks[(i + 1), 
                                                                                             2])) {
        peaks[(i + 1), ] <- NA
        peaks <- na.omit(peaks)
        res_rownr <- res_rownr - 1
      }
      if ((i < res_rownr) & (peaks[i, 3] == -1 & peaks[(i + 
                                                        1), 3] == -1) & (peaks[i, 2] > peaks[(i + 1), 
                                                                                             2])) {
        peaks[i, ] <- NA
        peaks <- na.omit(peaks)
        res_rownr <- res_rownr - 1
      }
    }
    if ((peaks[i, 3] != peaks[(i + 1), 3]) | is.na(peaks[i, 
                                                         3] != peaks[(i + 1), 3])) {
      i <- i + 1
    }
  }
  if (tune_opt == "min") {
    peaks_min <- peaks[peaks[, 3] < 0, ]
    dist <- peaks_min[2:(nrow(peaks_min)), ] - peaks_min[1:(nrow(peaks_min) - 
                                                              1), ]
    sed_rate <- (dist[, 1] * 100)/cycle
    sed_rate <- cbind(sed_rate, peaks_min[1:(nrow(peaks_min) - 
                                               1), 1], peaks_min[2:(nrow(peaks_min)), 1])
  }
  if (tune_opt == "max") {
    peaks_min <- peaks[peaks[, 3] > 0, ]
    dist <- peaks_min[2:(nrow(peaks_min)), ] - peaks_min[1:(nrow(peaks_min) - 
                                                              1), ]
    sed_rate <- (dist[, 1] * 100)/cycle
    sed_rate <- cbind(sed_rate, peaks_min[1:(nrow(peaks_min) - 
                                               1), 1], peaks_min[2:(nrow(peaks_min)), 1])
  }
  if (tune_opt == "minmax") {
    peaks_min <- peaks
    dist <- peaks_min[2:(nrow(peaks_min)), ] - peaks_min[1:(nrow(peaks_min) - 
                                                              1), ]
    sed_rate <- (dist[, 1] * 100)/(cycle/2)
    sed_rate <- cbind(sed_rate, peaks_min[1:(nrow(peaks_min) - 
                                               1), 1], peaks_min[2:(nrow(peaks_min)), 1])
  }
  top <- c(sed_rate[1, 1], data[1, 1], sed_rate[1, 2])
  bot <- c(sed_rate[nrow(sed_rate), 1], sed_rate[nrow(sed_rate), 
                                                 3], data[nrow(data), 1])
  sed_rate <- rbind(top, sed_rate, bot)
  data[, 3] <- NA
  p <- 1
  for (i in 1:nrow(data)) {
    if (data[i, 1] < sed_rate[p, 2]) {
      data[i, 3] <- sed_rate[p, 1]
    }
    if (data[i, 1] == sed_rate[p, 2] & p + 1 <= nrow(sed_rate)) {
      data[i, 3] <- (sed_rate[p, 1] + sed_rate[(p + 1), 
                                               1])/2
    }
    if (p > nrow(sed_rate)) {
      p <- nrow(sed_rate)
    }
    if (p == nrow(sed_rate)) {
      data[i, 3] <- sed_rate[nrow(sed_rate), 1]
    }
    if (data[i, 1] > sed_rate[p, 2]) {
      data[i, 3] <- sed_rate[p, 1]
    }
    if (data[i, 1] == sed_rate[p, 3] & p + 1 <= nrow(sed_rate)) {
      p <- p + 1
    }
  }
  tracked_cycle_curve <- data[, c(1, 3)]
  sedrates <- data.frame(tracked_cycle_curve)
  dat <- as.matrix(tracked_cycle_curve)
  dat <- na.omit(dat)
  dat <- dat[order(dat[, 1], na.last = NA, decreasing = F), 
  ]
  npts <- length(dat[, 1])
  start <- dat[1, 1]
  end <- dat[length(dat[, 1]), 1]
  x1 <- dat[1:(npts - 1), 1]
  x2 <- dat[2:(npts), 1]
  dx = x2 - x1
  dt = median(dx)
  xout <- seq(start, end, by = dt)
  npts <- length(xout)
  interp <- approx(dat[, 1], dat[, 2], xout, method = "linear", 
                   n = npts)
  sedrates <- as.data.frame(interp)
  npts <- length(sedrates[, 1])
  sedrates[1] = sedrates[1] * 100
  sedrates[2] = 1/sedrates[2]
  dx = sedrates[2, 1] - sedrates[1, 1]
  midptx = (sedrates[2:npts, 1] + sedrates[1:(npts - 1), 1])/2
  slope = (sedrates[2:npts, 2] - sedrates[1:(npts - 1), 2])/dx
  yint = sedrates[2:npts, 2] - (slope * sedrates[2:npts, 1])
  midpty = (slope * midptx) + yint
  hsum = cumsum(midpty * dx)
  hsum = append(0, hsum)
  out = data.frame(cbind(sedrates[, 1]/100, hsum))
  data[, 4] <- out[, 2]
  colnames(data) <- c("depth", "proxy", "cm/kyr", "time")
  if (output == 0) {
    data <- data
    if (genplot == TRUE) {
      if (keep_editable == FALSE) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
      }
      layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 4, 
                              ncol = 1)
      graphics::layout(mat = layout.matrix, heights = c(1), 
                       widths = c(1))
      par(mar = c(4, 4, 1, 1))
      plot(x = data[, 1], y = data[, 2], type = "l", main = "Data depth domain", 
           xlab = "meters", ylab = "proxy")
      plot(x = data[, 1], y = data[, 3], type = "l", xlab = "meters", 
           ylab = "cm/kyr (ka)", main = "sedimentation rate plot")
      plot(data[, 1], data[, 4], type = "l", xlab = "meters", 
           ylab = "Time (ka)", main = "Data time domain")
      plot(data[, 4], data[, 2], type = "l", xlab = "time (ka)", 
           ylab = "proxy", main = "Data time domain")
    }
  }
  if (output == 1) {
    data <- data[, c(1, 3)]
    if (genplot == TRUE) {
      if (keep_editable == FALSE) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
      }
      layout.matrix <- matrix(c(1), nrow = 1, ncol = 1)
      graphics::layout(mat = layout.matrix, heights = c(1), 
                       widths = c(1))
      par(mar = c(4, 4, 1, 1))
      plot(x = data[, 1], y = data[, 2], type = "l", xlab = "meters", 
           ylab = "cm/kyr (ka)", main = "sedimentation rate plot")
    }
  }
  if (output == 2) {
    data <- data[, c(1, 4)]
    if (genplot == TRUE) {
      if (keep_editable == FALSE) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
      }
      layout.matrix <- matrix(c(1), nrow = 1, ncol = 1)
      graphics::layout(mat = layout.matrix, heights = c(1), 
                       widths = c(1))
      par(mar = c(4, 4, 1, 1))
      plot(data[, 1], data[, 2], type = "l", xlab = "meters", 
           ylab = "Time (ka)", main = "Data time domain")
    }
  }
  return(data)
}



#** 11.1.2. Taner filtering of the detrital elements in depth domain ----

#*** 11.1.2.1. Titanium ----

alb_Ti <- cbind(Alb$DepthAdj,Alb$Ti)
alb_Ti <- na.omit(alb_Ti)
alb_Ti[!is.finite(alb_Ti)] <- NA
alb_Ti <- na.omit(alb_Ti)
alb_Ti <- iso(dat=alb_Ti, xmin=73.453, xmax=93.609) 
alb_Ti <- linterp(alb_Ti, dt=0.005, genplot=T)

Ti_173R <- taner(alb_Ti,xmax=2,fhigh=1/0.97,flow=1/0.75,demean=TRUE)

Ti_173R_min_tun <-    minimal_tuning(data = Ti_173R,
                                     pts = 20,
                                     cycle = 173,
                                     tune_opt = "min", 
                                     output = 0, 
                                     genplot = TRUE,
                                     keep_editable = FALSE
                                     )

colnames(Ti_173R_min_tun)<-c("Ti Depth","Ti","Ti SR","Ti time")



#*** 11.1.2.2. Silicon ----

alb_Si <- cbind(Alb$DepthAdj,Alb$Si)
alb_Si <- na.omit(alb_Si)
alb_Si[!is.finite(alb_Si)] <- NA
alb_Si <- na.omit(alb_Si)
alb_Si <- iso(dat=alb_Si, xmin=73.453, xmax=93.609) 
alb_Si <- linterp(alb_Si, dt=0.005, genplot=T)

Si_173R <- taner(alb_Si,xmax=2,fhigh=1/0.97,flow=1/0.75,demean=TRUE)

Si_173R_min_tun <-    minimal_tuning(data = Si_173R,
                                     pts = 20,
                                     cycle = 173,
                                     tune_opt = "min",
                                     output = 0,
                                     genplot = TRUE,
                                     keep_editable = FALSE
                                     )

colnames(Si_173R_min_tun)<-c("Si Depth","Si","Si SR","Si time")



#*** 11.1.2.3. Aluminium ----

alb_Al <- cbind(Alb$DepthAdj,Alb$Al)
alb_Al <- na.omit(alb_Al)
alb_Al[!is.finite(alb_Al)] <- NA
alb_Al <- na.omit(alb_Al)
alb_Al <- iso(dat=alb_Al, xmin=73.453, xmax=93.609) 
alb_Al <- linterp(alb_Al, dt=0.005, genplot=T)

Al_173R <- taner(alb_Al,xmax=2,fhigh=1/0.97,flow=1/0.75,demean=TRUE)

Al_173R_min_tun <-    minimal_tuning(data = Al_173R,
                                     pts = 20,
                                     cycle = 173,
                                     tune_opt = "min",
                                     output = 0,
                                     genplot = TRUE,
                                     keep_editable = FALSE
                                     )

colnames(Al_173R_min_tun)<-c("Al Depth","Al","Al SR","Al time")




#*** 11.1.2.4. Potassium ----

alb_K <- cbind(Alb$DepthAdj,Alb$K)
alb_K <- na.omit(alb_K)
alb_K[!is.finite(alb_K)] <- NA
alb_K <- na.omit(alb_K)
alb_K <- iso(dat=alb_K, xmin=73.453, xmax=93.609) 
alb_K <- linterp(alb_K, dt=0.005, genplot=T)

K_173R <- taner(alb_K,xmax=2,fhigh=1/0.97,flow=1/0.75,demean=TRUE)

K_173R_min_tun <-    minimal_tuning(data = K_173R,
                                    pts = 20,
                                    cycle = 173,    
                                    tune_opt = "min",
                                    output = 0,
                                    genplot = TRUE,
                                    keep_editable = FALSE
                                    )

colnames(K_173R_min_tun)<-c("K Depth","K","K SR","K time")



#*** 11.1.2.5. Zirconium ----

alb_Zr <- cbind(Alb$DepthAdj,Alb$Zr)
alb_Zr <- na.omit(alb_Zr)
alb_Zr[!is.finite(alb_Zr)] <- NA
alb_Zr <- na.omit(alb_Zr)
alb_Zr <- iso(dat=alb_Zr, xmin=73.453, xmax=93.609) 
alb_Zr <- linterp(alb_Zr, dt=0.005, genplot=T)

Zr_173R <- taner(alb_Zr,xmax=2,fhigh=1/0.97,flow=1/0.75,demean=TRUE)

Zr_173R_min_tun <-    minimal_tuning(data = Zr_173R,pts = 20,cycle = 173,
                                     tune_opt = "min",output = 0,genplot = TRUE,keep_editable = FALSE)

colnames(Zr_173R_min_tun)<-c("Zr Depth","Zr","Zr SR","Zr time")




#** 11.1.3. File containing the filtering of all the detrital element ----

#All detrital
Det_173_min_tun <- cbind(Ti_173R_min_tun,
                         Si_173R_min_tun,
                         Al_173R_min_tun,
                         K_173R_min_tun,
                         Zr_173R_min_tun
                         )
head((Det_173_min_tun))
graphics.off()



#** 11.1.4. Sedimentation rate ----

plot(Det_173_min_tun[,1],
     Det_173_min_tun[,3],
     type="l",
     xlab = "Depth adjusted (m)", 
     ylab = "Sed rate (cm/kyr)",
     ylim=c(0.38,0.55),
     lwd = 2
     )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,7],
      col="red",
      lwd=2
      )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,11],
      col="green",
      lwd=2
      )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,15],
      col="blue",
      lwd=2
      )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,19],
      col="grey",
      lwd=2
      )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,19]-Det_173_min_tun[,19]+0.493,
      col="red",
      lwd=2, 
      lty=3
      )

#Zr sedimentation rate is significantly different compared to the other detrital
#element and will not be used.

plot(Det_173_min_tun[,1],
     Det_173_min_tun[,3],
     type="l",
     xlab = "Depth adjusted (m)", 
     ylab = "Sed rate (cm/kyr)",
     ylim=c(0.38,0.55),
     lwd = 2
     )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,7],
      col="red",
      lwd=2
      )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,11],
      col="green",
      lwd=2
      )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,15],
      col="blue",
      lwd=2
      )

lines(Det_173_min_tun[,1],
      Det_173_min_tun[,19]-Det_173_min_tun[,19]+0.493,
      col="red",
      lwd=2, 
      lty=3
      )

legend("bottomleft",
       legend = c("Ti (0.75-0.97 m filter)", 
                  "Si (0.75-0.97 m filter)", 
                  "Al (0.75-0.97 m filter)", 
                  "K (0.75-0.97 m filter)", 
                  "ASM = 0.493 cm/kyr (theoretical)"),
       col = c("black", "red", "green", "blue", "red"),
       lty = c(1, 1, 1, 1,3),
       lwd = c(2, 2, 2, 2,2),
       pch = c(NA, NA, NA, NA,NA),
       pt.cex = 2,
       bty = "n",         # No box around legend
       cex = 0.7)  



#** 11.1.5. Mean and standard deviation calculation ----

#Smoothing the sed rate and turning it into frequency to calculate the mean and the
#standard deviation of the age model
Ti_173_freq<-noLow(cbind(Det_173_min_tun[,c(1)],1/(Det_173_min_tun[,c(3)]*173)),
                   0.05, 
                   output = 2
                   )

Si_173_freq<-noLow(cbind(Det_173_min_tun[,c(1)],1/(Det_173_min_tun[,c(7)]*173)),
                   0.05,
                   output = 2
                   )

Al_173_freq<-noLow(cbind(Det_173_min_tun[,c(1)],1/(Det_173_min_tun[,c(11)]*173)),
                   0.05,
                   output = 2
                   )

K_173_freq<-noLow(cbind(Det_173_min_tun[,c(1)],1/(Det_173_min_tun[,c(15)]*173)),
                  0.05,
                  output = 2
                  )

Det_173_freq <- cbind(Det_173_min_tun[,1],
                      Ti_173_freq[,2],
                      Si_173_freq[,2],
                      Al_173_freq[,2], 
                      K_173_freq[,2]
                      )

head(Det_173_freq)

Det_173_freq<-as.data.frame(Det_173_freq) 
Det_173_freq$mean<-rowMeans(Det_173_freq[,c(2:5)], na.rm = FALSE)

library(matrixStats)
Det_173_freq$sds<-rowSds(as.matrix(Det_173_freq[,c(2:5)]), na.rm = FALSE)
colnames(Det_173_freq)<- cbind("Depth",
                               "Ti freq",
                               "Si freq",
                               "Al freq", 
                               "K freq", 
                               "mean", 
                               "sds"
                               )
graphics.off()

plot(Det_173_freq[,1],
     1/(Det_173_freq[,2]),
     type="l",
     col="black",
     lwd= 2,
     xlab = "Depth adjusted (m)", 
     ylab = "Period (cm)", 
     ylim = c(65,95)
     )

lines(Det_173_freq[,1],
      1/(Det_173_freq[,3]),
      col="red",
      lwd=2
      )

lines(Det_173_freq[,1],
      1/(Det_173_freq[,4]),
      col="green",
      lwd=2
      )

lines(Det_173_freq[,1],
      1/(Det_173_freq[,5]),
      col="blue",
      lwd=2
      )

#Graphical visualisation of the uncertaineties
depth<-Det_173_freq[,1]
curve_mean <- Det_173_freq[,6]
curve_sd <- Det_173_freq[,7]

#Add polygon for mean ± SD 
polygon(c(depth, rev(depth)),
        1/c(curve_mean + curve_sd, rev(curve_mean - curve_sd)),
        col = rgb(0.7, 0.7, 0.7, 0.4), 
        border = NA
        )

#Optionally overlay mean line
lines(depth, 1/curve_mean, col = "yellow", lwd = 2)
legend("bottomleft",
       legend = c("Ti", "Si", "K", "Al", "Mean", "Mean ± SD"),
       col = c("black", "red", "blue", "green","yellow", rgb(0.7, 0.7, 0.7, 0.8)),
       lty = c(1, 1, 1, 1, 1,NA),
       lwd = c(1, 1, 1, 1, 1,NA),
       pch = c(NA, NA, NA, NA,NA, 15),
       pt.cex = 2,
       bty = "n",         # No box around legend
       cex = 1
       )         # Larger legend text




#** 11.1.6. Building of the Age model ----

#Calculating the minimal, mean and maximal sed rates
mean_sed<-freq2sedrate(cbind(Det_173_freq[,1],
                             (Det_173_freq[,6])),
                       period=17300
                       )

max_sed<-freq2sedrate(cbind(Det_173_freq[,1],
                            (Det_173_freq[,6]+2*Det_173_freq[,7])),
                      period=17300
                      )

min_sed<-freq2sedrate(cbind(Det_173_freq[,1],
                            (Det_173_freq[,6]-2*Det_173_freq[,7])),
                      period=17300
                      )

#Sed rate to time and creation fo the Age Model
Mean_AgeModel<-sedrate2time(mean_sed)
Max_AgeModel<-sedrate2time(max_sed)
Min_AgeModel<-sedrate2time(min_sed)

plot(Mean_AgeModel,type="l", lwd=2)
lines(Max_AgeModel,col="green",lwd=2)
lines(Min_AgeModel,col="blue",lwd=2)

Age_Model <- cbind(Mean_AgeModel,Max_AgeModel[,2],Min_AgeModel[,2])
colnames(Age_Model)<- cbind("Adjusted depth (m)",
                            "Mean age (ka)",
                            "Max age (ka)",
                            "Min age (ka)"
                            )

#Interpolation of the Age-Model every 1 mm to have an age-depth correlation for 
#the entire Ti dataset with the original sampling step of 1mm
Mean_Age_Model_1mm <- linterp(Mean_AgeModel, dt=0.001, genplot=T)
Max_Age_Model_1mm <- linterp(Max_AgeModel, dt=0.001, genplot=T)
Min_Age_Model_1mm <- linterp(Min_AgeModel, dt=0.001, genplot=T)

Age_Model_1mm <- cbind(Mean_Age_Model_1mm,Max_Age_Model_1mm[,2],Min_Age_Model_1mm[,2])
colnames(Age_Model)<- cbind("Adjusted depth (m)",
                            "Mean age (ka)",
                            "Max age (ka)",
                            "Min age (ka)"
                            )
graphics.off()
#Saving of the age models as .csv files
#  write.csv(Age_Model_1mm,"AgeModel_rsp1mm_filter0.75-0.97m.csv")



#* 11.2. Age-Depth model (Wavelet tracking; preferred) ----

#** 11.2.1. Tracking of the 173 kyr for Ti, Al, Si and K ----

#As shown in the time series analysis (MTM, EHA and CWT), the 173 kyr cycle is 
#comprise between 0.75 and 0.9 m. To avoid overpass this range when tracking the 
#Wavelet, especially with periodicities lower than 0.74 m, that could be linked 
#to the short eccentricity (100-135 kyr) band, we used a lowpass filtering of 
#the tracked Wavelet in accordance with the general trend observed in the CWT.

#*** 11.2.1.1. Titanium ----

alb_Ti <- cbind(Alb$DepthAdj,Alb$Ti)
alb_Ti <- na.omit(alb_Ti)
alb_Ti[!is.finite(alb_Ti)] <- NA
alb_Ti <- na.omit(alb_Ti)
alb_Ti <- iso(dat=alb_Ti, xmin=73.453, xmax=100)
alb_Ti_ndet <-  linterp(alb_Ti, dt=0.005, genplot=T)

#Waverider
alb_Ti_wt <- analyze_wavelet(alb_Ti_ndet, 
                             dj = 1/100, 
                             lowerPeriod = 0.5, 
                             upperPeriod = 5, 
                             verbose = FALSE, 
                             omega_nr = 10
                             )

plot_wavelet(wavelet = alb_Ti_wt,  
             lowerPeriod = NULL,
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE, 
             periodlab = "Period (metres)", 
             x_lab = "Depth (metres)", 
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = NULL, 
             add_points = NULL, 
             add_abline_h = NULL, 
             add_abline_v = NULL, 
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE, 
             demean_mtm = TRUE, 
             detrend_mtm = TRUE, 
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )

#track the period (m) of the 173 kyr cycle
#    alb_Ti_track_WR <- track_period_wavelet(astro_cycle = 173,  
#                                            wavelet = alb_Ti_wt, 
#                                            n.levels = 100, 
#                                            periodlab = "Period (metres)", 
#                                            x_lab = "depth (metres)", 
#                                            palette_name = "rainbow",
#                                            color_brewer = "grDevices",
#                                            plot_horizontal = TRUE
#                                            )

#   alb_Ti_track_WR_comp <- completed_series(wavelet = alb_Ti_wt, 
#                                            tracked_curve = alb_Ti_track_WR, 
#                                            period_up = 1, 
#                                            period_down = 0.70, 
#                                            extrapolate = TRUE, 
#                                            genplot = FALSE, 
#                                            keep_editable = FALSE
#                                            )

#    alb_Ti_track_WR_comp[alb_Ti_track_WR_comp[,2]<0.74,2]<-0.74 

#    alb_Ti_track_WR_comp <- loess_auto(alb_Ti_track_WR_comp)

#    write.csv(alb_Ti_track_WR_comp,"alb_Ti_track_WR.csv")

# Loading the tracked curve
alb_Ti_track_WR_comp <- read.csv("alb_Ti_track_WR.csv")
alb_Ti_track_WR_comp <- alb_Ti_track_WR_comp[,c(2,3)]


plot_wavelet(wavelet = alb_Ti_wt, 
             lowerPeriod = NULL, 
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE, 
             periodlab = "Period (metres)", 
             x_lab = "Depth (metres)", 
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = cbind(alb_Ti_track_WR_comp[,1],
                               alb_Ti_track_WR_comp[,2]), 
             add_points = NULL, 
             add_abline_h = NULL, 
             add_abline_v = NULL, 
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE, 
             demean_mtm = TRUE, 
             detrend_mtm = TRUE, 
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )



#*** 11.2.1.2. Silicon ----

alb_Si <- cbind(Alb$DepthAdj,Alb$Si)
alb_Si <- na.omit(alb_Si)
alb_Si[!is.finite(alb_Si)] <- NA
alb_Si <- na.omit(alb_Si)
alb_Si <- iso(dat=alb_Si, xmin=73.453, xmax=100)
alb_Si_ndet <-  linterp(alb_Si, dt=0.005, genplot=T)

#Waverider
alb_Si_wt <- analyze_wavelet(alb_Si_ndet, 
                             dj = 1/100, 
                             lowerPeriod = 0.5, 
                             upperPeriod = 5, 
                             verbose = FALSE, 
                             omega_nr = 10
                             )

plot_wavelet(wavelet = alb_Si_wt,  
             lowerPeriod = NULL, 
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE,
             periodlab = "Period (metres)",
             x_lab = "Depth (metres)",
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = NULL,
             add_points = NULL, 
             add_abline_h = NULL, 
             add_abline_v = NULL,
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE, 
             demean_mtm = TRUE,
             detrend_mtm = TRUE, 
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )

#track the period (m) of the 173 kyr cycle
#    alb_Si_track_WR <- track_period_wavelet(astro_cycle = 173,
#                                            wavelet = alb_Si_wt,
#                                            n.levels = 100,
#                                            periodlab = "Period (metres)",
#                                            x_lab = "depth (metres)",
#                                            palette_name = "rainbow",
#                                            color_brewer = "grDevices",
#                                            plot_horizontal = TRUE
#                                            )

#   alb_Si_track_WR_comp <- completed_series(wavelet = alb_Si_wt, 
#                                            tracked_curve = alb_Si_track_WR,
#                                            period_up = 1, 
#                                            period_down = 0.70, 
#                                            extrapolate = TRUE, 
#                                            genplot = FALSE, 
#                                            keep_editable = FALSE
#                                            )

#    alb_Si_track_WR_comp[alb_Si_track_WR_comp[,2]<0.74,2]<-0.74 

#To keep a similar trend in the Si data in the 73-80 m interval, when tracking,
#we deliberately went for a higher band (above 0.9 m) and we implemented a
#highpass filter for all the values above 0.86m to keep a similar trend between
#Si and the other detrital elements

#     alb_Si_track_WR_comp[alb_Si_track_WR_comp[,2]>0.86,2]<-0.81 

#     alb_Si_track_WR_comp <- loess_auto(alb_Si_track_WR_comp)

#     write.csv(alb_Si_track_WR_comp,"alb_Si_track_WR.csv")


# Loading the tracked curve
alb_Si_track_WR_comp <- read.csv("alb_Si_track_WR.csv")
alb_Si_track_WR_comp <- alb_Si_track_WR_comp[,c(2,3)]

plot_wavelet(wavelet = alb_Si_wt, 
             lowerPeriod = NULL, 
             upperPeriod = NULL, 
             n.levels = 100,
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE, 
             periodlab = "Period (metres)",
             x_lab = "Depth (metres)", 
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = cbind(alb_Si_track_WR_comp[,1],
                               alb_Si_track_WR_comp[,2]), 
             add_points = NULL, 
             add_abline_h = NULL, 
             add_abline_v = NULL, 
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE,
             demean_mtm = TRUE, 
             detrend_mtm = TRUE,
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )



#*** 11.2.1.3. Aluminium ----

alb_Al <- cbind(Alb$DepthAdj,Alb$Al)
alb_Al <- na.omit(alb_Al)
alb_Al[!is.finite(alb_Al)] <- NA
alb_Al <- na.omit(alb_Al)
alb_Al <- iso(dat=alb_Al, xmin=73.453, xmax=100)
alb_Al_ndet <-  linterp(alb_Al, dt=0.005, genplot=T)

#Waverider
alb_Al_wt <- analyze_wavelet(alb_Al_ndet, 
                             dj = 1/100, 
                             lowerPeriod = 0.5, 
                             upperPeriod = 5, 
                             verbose = FALSE, 
                             omega_nr = 10
                             )

plot_wavelet(wavelet = alb_Al_wt,  
             lowerPeriod = NULL, 
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE, 
             periodlab = "Period (metres)", 
             x_lab = "Depth (metres)", 
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = NULL, 
             add_points = NULL, 
             add_abline_h = NULL,
             add_abline_v = NULL,
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE, 
             demean_mtm = TRUE, 
             detrend_mtm = TRUE, 
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )

#track the period (m) of the 173 kyr cycle
#   alb_Al_track_WR <- track_period_wavelet(astro_cycle = 173,  
#                                           wavelet = alb_Al_wt, 
#                                           n.levels = 100, 
#                                           periodlab = "Period (metres)", 
#                                           x_lab = "depth (metres)",
#                                           palette_name = "rainbow", 
#                                           color_brewer = "grDevices",
#                                           plot_horizontal = TRUE
#                                           )

#   alb_Al_track_WR_comp <- completed_series(wavelet = alb_Al_wt, 
#                                            tracked_curve = alb_Al_track_WR,
#                                            period_up = 1, 
#                                            period_down = 0.70, 
#                                            extrapolate = TRUE, 
#                                            genplot = FALSE, 
#                                            keep_editable = FALSE
#                                            )

#     alb_Al_track_WR_comp[alb_Al_track_WR_comp[,2]<0.74,2]<-0.74 

#     alb_Al_track_WR_comp <- loess_auto(alb_Al_track_WR_comp)

#     write.csv(alb_Al_track_WR_comp,"alb_Al_track_WR.csv")


# Loading the tracked curve
alb_Al_track_WR_comp <- read.csv("alb_Al_track_WR.csv")
alb_Al_track_WR_comp <- alb_Al_track_WR_comp[,c(2,3)]

plot_wavelet(wavelet = alb_Al_wt, 
             lowerPeriod = NULL, 
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices",
             useRaster = TRUE, 
             periodlab = "Period (metres)", 
             x_lab = "Depth (metres)",
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = cbind(alb_Al_track_WR_comp[,1],
                               alb_Al_track_WR_comp[,2]), 
             add_points = NULL, 
             add_abline_h = NULL,
             add_abline_v = NULL, 
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE, 
             demean_mtm = TRUE,
             detrend_mtm = TRUE, 
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )



#*** 11.2.1.4. Potassium ----

alb_K <- cbind(Alb$DepthAdj,Alb$K)
alb_K <- na.omit(alb_K)
alb_K[!is.finite(alb_K)] <- NA
alb_K <- na.omit(alb_K)
alb_K <- iso(dat=alb_K, xmin=73.453, xmax=100)
alb_K_ndet <-  linterp(alb_K, dt=0.005, genplot=T)

#Waverider
alb_K_wt <- analyze_wavelet(alb_K_ndet, 
                            dj = 1/100, 
                            lowerPeriod = 0.5, 
                            upperPeriod = 5, 
                            verbose = FALSE,
                            omega_nr = 10
                            )

plot_wavelet(wavelet = alb_K_wt,  
             lowerPeriod = NULL, 
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE, 
             periodlab = "Period (metres)", 
             x_lab = "Depth (metres)", 
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = NULL, 
             add_points = NULL, 
             add_abline_h = NULL,
             add_abline_v = NULL, 
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE,
             add_MTM = FALSE,
             demean_mtm = TRUE, 
             detrend_mtm = TRUE,
             padfac_mtm = 5, 
             tbw_mtm = 3,
             plot_horizontal = TRUE
             )

#track the period (m) of the 173 kyr cycle
#alb_K_track_WR <- track_period_wavelet(astro_cycle = 173, 
#                                       wavelet = alb_K_wt, 
#                                       n.levels = 100, 
#                                       periodlab = "Period (metres)",
#                                       x_lab = "depth (metres)",
#                                       palette_name = "rainbow", 
#                                       color_brewer = "grDevices",
#                                       plot_horizontal = TRUE
#                                       )

#   alb_K_track_WR_comp <- completed_series(wavelet = alb_K_wt, 
#                                           tracked_curve = alb_K_track_WR,
#                                           period_up = 1, 
#                                           period_down = 0.70, 
#                                           extrapolate = TRUE, 
#                                           genplot = FALSE, 
#                                           keep_editable = FALSE
#                                           )

#    alb_K_track_WR_comp[alb_K_track_WR_comp[,2]<0.74,2]<-0.74 

#To keep a similar trend in the K data in the 73-80 m interval, when tracking,
#we implemented a highpass filter for all the values above 0.875m to keep a
#similar trend between K and the other detrital elements

#     alb_K_track_WR_comp[alb_K_track_WR_comp[,2]>0.875,2]<-0.875 

#     alb_K_track_WR_comp <- loess_auto(alb_K_track_WR_comp)

#     write.csv(alb_K_track_WR_comp,"alb_K_track_WR.csv")


# Loading the tracked curve
alb_K_track_WR_comp <- read.csv("alb_K_track_WR.csv")
alb_K_track_WR_comp <- alb_K_track_WR_comp[,c(2,3)]

plot_wavelet(wavelet = alb_K_wt, 
             lowerPeriod = NULL,
             upperPeriod = NULL, 
             n.levels = 100, 
             palette_name = "rainbow",
             color_brewer = "grDevices", 
             useRaster = TRUE, 
             periodlab = "Period (metres)", 
             x_lab = "Depth (metres)",
             keep_editable = FALSE, 
             dev_new = F, 
             add_lines = cbind(alb_K_track_WR_comp[,1],
                               alb_K_track_WR_comp[,2]), 
             add_points = NULL, 
             add_abline_h = NULL, 
             add_abline_v = NULL, 
             add_MTM_peaks = FALSE, 
             add_data = TRUE, 
             add_avg = TRUE, 
             add_MTM = FALSE,
             demean_mtm = TRUE, 
             detrend_mtm = TRUE, 
             padfac_mtm = 5, 
             tbw_mtm = 3, 
             plot_horizontal = TRUE
             )
graphics.off()


#*** 11.2.1.5. Depth model resulting from the 173 kyr cycle tracking of the 4 detrital elements ----

alb_Ti_track_comp1 <- read.csv("alb_Ti_track_WR.csv")
alb_Si_track_comp1 <- read.csv("alb_Si_track_WR.csv")
alb_Al_track_comp1 <- read.csv("alb_Al_track_WR.csv")
alb_K_track_comp1 <- read.csv("alb_K_track_WR.csv")

plot(alb_Ti_track_comp1[,2],
     100*(alb_Ti_track_comp1[,3]),
     type="l",col="black",
     lwd= 2,
     xlab = "Depth adjusted (m)", 
     ylab = "Period (cm)", 
     ylim = c(70,92)
     )

lines(alb_Si_track_comp1[,2],
      100*(alb_Si_track_comp1[,3]),
      col="red",
      lwd=2
      )

lines(alb_Al_track_comp1[,2],
      100*(alb_Al_track_comp1[,3]),
      col="green",
      lwd=2
      )

lines(alb_K_track_comp1[,2],
      100*(alb_K_track_comp1[,3]),
      col="blue",
      lwd=2
      )

track_comp_1 <- cbind(alb_Ti_track_comp1[,3],
                      alb_Si_track_comp1[,3],
                      alb_Al_track_comp1[,3],
                      alb_K_track_comp1[,3]
                      )

depth<-alb_Ti_track_comp1[,2]
curve_mean <- rowMeans(100*(track_comp_1))
curve_sd <-  rowSds(100*(track_comp_1))

#Add polygon for mean ± SD 
polygon(c(depth, rev(depth)),
        c(curve_mean + curve_sd, rev(curve_mean - curve_sd)),
        col = rgb(0.7, 0.7, 0.7, 0.4), border = NA)

#Optionally overlay mean line
lines(depth, curve_mean, col = "yellow", lwd = 2)

legend("bottomleft",
       legend = c("Ti", "Si", "K", "Al", "Mean", "Mean ± SD"),
       col = c("black", "red", "blue", "green","yellow", rgb(0.7, 0.7, 0.7, 0.8)),
       lty = c(1, 1, 1, 1, 1,NA),
       lwd = c(2, 2, 2, 2, 2,NA),
       pch = c(NA, NA, NA, NA,NA, 15),
       pt.cex = 2,
       bty = "n",         # No box around legend
       cex = 0.75
       )         # Larger legend text



#** 11.2.2. Age model ----

#Age model
track_comp_1_mean <- rowMeans(1/track_comp_1)

track_comp_1_sds <- rowSds(1/track_comp_1)

track_comp_1_plus_2Sd <- cbind(alb_Ti_track_comp1[,2],
                               1/(track_comp_1_mean+(2*track_comp_1_sds)))

track_comp_1_min_2Sd <- cbind(alb_Ti_track_comp1[,2],
                              1/(track_comp_1_mean-(2*track_comp_1_sds)))

track_comp_1_mean <- cbind(alb_Ti_track_comp1[,2],1/track_comp_1_mean)

time_track_plus_2Sd <- curve2time(track_comp_1_plus_2Sd,tracked_cycle_period = 173)

time_track_min_2Sd <- curve2time(track_comp_1_min_2Sd,tracked_cycle_period = 173)

time_mean <- curve2time(track_comp_1_mean,tracked_cycle_period = 173)

graphics.off()

Mean_AgeModel_WR<-time_mean
Max_AgeModel_WR<-time_track_plus_2Sd
Min_AgeModel_WR<-time_track_min_2Sd

plot(Mean_AgeModel_WR,type="l", xlab="Depth adjusted (m)", ylab="Time (ka)")
lines(Max_AgeModel_WR,col="green")
lines(Min_AgeModel_WR,col="blue")

legend("bottomright",
       legend = c("Maximum age (+2 SD)","Mean age", "Minimum age (-2 SD)"),
       col = c("green","black","blue"),
       lty = c(1, 1, 1),
       lwd = c(2, 2, 2),
       pch = c(NA, NA, NA),
       pt.cex = 2,
       bty = "n",         # No box around legend
       cex = 1
       )         # Larger legend text

Age_Model_WR <- cbind(Mean_AgeModel_WR,Max_AgeModel_WR[,2],Min_AgeModel_WR[,2])
colnames(Age_Model_WR)<- cbind("Adjusted depth (m)",
                               "Mean age (ka)",
                               "Max age (ka)",
                               "Min age (ka)"
                               )

#Interpolation every 1mm of the Age-depth model to have a direct comparison
#between the age model and the original sampling rate of 1mm 
Mean_Age_Model_WR_1mm <- linterp(Mean_AgeModel_WR, dt=0.001, genplot=T)
Max_Age_Model_WR_1mm <- linterp(Max_AgeModel_WR, dt=0.001, genplot=T)
Min_Age_Model_WR_1mm <- linterp(Min_AgeModel_WR, dt=0.001, genplot=T)

Age_Model_WR_1mm <- cbind(Mean_Age_Model_WR_1mm,
                          Max_Age_Model_WR_1mm[,2],
                          Min_Age_Model_WR_1mm[,2]
                          )

colnames(Age_Model_WR)<- cbind("Adjusted depth (m)",
                               "Mean age (ka)",
                               "Max age (ka)",
                               "Min age (ka)"
                               )
graphics.off()

#Saving of the age model as .csv files
#  write.csv(Age_Model_WR_1mm,"AgeModel_WR_rsp1mm.csv")



#** 11.2.3. Sedimentation rate ----

plot(alb_Ti_track_comp1[,2],
     (100*(alb_Ti_track_comp1[,3]))/173,
     type="l",
     col="black",
     lwd= 2,
     xlab = "Depth adjusted (m)", 
     ylab = "Sed rate (cm/kyr)", 
     ylim = c(0.4,0.55)
     )

lines(alb_Si_track_comp1[,2],
      (100*(alb_Si_track_comp1[,3]))/173,
      col="red",
      lwd=2
      )

lines(alb_Al_track_comp1[,2],
      (100*(alb_Al_track_comp1[,3]))/173,
      col="green",
      lwd=2
      )

lines(alb_K_track_comp1[,2],
      (100*(alb_K_track_comp1[,3]))/173,
      col="blue",
      lwd=2
      )

lines(alb_K_track_comp1[,2],
      (alb_K_track_comp1[,3])-(alb_K_track_comp1[,3])+0.493,
      col="red",
      lwd=2, 
      lty=3
      )

legend("bottomleft",
       legend = c("Ti", "Si", "K", "Al", "ASM = 0.493 cm/kyr (theoretical)" ),
       col = c("black", "red", "blue", "green","red"),
       lty = c(1, 1, 1, 1, 3),
       lwd = c(2, 2, 2, 2, 2),
       pch = c(NA, NA, NA, NA),
       pt.cex = 2,
       bty = "n",         # No box around legend
       cex = 0.75)         # Larger legend text

--------------------------------------------------------
  
#Step 12: Astronomical Time Scale ----

#This step is performed in Microsoft Excel. The ATS is presented in 
#Supplementary Materials 6.

#The built age-model data are set along the composite ACL-ACU core at 
#the corresponding adjusted depth (from 73.452 to 93.609 m).

#The data from 73.452 to 73.72 are overlapping with the one of Zhao et al (2022b).

#From 62.111 m to 73.72 m the absolute age of Zhao et al (2022b) are used. 
#The absolute age of 499.9115 ± 0.9 Ma located at 73.72 m is used to anchor 
#both ACU and ACL in time domain.

#From 73.721 to 93.609 m, the time data used are from our newly built age-model.

--------------------------------------------------------
  
#Step 13: Milankovitch filters of the 1 kyr resampled Miaolingian composite dataset ----

#* 13.1. Creation of a ".csv" file ----

#This step is conducted in Microsoft Excel.

#The created file contains the Adjusted depth, the Time in kyr and the Ti 
#Z-score values as follows:
  
#AdjDepth; Time; Ti

#62.111; 497083.4; 1.57549322
#62.112; 497083.6; 1.10323208
#62.113; 497083.7; 0.63097093
#…; …; …
#93.609; 504154.084; -0.5940352
#93.610; 504154.296; -2.3646753

#I called this file alb_Ti_Mia_WR.csv as this includes all the Ti data for 
#the entire Miaolingian Epoch recorded in the core tuned using the WaverideR 
#age-depth model.



#* 13.2. Reading of the file in R and resampling every 1 kyr of the time interval ----

alb_Ti_Mia <- read.csv("alb_Ti_Mia_WR.csv",sep=";")
alb_Ti_Mia <- cbind(alb_Ti_Mia$Time,alb_Ti_Mia$Ti)
alb_Ti_Mia <- na.omit(alb_Ti_Mia)
alb_Ti_Mia[!is.finite(alb_Ti_Mia)] <- NA
alb_Ti_Mia <- na.omit(alb_Ti_Mia)

MIA<- iso(dat=alb_Ti_Mia, xmin=497000, xmax=505000) 
MIA <- linterp(MIA,dt=1, genplot=T)
MIA_Ti_tuned <- linterp(MIA, genplot=T)



#* 13.3. Resampling every 1 kyr and Milankovitch filters ----

#Resampling every 1 kyr
MIA<- iso(dat=alb_Ti_Mia, xmin=497000, xmax=505000) 
MIA <- linterp(MIA,dt=1, genplot=T)
MIA_Ti_tuned <- linterp(MIA, genplot=T)

#Filter ou the long obliquity cycle
MIA_Ti_1200 <- taner(MIA_Ti_tuned,xmax=0.01, fhigh=1/1150, flow=1/1300, demean = T)

#Filter out the long eccentricity cycle
MIA_Ti_405 <- taner(MIA_Ti_tuned,xmax=0.01, fhigh=1/450, flow=1/350, demean = T)

#Filter out the 173-kyr obliquity cycle  
MIA_Ti_173 <- taner(MIA_Ti_tuned, xmax=0.01, fhigh=1/193, flow=1/150, demean=T) 

#Filter out the short eccentricity cycle
MIA_Ti_100 <- taner(MIA_Ti_tuned, xmax=0.02, fhigh=1/85, flow=1/140, demean = T)

#Filter out the entire obliquity cycle  
MIA_Ti_31_all <- taner(MIA_Ti_tuned, xmax=0.05, fhigh=1/24, flow=1/40, demean=T)

#Filter out the obliquity cycle   that generate the 173 kyr modulation
MIA_Ti_31 <- taner(MIA_Ti_tuned, xmax=0.05, fhigh=1/27, flow=1/40, demean=T)

#Filter out the precession cycle  
MIA_Ti_20 <- taner(MIA_Ti_tuned,xmax=0.07, fhigh=1/15, flow=1/21, demean=T)

#Saving of the filtered Milankovitch periodicities
#write.csv(MIA_Ti_tuned,"MIA_WR_rsp_1kyr.csv") 
#write.csv(MIA_Ti_1200,"MIA_WR_1200_rsp_1kyr.csv")
#write.csv(MIA_Ti_405,"MIA_WR_405_rsp_1kyr.csv")
#write.csv(MIA_Ti_173,"MIA_WR_173_rsp_1kyr.csv")
#write.csv(MIA_Ti_100,"MIA_WR_100_rsp_1kyr.csv")
#write.csv(MIA_Ti_31_all,"MIA_WR_31_all_rsp_1kyr.csv")
#write.csv(MIA_Ti_31,"MIA_WR_31_rsp_1kyr.csv")
#write.csv(MIA_Ti_20,"MIA_WR_20_rsp_1kyr.csv")

--------------------------------------------------------
  
#Step 14: Lag-1 autcorrelation of the Miaolingian Series ----

#* 14.1. Resampling of the dataset every 5 kyr ----

alb_Ti_Mia <- read.csv("alb_Ti_Mia_WR.csv",sep=";")
alb_Ti_Mia <- cbind(alb_Ti_Mia$Time,alb_Ti_Mia$Ti)
alb_Ti_Mia <- na.omit(alb_Ti_Mia)
alb_Ti_Mia[!is.finite(alb_Ti_Mia)] <- NA
alb_Ti_Mia <- na.omit(alb_Ti_Mia)

MIA<- iso(dat=alb_Ti_Mia, xmin=497000, xmax=505000) 
MIA <- linterp(MIA,dt=5, genplot=T)
MIA_Ti_tuned <- linterp(MIA, genplot=T)

#* 14.2. Lag-1 autocorrelation ----

#Lag-1 window set from 20 to 90 kyr to avoid the 5-20 kyr high frequency noise 
#and limit to 90 kyr that correspond to the lower boundary of the 
#short eccentricity band

win_min <- 20
win_max <- 90
n_sim <- 100
dat <- data.frame(MIA_Ti_tuned)
dt2 <- dat[2, 1] - dat[1, 1]
dat <- dat[order(dat[, 1], na.last = NA, decreasing = F), ]
npts <- length(dat[, 1])
start <- dat[1, 1]
end <- dat[length(dat[, 1]), 1]
x1 <- dat[1:(npts - 1), 1]
x2 <- dat[2:(npts), 1]
dx = x2 - x1
dt = mean(dx)
sdt = sd(dx)
xout <- seq(start, end, by = dt)
npts <- length(xout)
interp <- approx(dat[, 1], dat[, 2], xout, method = "linear", n = npts)
d <- as.data.frame(interp)
mat_sim <- matrix(data = NA,
                  nrow = nrow(d),
                  ncol = n_sim)


xout_vals <- d[, 1]
i <- 1
npts <- length(xout_vals)
new_sampling_rate <- NULL

fit <- for(i in 1:n_sim){
  
  if (sdt == 0) {
    new_sampling_rate <- runif(n=1, min = dt, max = win_min)}
  else if ((dt - (2 * sdt)) > 0) {
    new_sampling_rate <- truncnorm::rtruncnorm(
      n = 1,
      a = dt - 2 * sdt,
      b = dt + (2 * sdt),
      mean = dt,
      sd = sdt
    )}else {
      new_sampling_rate <- truncnorm::rtruncnorm(
        n = 1,
        a = dt / 2,
        b = dt + (2 * sdt),
        mean = dt,
        sd = sdt
      )
    }
  
  
  win_size <- stats::runif(n = 1,
                           min = min(c(win_min, win_max)),
                           max = max(c(win_min, win_max)))
  
  dt_new <- new_sampling_rate
  
  xout_new <- seq(from = start, to = end, by = dt_new)
  
  npts_new <- length(xout_new)
  
  interp_new <- approx(dat[, 1], 
                       dat[, 2], 
                       xout_new,
                       method = "linear", 
                       n = npts
  )
  
  d_new <- as.data.frame(interp_new)
  d_new[, 3] <- NA
  for (k in 1:nrow(d_new)) {
    row_nr_1 <- DescTools::Closest(d_new[, 1], 
                                   d_new[k, 1] - (win_size / 2),
                                   which = TRUE)
    row_nr_2 <- DescTools::Closest(d_new[, 1], 
                                   d_new[k, 1] + (win_size / 2), 
                                   which = TRUE)
    data_sel <- d_new[row_nr_1[1]:row_nr_2[1], ]
    corr <- acf(data_sel[, 2], plot = F)
    a <- as.numeric(unlist(corr[1])[1])
    d_new[k, 3] <- a
  }
  yleft_comp <- d_new[1, 3]
  yright_com <- d_new[nrow(d_new), 3]
  app <- approx(
    d_new[, 1],
    d_new[, 3],
    xout_vals,
    method = "linear",
    n = npts,
    yleft = yleft_comp,
    yright = yright_com
  )
  app_res <- cbind(app$y)
  app_res_norm <- (app_res - min(app_res, na.rm = TRUE)) /
    (max(app_res, na.rm = TRUE) - min(app_res, na.rm = TRUE))
  mat_sim[,i] <- app_res_norm
}

mat_sim_mean <- rowMeans(mat_sim)
mat_sim_sd <- rowSds(mat_sim)
results <- cbind(xout_vals, mat_sim_mean, mat_sim_sd)

lag_1 <- results
graphics.off()

plot(lag_1, 
     type = "l", 
     main = "Lag 1 over time",
     xlab = "Age Before Present (kyr)", 
     ylab = "Lag 1 power"
     )

#    write.csv(lag_1,"lag1-rsp5kyr-window20-90kyr.csv") 



#* 14.3. Filtering of the Milankovitch cycles and comparison with lag-1 results ----

MIA_Ti__405 <- taner(MIA_Ti_tuned[,c(1,2)],flow=1/350,fhigh=1/450,roll=10^20,xmax=1/50)
MIA_Ti__173<- taner(MIA_Ti_tuned[,c(1,2)],flow=1/150,fhigh=1/193,roll=10^20,xmax=1/50)
MIA_Ti__100<- taner(MIA_Ti_tuned[,c(1,2)],flow=1/140,fhigh=1/85,roll=10^20,xmax=1/50)
MIA_Ti__31 <- taner(MIA_Ti_tuned[,c(1,2)],flow=1/27,fhigh=1/40,roll=10^20,xmax=1/50)

graphics.off()

lag_1_background<- taner(lag_1[,c(1,2)],flow=1/500,fhigh=0,roll=10^20,xmax=1/50)

MIA_Ti__405_hilbert <- hilbert(MIA_Ti__405)
MIA_Ti__173_hilbert <- hilbert(MIA_Ti__173)
MIA_Ti__100_hilbert <- hilbert(MIA_Ti__100)
MIA_Ti__31_hilbert <- hilbert(MIA_Ti__31)

graphics.off()

#Creation of a matrix to display the lag 1 results and the Milankovitch filters
#in one single figure    
layout.matrix <- matrix(c(1,2,3,4,5), nrow = 5, ncol = 1)
graphics::layout(mat = layout.matrix,
                 heights = c(1), widths = c(1)) 

#Lag 1 results
par(mar = c(0,3,2,1)) 

plot(lag_1[,c(1,2)],
     type="l",
     ylab="", 
     xlab = "", 
     ylim=c(0.3,0.9), 
     bty = "n", 
     yaxt = "n", 
     xaxt = "n"
     )

axis(2, seq(from = 0.3, to = 0.9, by = 0.2), cex.axis = 0.5, las = 1)
title(ylab="Lag 1 power", line=2.2, cex.lab=0.75)
title(main ="Lag 1 and sea level trend", line=0, cex.main=0.75)
lines(lag_1_background[,1],(lag_1_background[,2]),col="blue", lwd=1.5)


# 405 kyr cycle and envelope
par(mar = c(0,3,0,1))

plot(MIA_Ti__405_hilbert[,1],
     (MIA_Ti__405_hilbert[,2]),
     type = "l",
     col="red", 
     lwd=1.5, 
     ylim = c(-0.6,0.6), 
     xlab = "", 
     ylab = "",
     xaxt ="n", 
     yaxt = "n",
     bty = "n"
     )

lines(MIA_Ti__405[,1],(MIA_Ti__405[,2]),col="grey", lwd=0.5)
axis(2, seq(from = -0.3, to = 0.3, by = 0.3), cex.axis = 0.5, las=1)
title(ylab="filter output", line=2.2, cex.lab=0.75)
title(main ="405 kyr cycle and envelope", line=-2.5, cex.main=0.75)

#173 kyr cycle and envelope
par(mar = c(0,3,0,1))

plot(MIA_Ti__173_hilbert[,1],
     (MIA_Ti__173_hilbert[,2]),
     col="purple",
     lwd=1.5, 
     type = "l",
     ylab = "", 
     xlab = "",
     yaxt ="n",
     xaxt = "n", 
     bty = "n", 
     ylim = c(-0.6,0.7)
     )

lines(MIA_Ti__173[,1],
      (MIA_Ti__173[,2]-mean(MIA_Ti__173[,2])),
      col="grey",
      lwd=0.5
      )

axis(2, seq(from = -0.6, to = 0.6, by = 0.3), cex.axis = 0.5, las=1)
title(ylab="filter output", line=2.2, cex.lab=0.75)
title(main ="173 kyr cycle and envelope", line=-2, cex.main=0.75)

#100 kyr cycle and envelope
par(mar = c(0,3,0,1))

plot(MIA_Ti__100_hilbert[,1],
     (MIA_Ti__100_hilbert[,2]),
     col="green", 
     lwd=1.5, 
     type = "l", 
     ylim = c(-1.1,1.1), 
     bty = "n", 
     ylab ="", 
     xlab = "", 
     yaxt = "n", 
     xaxt = "n"
     )

lines(MIA_Ti__100[,1],(MIA_Ti__100[,2]),col="grey", lwd=0.5)
axis(2, seq(from = -0.9, to = 0.9, by = 0.45), cex.axis = 0.5, las=1)
title(ylab="filter output", line=2.2, cex.lab=0.75)
title(main ="100 kyr cycle and envelope", line=-1, cex.main=0.75)

#31 kyr cycle and envelope
par(mar = c(3,3,0,1))

plot(MIA_Ti__31_hilbert[,1],
     (MIA_Ti__31_hilbert[,2]),
     col="orange", 
     lwd=1.5, 
     type = "l", 
     ylim = c(-1.6,1.6),
     bty = "n",
     ylab ="", 
     xlab ="", 
     yaxt = "n", 
     xaxt = "n"
     )

lines(MIA_Ti__31[,1],(MIA_Ti__31[,2]),col="grey", lwd=0.5)
axis(2, seq(from = -1.5, to = 1.5, by = 0.75), cex.axis = 0.5, las=1)
title(ylab="filter output", line=2.2, cex.lab=0.75)
axis(1, seq(from = 497000, to = 504000, by = 1000), cex.axis = 0.5)
title(xlab="Age Before Present (kyr)", line=2, cex.lab=0.75)
title(main ="31 kyr cycle and envelope", line=-0.5, cex.main=0.75)

#Saving the lag-1 results and Milankovitch filters
#    write.csv(lag_1_background,"Lag1_20-90window_background_0-1over500kyr.csv")
#    write.csv(MIA_Ti__405_hilbert,"MIA_TI_rsp_5kyr_405kyr_hilbert.csv")
#    write.csv(MIA_Ti__173_hilbert,"MIA_TI_rsp_5kyr_173kyr_hilbert.csv")
#    write.csv(MIA_Ti__100_hilbert,"MIA_TI_rsp_5kyr_100kyr_hilbert.csv")
#    write.csv(MIA_Ti__31_hilbert,"MIA_TI_rsp_5kyr_31kyr_hilbert.csv")

