remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data ####

source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

#Seleccionamos la ventana temporal. Enero de 2005 a Diciembre de 2019.
Yl <- window(Yl.f, start = c(2005, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2005, 01), end = c(2019, 12))


# Inciso 1 ####

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
#MARTIN: No entiendo muy bien por qué, pero solo eliminamos 10 valores, no 12. 
#MARTIN: ¿No deberíamos comernos únicamente los primeros 3 valores?
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation
VAR <- VAR(Ydt, p = p, type = "const")

m <- VAR$K # No. of variables in the VAR
N <- VAR$obs # No. of effective sample observations, excluding "p" starting values

# Ad hoc function
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}

# Re-estimate VAR (no feedback from local vars. to pcom)
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
VAR


# Model checking
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation

# A Matrix
#esta representa a la matriz A0^(-1)
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# SVAR estimation (AB model configuration)
#aca simplemente estimamos el SVAR imponiendo ciertas assumptions embebidas en A y en B
#los NA en A y en B se van a estimar como en la slide 24 de la lecture 4
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR


H <- 18 # Horizon
H_ERPT <- 120 # Horizon for ERPT

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")


#Reporte de resultados inciso 1
#(quiza podriamos hacer graficos mas lindos con ggplot)

# IRF
SVAR.SIRF <- SVAR.sirf(SVAR, H)
plot.sirf(SVAR.SIRF, m, H)

# Cumulative IRF
SVAR.SIRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
plot.sirf(SVAR.SIRF.c, m, H)

# FEVD
SVAR.FEVD <- SVAR.fevd(SVAR, H)
plot.fevd(SVAR.FEVD, m, H)

# HD
SVAR.HD <- SVAR.hd(SVAR)
plot.hd(Yd, SVAR.HD, m, pmax)

# ERPT
SVAR.ERPT <- SVAR.erpt(SVAR, H_ERPT, 3, 2, cumulative = TRUE)
plot.erpt(SVAR.ERPT, H_ERPT)

#PENDIENTE: Revisar cuestiones de estacionariedad de las series.


# Inciso 2 ####
#Estimamos el modelo sin el índice de precios internos.

#Duda: ¿hay que hacer esto de nuevo o con los resultados que obtuvimos con las 3 variables alcanza?
Yd_2 <- Yd[, 2:3]
popt_2 <- VARselect(Yd_2, lag.max = pmax, type = "const")
popt_2
p_2 <- popt_2$selection[1] # AIC

Yd0_2 <- Yd_2[1:pmax, ] # Initial values
Ydt_2 <- Yd_2[(pmax - p_2 + 1):nrow(Yd_2), ] # Starting in Jan-04

# Estimation
VAR_2 <- VAR(Ydt_2, p = p_2, type = "const")

m_2 <- VAR_2$K # No. of variables in the VAR
N_2 <- VAR_2$obs

roots(VAR_2, modulus = TRUE)
serial.test(VAR_2, lags.bg = 12, type = "ES")

# SVAR estimation

# A Matrix
Amat_2 <- diag(m_2)
for (i in 2:m_2) {
  for (j in 1:(i - 1)) {
    Amat_2[i, j] <- NA
  }
}

# B Matrix
Bmat_2 <- matrix(0, m_2, m_2)
for (i in 1:m_2) {
  Bmat_2[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_2 <- SVAR(VAR_2, Amat = Amat_2, Bmat = Bmat_2, lrtest = FALSE)
SVAR_2

#Reportamos resultados modelo inciso 2.

#hay que volver a llamarla m porque las funciones usadas a continuación están definidas usando m, y la cantidad de variables no se puede especificar de otra manera
#después de haber calculado las cosas que queremos volvemos a darle el valor inicial.
m <- m_2

# IRF
SVAR_2.SIRF <- SVAR.sirf(SVAR_2, H)
plot.sirf(SVAR_2.SIRF, m, H)

# Cumulative IRF
SVAR_2.SIRF.c <- SVAR.sirf(SVAR_2, H, cumulative = TRUE)
plot.sirf(SVAR_2.SIRF.c, m, H)

# FEVD
SVAR_2.FEVD <- SVAR.fevd(SVAR_2, H)
plot.fevd(SVAR_2.FEVD, m, H)

# HD
SVAR_2.HD <- SVAR.hd(SVAR_2)
plot.hd(Yd_2, SVAR_2.HD, m, pmax)

# ERPT
SVAR_2.ERPT <- SVAR.erpt(SVAR_2, H_ERPT, 2, 1, cumulative = TRUE)
plot.erpt(SVAR_2.ERPT, H_ERPT)

#Pendiente: comparar con resultados modelo inciso 1.

######## (3) ### run each subsample and then the same with 3 lags

######## from January 2005 to July 2011

remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("PS1_Data.R")
library(vars)
Yl.f <- cbind(pcom, er, pc)
Yl.f <- log(Yl.f) 
Yd.f <- 100 * diff(Yl.f)

Yl <- window(Yl.f, start = c(2003, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2003, 01), end = c(2019, 12))

Yone <- window(Yd.f, start = c(2004, 01), end = c(2011, 07))

pmax <- 12

popt <- VARselect(Yone, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]

Yone0 <- Yone[1:pmax, ] # Initial values, 12M
Yonet <- Yone[(pmax - p + 1):nrow(Yone), ]

VAR <- VAR(Yonet, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Yone0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

######## from August 2011 to September 2015

remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc)
Yl.f <- log(Yl.f) 
Yd.f <- 100 * diff(Yl.f)

Yl <- window(Yl.f, start = c(2003, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2003, 01), end = c(2019, 12))

pmax <- 12

Ytwo <- window(Yd.f, start = c(2010, 08), end = c(2015, 09))

popt <- VARselect(Ytwo, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]

Ytwo0 <- Ytwo[1:pmax, ] # Initial values, 12M
Ytwot <- Ytwo[(pmax - p + 1):nrow(Ytwo), ]

VAR <- VAR(Ytwot, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Ytwo0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

######## from October 2015 to August 2019

remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc)
Yl.f <- log(Yl.f) 
Yd.f <- 100 * diff(Yl.f)

Yl <- window(Yl.f, start = c(2003, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2003, 01), end = c(2019, 12))

pmax <- 12

Ythree <- window(Yd.f, start = c(2014, 10), end = c(2019, 08))

popt <- VARselect(Ythree, lag.max = pmax, type = "const")
popt
p <- popt$selection[2] # [check and change for each sample} /// [change again if autocorrelated errors]

Ythree0 <- Ythree[1:pmax, ] # Initial values, 12M
Ythreet <- Ythree[(pmax - p + 1):nrow(Ythree), ]

VAR <- VAR(Ythreet, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Ythree0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# Punto 7 - M2 private
###
library(ggplot2)
df <- data.frame(as.matrix(Yd), date=time(Yd))
plot.ts(diff(log(M2private)))
plot.ts(diff(log(M2adj)))
###

remove(list = ls(all.names = TRUE))
gc()
library(vars)
library(seasonal)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

M2private <- read_excel("m2.xlsx")
M2private <- ts(M2private, start = c(2003, 12), frequency = 12)
seas.adj <- seas(M2private)
M2adj <- seas.adj$series$s11
source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc, M2adj)
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences
Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

pmax <- 12
popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]
p

Yd0 <- Yd[1:pmax, ] # Initial values, 12M
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # [check as it depends upon the lags]


VAR <- VAR(Ydt, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM
pars.S

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# Punto 8 - Time deposits
remove(list = ls(all.names = TRUE))
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pf <- read_excel("pftna.xlsx")
pf <- ts(pf, start = c(2003, 01), frequency = 12)
source("PS1_Data.R")

pf <- pf/(100*12)
pf <- log(1+pf)
pf.d <- diff(pf)

#plot(pf)

Yl.f <- cbind(pcom, er, pc)
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl.f <- cbind(Yl.f, pf)
Yd.f <- cbind(Yd.f, pf.d)

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

pmax <- 12
popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]
p

Yd0 <- Yd[1:pmax, ] # Initial values, 12M
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # [check as it depends upon the lags]


VAR <- VAR(Ydt, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM
pars.S

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# Punto 9 - M2, rates, activity, wages ####

remove(list = ls(all.names = TRUE))
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("PS1_Data.R")

wages.file <- paste(tempfile(), ".csv", sep = "")
download.file("https://infra.datos.gob.ar/catalog/sspm/dataset/157/distribution/157.2/download/remuneracion-bruta-promedio-asalariados-registrados-sector-privado-segun-actividad-sin-estacionalidad.csv
", wages.file, mode = "wb")

wages <- read.csv(wages.file)
wages <- wages[c(1,16)]
wages_num <- as.numeric(wages$rem_bruta_asal_sin_est_total)
remove(wages.file)

wages_num <- wages_num[complete.cases(wages_num)] # delete NAs
wages_num <- ts(wages_num, start = c(1995, 01), end = c(2020, 03), frequency = 12)
wages_num  <- window(wages_num, start = c(2003, 12))

wages_real <- wages_num * tail(pc, n=1)/pc
wages_real <- wages_real / tail(wages_real, n=1)
wages_real  <- window(wages_real, start = c(2003, 12))

M2private <- read_excel("m2.xlsx")
M2private <- ts(M2private, start = c(2003, 12), frequency = 12)
seas.adj <- seas(M2private)
M2adj <- seas.adj$series$s11

pf <- read_excel("pftna.xlsx")
pf <- ts(pf, start = c(2003, 01), frequency = 12)
pf <- window(pf, start = c(2003, 12), end = c(2019, 12))
pf <- pf/(100*12)
pf <- log(1+pf)
pf.d <- diff(pf)

emae.file <- paste(tempfile(), ".csv", sep = "")
download.file("https://infra.datos.gob.ar/catalog/sspm/dataset/143/distribution/143.3/download/emae-valores-anuales-indice-base-2004-mensual.csv", emae.file, mode = "wb")
emae <- read.csv(emae.file)
emae <- emae[c(1,6)]
emae_num <- as.numeric(emae$emae_desestacionalizada_vm)
emae_num <- emae_num[complete.cases(emae_num)] # delete NAs
emae_num <- ts(emae_num, start = c(2004, 01), end = c(2020, 07), frequency = 12)
emae_num <- window(emae_num, start = c(2003, 12), end = c(2019, 12))*100

library(tstools)

anio <-2009

d2009_1<- create_dummy_ts(end_basic = c(anio,06), dummy_start = c(anio,03), dummy_end = c(anio,06), sp = TRUE, start_basic = c(2004,01),  dummy_value = -1, frequency = 12)
d2009_2<- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(anio,08), dummy_end = c(anio,8), sp = TRUE, start_basic = c(anio,7),  dummy_value = 1, frequency = 12)
d2009_2[1] = 1
d2009 <- ts(c(d2009_1,d2009_2),  start = start(d2009_2),frequency = frequency(121))
d2009 <- ts(d2009, start = c(2004,01), end = c(2019, 12), frequency = 12)

anio <-2012

d2012_1<- create_dummy_ts(end_basic = c(anio,06), dummy_start = c(anio,03), dummy_end = c(anio,06), sp = TRUE, start_basic = c(2004,01),  dummy_value = -1, frequency = 12)
d2012_2<- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(anio,08), dummy_end = c(anio,8), sp = TRUE, start_basic = c(anio,7),  dummy_value = 1, frequency = 12)
d2012_2[1] = 1
d2012 <- ts(c(d2012_1,d2012_2),  start = start(d2012_2),frequency = frequency(121))
d2012 <- ts(c(d2012_1,d2012_2),  start = start(d2012_2),frequency = frequency(121))
d2012 <- ts(d2012, start = c(2004,01), end = c(2019, 12), frequency = 12)

anio <-2018

d2018_1<- create_dummy_ts(end_basic = c(anio,06), dummy_start = c(anio,03), dummy_end = c(anio,06), sp = TRUE, start_basic = c(2004,01),  dummy_value = -1, frequency = 12)
d2018_2<- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(anio,08), dummy_end = c(anio,8), sp = TRUE, start_basic = c(anio,7),  dummy_value = 1, frequency = 12)
d2018_2[1] = 1
d2018 <- ts(c(d2018_1,d2018_2),  start = start(d2018_2),frequency = frequency(121))
d2018 <- ts(d2018, start = c(2004,01), end = c(2019, 12), frequency = 12)

dlog.emae.reg <- lm(emae_num ~ -1 + d2009 + d2012 + d2018)
dlog.emae.adj <- dlog.emae.reg$residuals

dlog_residuos <- as.numeric(dlog.emae.adj)
dlog_residuos <- ts(dlog_residuos, start = c(2004, 01), end = c(2019, 12), frequency = 12)


Yl.f <- cbind(pcom, er, pc, wages_real, M2adj)
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences
Yl.f <- window(Yl.f, start = c(2003, 12), end = c(2019, 12))
Yd.f <- window(Yd.f, start = c(2003, 12), end = c(2019, 12))

Yl.f <- cbind(Yl.f, pf, dlog_residuos)
Yd.f <- cbind(Yd.f, pf.d, dlog_residuos)

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

pmax <- 12
popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]
p

Yd0 <- Yd[1:pmax, ] # Initial values, 12M
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # [check as it depends upon the lags]

VAR <- VAR(Ydt, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM
pars.S

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# Punto 10 ####
#Forma 1: fuente oficial

#Descargamos los valores de la fuente, tomando el promedio mensual:
dolar_ccl <- read.csv(url("https://apis.datos.gob.ar/series/api/series/?collapse=month&collapse_aggregation=avg&ids=168.1_T_CAMBIDRS_D_0_0_29&limit=5000&format=csv"))

#Construimos la serie y la restringimos
dolar_ccl <- ts(dolar_ccl$tipo_cambio_implicito_en_adrs, start = c(2002, 04), frequency = 12)
dolar_ccl <- window(dolar_ccl, start = c(2004, 01), end = c(2019, 12))

#Graficamos para verificar
plot(dolar_ccl, type = "l",lwd=2, col="black", xlab="",ylab="",bty="n", main="CCL implícito en ADRs", ylim=c(0,80))


#Forma 2: tomando valores de la acción del Banco Galicia
#Para obtener tipo de cambio implícito en los ADR (Contado Contra Liquidación)
#es necesario hacer la siguiente cuenta:
#Dólar CCL = (Precio Local de Acción/Precio del ADR)* Factor de conversión
#En el caso del Grupo Financiero Galicia S.A. (Banco Galicia), el factor correspondiente es 10.


#Utilizamos una descarga automática basada en Yahoo Finance de los valores en NY y en BA, tomando precio de cierre mensual ajustado:
library(tseries)
P_adr <- get.hist.quote(instrument = "GGAL", start = "2004-01-01", end = "2019-12-31", 
                    quote = "AdjClose", provider = "yahoo", compression = "m", 
                    retclass = "zoo")
plot(P_adr,type = "l",lwd=2, col="black", xlab="",ylab="",bty="n", main="GGAL Adjusted Price", ylim=c(0,70))

P_local <- get.hist.quote(instrument = "GGAL.BA", start = "2004-01-01", end = "2019-12-31", 
                          quote = "AdjClose", provider = "yahoo", compression = "m", 
                          retclass = "zoo")
plot(P_local,type = "l",lwd=2, col="black", xlab="",ylab="",bty="n", main="GGAL.BA Adjusted Price", ylim=c(0,170))

#Construimos la serie
dolar_implicito <- (P_local/P_adr)*10
dolar_ccl2 = ts(dolar_implicito, start = c(2004, 01), frequency = 12)

#Graficamos para verficar
plot(dolar_ccl2, type = "l",lwd=2, col="black", xlab="",ylab="",bty="n", main="CCL a partir de $GGAL", ylim=c(0,140))


#Analizamos y graficamos las dos alternativas para comparar:
discrepancia = (dolar_ccl-dolar_ccl2)/dolar_ccl
mean(discrepancia)
sd(discrepancia)
ts.plot(discrepancia)
#Conclusión: aunque en promedio no son muy distintas, las medidas difieren bastante en algunos meses. Parece ser mejor usar el oficial, ante la duda.



# Comenzamos el análisis VAR
# Definimos las nuevas variables, con ventana temporal enero2005-diciembre2019
Yl.f_10 <- cbind(pcom, dolar_ccl, pc) # Variables en log
Yl.f_10 <- log(Yl.f_10) # log transformation
Yd.f_10 <- 100 * diff(Yl.f_10) # Variables en log-differences
Yl_10 <- window(Yl.f_10, start = c(2005, 01), end = c(2019, 12))
Yd_10 <- window(Yd.f_10, start = c(2005, 01), end = c(2019, 12))


# Volvemos a llamar al paquete, para evitar <<enmascaramiento>> del paquete ts:
library(vars)

# Lag order selection

SVAR_10popt_10 <- VARselect(Yd_10, lag.max = pmax, type = "const")
popt_10
p_10 <- popt_10$selection[1] # AIC

# Valores iniciales
Yd0_10 <- Yd_10[1:pmax, ] # Initial values
Ydt_10 <- Yd_10[(pmax - p_10 + 1):nrow(Yd_10), ] 

# Estimation
VAR_10 <- VAR(Ydt_10, p = p_10, type = "const")

# Control
m_10 <- VAR_10$K # No. of variables in the VAR
N_10 <- VAR_10$obs
roots(VAR_10, modulus = TRUE)
serial.test(VAR_10, lags.bg = 2, type = "ES") #DUDA: me lo estima con un solo lag, ¿les parece OK hacerlo así?

# Re-estimación con restricciones:
# Re-estimate VAR (no feedback from local vars. to pcom)
VAR_10 <- restrict(VAR_10, method = "man", resmat = matC(m_10, p_10, 1))
VAR_10

# SVAR estimation

# A Matrix
Amat_10 <- diag(m_10)
for (i in 2:m_10) {
  for (j in 1:(i - 1)) {
    Amat_10[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat_10 <- matrix(0, m_10, m_10)
for (i in 1:m_10) {
  Bmat_10[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_10 <- SVAR(VAR_10, Amat = Amat_10, Bmat = Bmat_10, lrtest = FALSE)

# Reportamos los resultados del modelo
m <- m_10
# IRF
SVAR_10.SIRF <- SVAR.sirf(SVAR_10, H)
plot.sirf(SVAR_10.SIRF, m=m_10, H)

# Cumulative IRF
SVAR_10.SIRF.c <- SVAR.sirf(SVAR_10, H, cumulative = TRUE)
plot.sirf(SVAR_10.SIRF.c, m = m_10, H)

# FEVD
SVAR_10.FEVD <- SVAR.fevd(SVAR_10, H)
plot.fevd(SVAR_10.FEVD, m = m_10, H)

# HD 
N <- N_10 #Hay que redifinir para que funcione la función SVAR.hd
SVAR_10.HD <- SVAR.hd(SVAR_10)
plot.hd(Yd_10, SVAR_10.HD, m = m_10, pmax)

# ERPT
SVAR_10.ERPT <- SVAR.erpt(SVAR_10, H_ERPT, 3, 2, cumulative = TRUE)
plot.erpt(SVAR_10.ERPT, H_ERPT)

# Replicación con bootstrap y gráfico final con bandas de confianza
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications
Yb_10 <- boot.rb.replicate(VAR_10, Yd0_10, pmax, R)
SVAR.SIRF.boot_10 <- SVAR.sirf.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot_10, m = m_10, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot_10 <- SVAR.sirf.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_10, m = m_10 , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_10 <- SVAR.fevd.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_10, m = m_10, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_10 <- SVAR.erpt.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot_10, H_ERPT)


# Punto 11: primer ordenamiento ####
#Construimos la brecha  #Esto está por si quieren ver
#er_para_brecha = window(er, start = c(2004, 01), end = c(2019, 12))
#brecha <- (dolar_ccl - er_para_brecha)
#plot(brecha) 


# Comenzamos el análisis VAR
# Definimos las nuevas variables, con ventana temporal enero2005-diciembre2019
pcom_log <- log(pcom)
er_log <- log(er)
brecha_log <- log(dolar_ccl)-log(er)
pc_log <- log(pc)
Yl.f_11 <- cbind(pcom_log, er_log, brecha_log, pc_log) 
Yd.f_11 <- 100 * diff(Yl.f_11) # Variables en log-differences
Yl_11 <- window(Yl.f_11, start = c(2005, 01), end = c(2019, 12))
Yd_11 <- window(Yd.f_11, start = c(2005, 01), end = c(2019, 12))

# Comenzamos análisis VAR
popt_11 <- VARselect(Yd_11, lag.max = pmax, type = "const")
popt_11
p_11 <- popt_11$selection[1] # AIC

# Valores iniciales
Yd0_11 <- Yd_11[1:pmax, ] # Initial values
Ydt_11 <- Yd_11[(pmax - p_11 + 1):nrow(Yd_11), ] 

# Estimation
VAR_11 <- VAR(Ydt_11, p = p_11, type = "const")

# Control
m_11 <- VAR_11$K # No. of variables in the VAR
N_11 <- VAR_11$obs
roots(VAR_11, modulus = TRUE)
serial.test(VAR_11, lags.bg = 1, type = "ES") #DUDA: me lo estima con un solo lag, ¿les parece OK hacerlo así?

# Re-estimación con restricciones:
# Re-estimate VAR (no feedback from local vars. to pcom)
VAR_11 <- restrict(VAR_11, method = "man", resmat = matC(m_11, p_11, 1))
VAR_11

# SVAR estimation

# A Matrix
Amat_11 <- diag(m_11)
for (i in 2:m_11) {
  for (j in 1:(i - 1)) {
    Amat_11[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat_11 <- matrix(0, m_11, m_11)
for (i in 1:m_11) {
  Bmat_11[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_11 <- SVAR(VAR_11, Amat = Amat_11, Bmat = Bmat_11, lrtest = FALSE)


# Replicación con bootstrap y gráfico final con bandas de confianza
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications

Yb_11 <- boot.rb.replicate(VAR_11, Yd0_11, pmax, R)
N <- N_11
SVAR.SIRF.boot_11 <- SVAR.sirf.boot(SVAR_11, Amat_11, Bmat_11, Yb_11, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot_11, m = m_11, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot_11 <- SVAR.sirf.boot(SVAR_11, Amat_11, Bmat_11, Yb_11, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_11, m = m_11 , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_11 <- SVAR.fevd.boot(SVAR_11, Amat_11, Bmat_11, Yb_11, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_11, m = m_11, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_11 <- SVAR.erpt.boot(SVAR_11, Amat_11, Bmat_11, Yb_11, pmax, H_ERPT, 4, 2, a, R, cumulative = TRUE) # DUDA: Acá puse el 4 en vez del 3 porque sino había un problema con las dimensiones y no estimaba, pero NO ESTOY muy seguro. No me termina de quedar claro qué representan estos argumentos (los números) en la función.
plot.erpt.boot(SVAR.ERPT.boot_11, H_ERPT)


#Punto 12 ####
#Tomamos como períodos con controles de capitales al período octubre 2011-diciembre2015 y desde septiembre 2019.
#Fuentes: AFIP con la Resolución General 3210 y la Resolución General 3819 , Poder Ejecutivo Nacional con DNU 19/609

#Creamos dos variables dummies:
library(tstools) 

dummy_cepo1 <- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(2011,10), dummy_end =c(2015,12), sp= NULL, start_basic = c(2004, 01), frequency = 12)
dummy_cepo2 <- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(2019,09), dummy_end =c(2019,12), sp= NULL, start_basic = c(2004, 01), frequency = 12)

#Creamos la variable de la brecha que toma valor 0 cuando no hay controles de capitales:
brecha_con_cepo1 <- brecha_log*dummy_cepo1
brecha_con_cepo1 <- window(brecha_con_cepo1, end = c(2019, 08))
brecha_con_cepo2 <- brecha_log*dummy_cepo2
brecha_con_cepo2 <- window(brecha_con_cepo2, start = c(2019, 09))

brecha_con_cepo_log <- concat_ts(brecha_con_cepo1, brecha_con_cepo2)
plot(brecha_con_cepo_log)

#Eliminamos variables intermedias:
remove(brecha_con_cepo1, brecha_con_cepo2)

