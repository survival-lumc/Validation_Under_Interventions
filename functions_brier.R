#----------------------
#function for Brier score

Brier <- function(time, status, risk, seq.time, weights=rep(1, length(time)))
{
  uncensored <- time>seq.time | (status>0 & !is.na(risk))
  status[time>seq.time] <- 0  # only evaluate event status up to seq.time
  return(1/(length(time))*sum((risk[uncensored]-status[uncensored])^2*weights[uncensored]))
}

#----------------------
#function for scaled Brier score (IPA)

ipa <- function(time, status, risk, seq.time, weights=rep(1, length(time)))
{
  #sf <- survfit(Surv(time, status)~1,weights=weights)
  #nullrisk <- 1-min(sf$surv[sf$time<=seq.time])
  uncensored <- time>seq.time | (status>0 & !is.na(risk))
  nullrisk <- 1/(length(time))*sum(status[uncensored]*weights[uncensored])
  brier1 <- Brier(time, status, risk, seq.time, weights=weights)
  brier0 <- Brier(time, status, rep(nullrisk,length(time)), seq.time, weights=weights)  
  return((brier0-brier1)/brier0*100)
}


