library(lubridate)
library(dplyr)
library(tidyr)
library(car)

readin <- function(dyear) {
  
  lines <- readLines(paste("bigxii", as.character(dyear), ".csv", sep = ""), warn = FALSE)
  
  starts <- which(grepl("#", lines))
  ends <- c((starts[2:length(starts)]-1), length(lines))
  
  new_frames <- mapply(function(start, end) {
    read.csv(text=paste0(lines[start:end], collapse="\n"), header=TRUE)
  }, starts, ends, SIMPLIFY=FALSE)
  
  time <- ms(as.character(substr(new_frames[[12]][,4],1,5)))
  new_frames[[12]]$top <- minute(time) + second(time)/60
  
  new_frames <- lapply(new_frames, function(x) janitor::remove_empty(x, which = "cols"))
  
  vars <- c("ScorOff", "ScorDef", "RushOff", "RushDef", "PassOff",
            "PassDef", "FirstDown", "OppFirstDown", "3rdDownPct", 
            "Opp3rdDownPct", "PenYds", "TimeOfPoss", "TOMargin")
  
  team <- character()
  year <- numeric()
  stat <- character()
  value <- numeric()
  for (i in 1:length(vars)) {
    tteam <- as.character(new_frames[[i]][,1]); team <- c(team, tteam)
    tyear <- rep.int(dyear, 10); year <- c(year, tyear)
    tstat <- rep.int(vars[i], 10); stat <- c(stat, tstat)
    tvalue <- new_frames[[i]][,ncol(new_frames[[i]])]; value <- c(value, tvalue)
  }
  
  bigxii <- data.frame(team, year, stat, value)
  return(bigxii)
  
}

bigxiilist <- lapply(2012:2019, readin)
bigxii_long <- bind_rows(bigxiilist)

bigxii <- bigxii_long %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(WinMargin = ScorOff - ScorDef) %>%
  rename(ThirdDownPct = `3rdDownPct`) %>%
  select(team, year, WinMargin, everything(), -ScorOff, -ScorDef)

mod1 <- lm(WinMargin ~ RushOff + RushDef + PassOff + PassDef + FirstDown + OppFirstDown +
             ThirdDownPct + Opp3rdDownPct + PenYds + TimeOfPoss + TOMargin, data = bigxii)

par(mfrow=c(2,2))
plot(mod1)
par(mfrow=c(1,1))

mmps(mod1)

summary(mod1)
vif(mod1)

X <- as.matrix(bigxii[4:14])
tranx <- powerTransform(X, family = "yjPower") 
summary(tranx)
testTransform(tranx, lambda = rep.int(1,11))

bigxii$tRushOff = sqrt(bigxii$RushOff)
bigxii$lRushDef = log(bigxii$RushDef)
bigxii$lPassDef = log(bigxii$PassDef)
bigxii$tOppFirstDown = sqrt(bigxii$OppFirstDown)
bigxii$tThirdDownPct = sqrt(bigxii$ThirdDownPct)
bigxii$tPenYds = sqrt(bigxii$PenYds)

temp_mod <- lm(WinMargin ~ tRushOff + lRushDef + PassOff + lPassDef + FirstDown + tOppFirstDown +
                 tThirdDownPct + Opp3rdDownPct + tPenYds + TimeOfPoss + TOMargin, data = bigxii)
tranmod <- powerTransform(temp_mod, family="yjPower")
summary(tranmod)

mod2 <- lm(WinMargin ~ tRushOff + lRushDef + PassOff + lPassDef + FirstDown + tOppFirstDown +
             tThirdDownPct + Opp3rdDownPct + tPenYds + TimeOfPoss + TOMargin, data = bigxii)

par(mfrow=c(2,2))
plot(mod2)
par(mfrow=c(1,1))

cd <- cooks.distance(mod2)
par(mfrow = c(1,1), cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.9)
plot(bigxii$tRushOff, cd, xlab="", ylab="", 
     main = "")
abline(h=(4/80), lty=2)

rbigxii <- bigxii[-which(cd>(4/80)),]

mod2_2 <- lm(WinMargin ~ tRushOff + lRushDef + PassOff + lPassDef + FirstDown + tOppFirstDown +
             tThirdDownPct + Opp3rdDownPct + tPenYds + TimeOfPoss + TOMargin, 
             data = rbigxii)

cd2 <- cooks.distance(mod2_2)
par(mfrow = c(1,1), cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.9)
plot(rbigxii$tRushOff, cd2, xlab="", ylab="", 
     main = "")
abline(h=(4/73), lty=2)

par(mfrow=c(2,2))
plot(mod2_2)
par(mfrow=c(1,1))

mmps(mod2_2)

summary(mod2_2)
vif(mod2_2)

step.model <- MASS::stepAIC(mod2_2, direction = "both", trace = FALSE)
summary(step.model)

library(ggplot2)
bigxii %>%
  ggplot(aes(x = as.factor(team), y = WinMargin)) + 
  geom_boxplot()

library(lme4)
mod3 <- lmer(WinMargin ~ tRushOff + lRushDef + PassOff + lPassDef + tThirdDownPct + 
            Opp3rdDownPct + tPenYds + TimeOfPoss + TOMargin + (1|team),
             data = rbigxii)
summary(mod3)
plot(mod3)

mod3_lme <- lme(WinMargin ~ tRushOff + lRushDef + PassOff + lPassDef + tThirdDownPct + 
               Opp3rdDownPct + tPenYds + TimeOfPoss + TOMargin,
               random = ~ 1|team,
             data = rbigxii)
summary(mod3_lme)

intervals(mod3_lme)

library(merTools)
df <- as.data.frame(rbigxii[which(rbigxii$team == "Texas" & rbigxii$year == 2019),-14])
df2 <- splitstackshape::expandRows(df, count=29, count.is.col=FALSE)
texas19to <- bind_cols(df2, TOMargin = seq(-1.42, 1.46, by = 0.1))

preds <- merTools::predictInterval(mod3, newdata = texas19to, level = 0.95, n.sims = 999)
texas19to_preds <- bind_cols(preds, TOMargin = texas19to$TOMargin)

ggplot(data = texas19to_preds, aes(x = TOMargin, y = fit)) + 
  geom_smooth(color = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.6, fill = "#bf5700") + 
  geom_segment(aes(x = 0.38, y = 2.11, xend = 0.38, yend = 11.6)) + 
  geom_point(aes(x=0.38, y=7.7), color="black") + 
  labs(x = "Average Season Turnover Margin", y = "Predicted Average Season Win Margin", 
       title = "Texas Football 2019 Reimagined", 
       subtitle = "Win Margin Prediction Interval for Varying Turnover Margin",
       caption = "Holds all other variables constant") + 
  theme(text=element_text(size=17))
ggsave("texas2019to_preds.png")
