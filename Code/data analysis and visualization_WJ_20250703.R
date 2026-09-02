# 1. readying --------------------------------------------------------------------------------------
# __1.1 loading packages ---------------------------------------------------------------------------
pacman::p_load(openxlsx, tidyverse, ggrepel, ggtree, ggimage, RColorBrewer, cowplot, rtrees, lme4, lmerTest, glmmTMB, emmeans, performance, MuMIn, piecewiseSEM)

# __1.2 loading functions ---------------------------------------------------------------------------
# function 01: z-transform for data
scale_z <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# function 02: exported the lm results
save_lm <- function(index, count) {
  model_list <- str_c("model_", index, "_", str_pad(1:count, 2, side = "left", pad = "0"), sep = "")
  model_save <- NULL
  for (i in 1:n_distinct(model_list)) {
    Std <- effectsize::standardize_parameters(get(model_list[i]), ci = 0.90) |> tibble() |> filter(Parameter != "(Intercept)")
    model <- summary(get(model_list[i]))
    call <- as.character(model$call)
    ID <- c(1, nrow(model$coefficients))
    model_save <- bind_rows(model_save, data.frame(index = index, 
                                                   formula = c(str_c(call[1], "(formula = ", call[2], ", data = ", call[3], ")"), rep(NA, max(ID) - 1)), 
                                                   fix_name = c(row.names(model$coefficients), rep(NA, max(ID) - ID[2])), 
                                                   fix_estimate = c(model$coefficients[, 1], rep(NA, max(ID) - ID[2])), 
                                                   fix_se = c(model$coefficients[, 2], rep(NA, max(ID) - ID[2])), 
                                                   fix_p_value = c(model$coefficients[, 4], rep(NA, max(ID) - ID[2])), 
                                                   Std_Coefficient = c(NA, Std$Std_Coefficient), 
                                                   CI  = c(NA, Std$CI), 
                                                   CI_low = c(NA, Std$CI_low), 
                                                   CI_high = c(NA, Std$CI_high)))
  }
  row.names(model_save) <- NULL
  return(model_save)
}

# function 03: exported the GLMM results
save_glmm <- function(index, count) {
  model_list <- str_c("model_", index, "_", str_pad(1:count, 2, side = "left", pad = "0"), sep = "")
  model_save <- NULL
  for (i in 1:n_distinct(model_list)) {
    model <- summary(get(model_list[i]))
    call <- as.character(model$call)
    mult_comp <- pairs(emmeans(get(model_list[i]), as.formula(str_c("~ ", as.character(model$call$formula[[3]][2]))), adjust = "tukey")) |> as.data.frame()
    ID <- c(1, nrow(model$coefficients$cond), length(as.character(model$varcor)[1]), nrow(mult_comp))
    model_save <- bind_rows(model_save, data.frame(index = index, 
                                                   formula = c(str_c(call[1], "(formula = ", call[2], ", data = ", call[3], ", family = ", call[4], ")"), rep(NA, max(ID) - 1)), 
                                                   # AIC = c(model$AICtab[1], rep(NA, max(ID) - 1)), 
                                                   # BIC = c(model$AICtab[2], rep(NA, max(ID) - 1)), 
                                                   # LogLik = c(model$AICtab[3], rep(NA, max(ID) - 1)), 
                                                   # deviance = c(model$AICtab[4], rep(NA, max(ID) - 1)), 
                                                   # df_resid = c(model$AICtab[5], rep(NA, max(ID) - 1)), 
                                                   fix_name = c(row.names(model$coefficients$cond), rep(NA, max(ID) - ID[2])), 
                                                   fix_estimate = c(model$coefficients$cond[, 1], rep(NA, max(ID) - ID[2])), 
                                                   fix_se = c(model$coefficients$cond[, 2], rep(NA, max(ID) - ID[2])), 
                                                   fix_z_value = c(model$coefficients$cond[, 3], rep(NA, max(ID) - ID[2])), 
                                                   fix_p_value = c(model$coefficients$cond[, 4], rep(NA, max(ID) - ID[2])), 
                                                   # rand_group = c("site", rep(NA, max(ID) - ID[3])), 
                                                   # rand_name = c("intercept", rep(NA, max(ID) - ID[3])), 
                                                   # rand_var = c(str_extract(as.character(model$varcor)[1], "\\d+\\.\\d+"), rep(NA, max(ID) - ID[3])), 
                                                   # rand_sd = c(NA, rep(NA, max(ID) - ID[3])), 
                                                   R2_conditional = c(r2(get(model_list[i]))[[1]], rep(NA, max(ID) - 1)), 
                                                   R2_marginal = c(r2(get(model_list[i]))[[2]], rep(NA, max(ID) - 1))
                                                   # emmeans_contrast = c(mult_comp[, 1], rep(NA, max(ID) - ID[4])), 
                                                   # emmeans_estimate = c(mult_comp[, 2], rep(NA, max(ID) - ID[4])), 
                                                   # emmeans_se = c(mult_comp[, 3], rep(NA, max(ID) - ID[4])), 
                                                   # emmeans_df = c(mult_comp[, 4], rep(NA, max(ID) - ID[4])), 
                                                   # emmeans_z_ratio = c(mult_comp[, 5], rep(NA, max(ID) - ID[4])), 
                                                   # emmeans_p_value = c(mult_comp[, 6], rep(NA, max(ID) - ID[4]))
    ))
  }
  row.names(model_save) <- NULL
  return(model_save)
}

# function 04: bootstrap for slope
boot_slope <- function(x1, x2, y, n, data) {
  set.seed(1234); slope <- data.frame(slope1 = as.numeric(), slope2 = as.numeric())
  for(i in 1:n) {
    data_1 <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    data_1 <- data[data_1, ]
    data_2 <- data_1[c(x1, x2, "alti_log", y, "n_total", "site")]
    names(data_2) <- c("x1", "x2", "alti_log", "y", "n_total", "site")
    model <- glmmTMB(cbind(y, n_total - y) ~ x1 + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_2, family = binomial)
    std_slope <- effectsize::standardize_parameters(model) |> tibble()
    slope[i, 1] <- std_slope[2, 2]
    model <- glmmTMB(cbind(y, n_total - y) ~ x2 + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_2, family = binomial)
    std_slope <- effectsize::standardize_parameters(model) |> tibble()
    slope[i, 2] <- std_slope[2, 2]
  }
  return(slope)
}

# function 05: exported the rlm results
save_rlm <- function(index, count) {
  model_list <- str_c("model_", index, "_", str_pad(1:count, 2, side = "left", pad = "0"), sep = "")
  model_save <- NULL
  for (i in 1:n_distinct(model_list)) {
    model <- summary(get(model_list[i]))
    p_value <- 2 * pt(-abs(model$coefficients[, 3]), model$df[2])
    # F_test <- sfsmisc::f.robftest(get(model_list[i]), var = index)
    call <- as.character(model$call)
    ID <- c(1, nrow(model$coefficients))
    data_1 <- data[index] |> as.vector()
    model_save <- bind_rows(model_save, data.frame(index = index, 
                                                   formula = c(str_c(call[1], "(formula = ", call[2], ", data = ", call[3], ")"), rep(NA, max(ID) - 1)), 
                                                   fix_name = c(row.names(model$coefficients), rep(NA, max(ID) - ID[2])), 
                                                   fix_estimate = c(model$coefficients[, 1], rep(NA, max(ID) - ID[2])), 
                                                   fix_se = c(model$coefficients[, 2], rep(NA, max(ID) - ID[2])), 
                                                   t_value = c(model$coefficients[, 3], rep(NA, max(ID) - ID[2])), 
                                                   t_p_value = c(p_value, rep(NA, max(ID) - ID[2])), 
                                                   Std_Coefficient = c(NA, model$coefficients[2, 1] * sd(data_1[[1]], na.rm = TRUE))
                                                   # F_value = c(F_test$statistic, rep(NA, max(ID) - ID[2])), 
                                                   # F_p_value = c(F_test$p.value, rep(NA, max(ID) - ID[2]))
    ))
  }
  row.names(model_save) <- NULL
  return(model_save)
}

# __1.3 loading datasets ---------------------------------------------------------------------------
load(file = "data/data_set1.rdata")
load(file = "data/data_set2_1.rdata")
load(file = "data/data_set2_2.rdata")
load(file = "data/data_set3.rdata")

# 2. Part 1: Diversity resistance to freezing rain (FR) --------------------------------------------
data_set1 <- data_set1 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude))
# total damage      (DR_total):  n_damaged/n_total
# uprooting         (DR_uproot): n_uproot/n_total
# clear-bole broken (DR_below):  n_below/n_total
# crown broken      (DR_up):     n_up/n_total
# branch broken     (DR_branch): n_branch/n_total

# __2.1 Q1: Can taxonomic diversity resist the damage caused by FR? --------------------------------
model_save <- data.frame(); model_list <- list()

# ____2.1.1 Test 1: TD_SR --------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

models <- str_c("model_", "TD_SR", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# ______(2) DR_uproot ------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

models <- str_c("model_", "TD_SR", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# ______(3) DR_below -------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

models <- str_c("model_", "TD_SR", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# ______(4) DR_up ----------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

models <- str_c("model_", "TD_SR", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# ______(5) DR_branch ------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

models <- str_c("model_", "TD_SR", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# save model --
# save(model_list, file = "save/model_Part1_Q1_list_20250703.rdata")
# save(model_save, file = "save/model_Part1_Q1_data_20250703.rdata")
# write.xlsx(model_save, "save/model_Part1_Q1_data_20250703.xlsx")
rm(model, models, model_list, model_save, i, model_TD_SR_01, model_TD_SR_02, model_TD_SR_03, model_TD_SR_04, model_TD_SR_05, model_TD_SR_06)

# __2.2 Q2: Can structural diversity resist the damage caused by FR? -------------------------------
model_save <- data.frame(); model_list <- list()

# ____2.2.1 Test 1: CV_TTH -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

models <- str_c("model_", "CV_TTH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ______(2) DR_uproot ------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

models <- str_c("model_", "CV_TTH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ______(3) DR_below -------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

models <- str_c("model_", "CV_TTH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ______(4) DR_up ----------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

models <- str_c("model_", "CV_TTH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ______(5) DR_branch ------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

models <- str_c("model_", "CV_TTH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ____2.2.2 Test 2: CV_DBH -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

models <- str_c("model_", "CV_DBH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ______(2) DR_uproot ------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

models <- str_c("model_", "CV_DBH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ______(3) DR_below -------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

models <- str_c("model_", "CV_DBH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ______(4) DR_up ----------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

models <- str_c("model_", "CV_DBH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ______(5) DR_branch ------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

models <- str_c("model_", "CV_DBH", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ____2.2.3 Test 3: CV_NMB -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

models <- str_c("model_", "CV_NMB", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# ______(2) DR_uproot ------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

models <- str_c("model_", "CV_NMB", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# ______(3) DR_below -------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

models <- str_c("model_", "CV_NMB", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# ______(4) DR_up ----------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

models <- str_c("model_", "CV_NMB", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# ______(5) DR_branch ------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

models <- str_c("model_", "CV_NMB", "_", str_pad(1:6, 2, side = "left", pad = "0")); for (i in 1:length(models)) {model <- get(models[i]); model_list <- append(model_list, list(model))}
model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# save model --
# save(model_list, file = "save/model_Part1_Q2_list_20250703.rdata")
# save(model_save, file = "save/model_Part1_Q2_data_20250703.rdata")
# write.xlsx(model_save, "save/model_Part1_Q2_data_20250703.xlsx")
rm(model, models, model_list, model_save, i, 
   model_CV_TTH_01, model_CV_TTH_02, model_CV_TTH_03, model_CV_TTH_04, model_CV_TTH_05, model_CV_TTH_06, 
   model_CV_DBH_01, model_CV_DBH_02, model_CV_DBH_03, model_CV_DBH_04, model_CV_DBH_05, model_CV_DBH_06, 
   model_CV_NMB_01, model_CV_NMB_02, model_CV_NMB_03, model_CV_NMB_04, model_CV_NMB_05, model_CV_NMB_06)

# __2.3 Q3: Do taxonomic and structural diversity differ in their resistance to FR? ----------------
# check colinearity --
model_save <- data.frame()

# ____2.3.1 Test 1: TD_SR vs CV_TTH ----------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_damaged", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_total_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_damaged", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_total_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_damaged", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_total_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_total", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(2) DR_uproot ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_uproot", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_uproot_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_uproot", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_uproot_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_uproot", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_uproot_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_uproot", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(3) DR_below -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_below", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_below_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_below", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_below_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_below", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_below_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_below", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(4) DR_up ----------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_up", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_up_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_up", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_up_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_up", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_up_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_up", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(5) DR_branch ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_branch", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_branch_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_branch", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_branch_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_branch", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_branch_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_branch", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ____2.3.2 Test 2: TD_SR vs CV_DBH ----------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_damaged", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_total_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_damaged", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_total_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_damaged", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_total_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_total", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(2) DR_uproot ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_uproot", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_uproot_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_uproot", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_uproot_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_uproot", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_uproot_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_uproot", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(3) DR_below -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_below", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_below_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_below", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_below_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_below", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_below_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_below", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(4) DR_up ----------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_up", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_up_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_up", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_up_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_up", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_up_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_up", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(5) DR_branch ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_branch", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_branch_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_branch", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_branch_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_branch", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_branch_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_branch", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ____2.3.3 Test 3: TD_SR vs CV_NMB ----------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_damaged", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_total_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_damaged", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_total_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_damaged", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_total_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_total", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(2) DR_uproot ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_uproot", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_uproot_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_uproot", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_uproot_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_uproot", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_uproot_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_uproot", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(3) DR_below -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_below", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_below_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_below", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_below_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_below", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_below_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_below", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(4) DR_up ----------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_up", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_up_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_up", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_up_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_up", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_up_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_up", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(5) DR_branch ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_branch", 1000, data_set1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/DR_branch_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_branch", 1000, data_set1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/DR_branch_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_branch", 1000, data_set1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/DR_branch_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_branch", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# save model --
save(model_save, file = "save/model_Part1_Q3_data_20250703.rdata")
write.xlsx(model_save, "save/model_Part1_Q3_data_20250703.xlsx")
rm(data_slope, slope, slope_test_1, slope_test_2, slope_test_3, model_save, model_1, model_2, model_3, model_4, model_5, model_6)

# 3. Part 2: Species-level structural diversity in resistance to FR --------------------------------
# __3.1 Q1: Does the structural diversity of species resist FR? ------------------------------------
data <- data_set2_1 |> bind_rows() |> group_by(species_CN) |> filter(DR_total != 0) |> mutate(n_plot = n()) |> filter(n_plot > 5) |> ungroup()
model_save_rlm <- data.frame()
list_sp <- data |> distinct(species_CN)

# ____3.1.1 Test 1: CV_TTH -------------------------------------------------------------------------
for (i in 1:nrow(list_sp)) {
  data_1 <- data |> filter(species_CN == as.character(list_sp[i, 1])) |> drop_na(CV_TTH)
  model_CV_TTH_01 <- MASS::rlm(DR_total ~ CV_TTH, data_1)
  model_CV_TTH_02 <- MASS::rlm(DR_uproot ~ CV_TTH, data_1)
  model_CV_TTH_03 <- MASS::rlm(DR_below ~ CV_TTH, data_1)
  model_CV_TTH_04 <- MASS::rlm(DR_up ~ CV_TTH, data_1)
  model_CV_TTH_05 <- MASS::rlm(DR_branch ~ CV_TTH, data_1)
  model_save_rlm <- bind_rows(model_save_rlm, save_rlm("CV_TTH", 5) |> mutate(species = as.character(list_sp[i, 1]), .after = index))
}

# ____3.1.2 Test 2: CV_DBH -------------------------------------------------------------------------
for (i in 1:nrow(list_sp)) {
  data_1 <- data |> filter(species_CN == as.character(list_sp[i, 1])) |> drop_na(CV_DBH)
  model_CV_DBH_01 <- MASS::rlm(DR_total ~ CV_DBH, data_1)
  model_CV_DBH_02 <- MASS::rlm(DR_uproot ~ CV_DBH, data_1)
  model_CV_DBH_03 <- MASS::rlm(DR_below ~ CV_DBH, data_1)
  model_CV_DBH_04 <- MASS::rlm(DR_up ~ CV_DBH, data_1)
  model_CV_DBH_05 <- MASS::rlm(DR_branch ~ CV_DBH, data_1)
  model_save_rlm <- bind_rows(model_save_rlm, save_rlm("CV_DBH", 5) |> mutate(species = as.character(list_sp[i, 1]), .after = index))
}

# ____3.1.3 Test 3: CV_NMB -------------------------------------------------------------------------
for (i in 1:nrow(list_sp)) {
  data_1 <- data |> filter(species_CN == as.character(list_sp[i, 1])) |> drop_na(CV_NMB)
  model_CV_NMB_01 <- MASS::rlm(DR_total ~ CV_NMB, data_1)
  model_CV_NMB_02 <- MASS::rlm(DR_uproot ~ CV_NMB, data_1)
  model_CV_NMB_03 <- MASS::rlm(DR_below ~ CV_NMB, data_1)
  model_CV_NMB_04 <- MASS::rlm(DR_up ~ CV_NMB, data_1)
  model_CV_NMB_05 <- MASS::rlm(DR_branch ~ CV_NMB, data_1)
  model_save_rlm <- bind_rows(model_save_rlm, save_rlm("CV_NMB", 5) |> mutate(species = as.character(list_sp[i, 1]), .after = index))
}

# save model --
# save(model_save_rlm, file = "save/model_Part2_Q1_data_20250703.rdata")
# write.xlsx(model_save_rlm, "save/model_Part2_Q1_data_20250703.xlsx")
rm(model_CV_TTH_01, model_CV_TTH_02, model_CV_TTH_03, model_CV_TTH_04, model_CV_TTH_05, 
   model_CV_DBH_01, model_CV_DBH_02, model_CV_DBH_03, model_CV_DBH_04, model_CV_DBH_05, 
   model_CV_NMB_01, model_CV_NMB_02, model_CV_NMB_03, model_CV_NMB_04, model_CV_NMB_05, 
   i, data, data_1, list_sp)

# __3.2 Q2: What species traits determine the response to FR? --------------------------------------
model_save <- data.frame(); model_multi <- data.frame()

# ____3.2.1 Test 1: CV_TTH -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_total") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_total") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_TTH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_total", index = "CV_TTH"))

# ______(2) DR_uproot ------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_uproot") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_uproot") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_TTH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_uproot", index = "CV_TTH"))

# ______(3) DR_below -------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_below") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_below") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_TTH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_below", index = "CV_TTH"))

# ______(4) DR_up ----------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_up") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_up") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_TTH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_up", index = "CV_TTH"))

# ______(5) DR_branch ------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_branch") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_branch") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_TTH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_branch", index = "CV_TTH"))

# ____3.2.2 Test 2: CV_DBH -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_total") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_total") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_DBH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_total", index = "CV_DBH"))

# ______(2) DR_uproot ------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_uproot") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_uproot") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_DBH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_uproot", index = "CV_DBH"))

# ______(3) DR_below -------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_below") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_below") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_DBH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_below", index = "CV_DBH"))

# ______(4) DR_up ----------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_up") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_up") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_DBH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_up", index = "CV_DBH"))

# ______(5) DR_branch ------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_branch") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_branch") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_DBH", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_branch", index = "CV_DBH"))

# ____3.2.3 Test 3: CV_NMB -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_total") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_total") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_NMB", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_total", index = "CV_NMB"))

# ______(2) DR_uproot ------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_uproot") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_uproot") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_NMB", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_uproot", index = "CV_NMB"))

# ______(3) DR_below -------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_below") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_below") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_NMB", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_below", index = "CV_NMB"))

# ______(4) DR_up ----------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_up") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_up") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_NMB", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_up", index = "CV_NMB"))

# ______(5) DR_branch ------------------------------------------------------------------------------
data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_branch") |> pivot_longer(cols = c(woody_density:niche_width), names_to = "trait", values_to = "value") |> 
  ggplot(aes(x = value, y = Std_Coefficient, color = trait)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color = "#ff0000") + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ trait, scales = "free") + 
  labs(x = NULL, y = "Standardized slopes", color = "Traits") + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none")

data <- data_set2_2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_branch") |> mutate(across(c(woody_density:niche_width), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ niche_width, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)
model_save <- bind_rows(model_save, save_lm("CV_NMB", 5))

model_all <- lm(Std_Coefficient ~ niche_width + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_multi <- bind_rows(model_multi, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_branch", index = "CV_NMB"))

# save model --
# save(model_save, file = "save/model_Part2_Q2_data_20250703.rdata")
# write.xlsx(model_save, "save/model_Part2_Q2_data_20250703.xlsx")
# model_multi |> rownames_to_column("fix_name") |> mutate(fix_name = str_remove(fix_name, "\\...[0-9]*")) |> write.xlsx("save/model_Part2_Q2_multi_data_20250703.xlsx")
rm(model_CV_TTH_01, model_CV_TTH_02, model_CV_TTH_03, model_CV_TTH_04, model_CV_TTH_05, 
   model_CV_DBH_01, model_CV_DBH_02, model_CV_DBH_03, model_CV_DBH_04, model_CV_DBH_05, 
   model_CV_NMB_01, model_CV_NMB_02, model_CV_NMB_03, model_CV_NMB_04, model_CV_NMB_05, 
   model_all, model_set, model_avg, model_save, model_multi, data)

# 4. Part 3: The potential mechanism of structural diversity in resisting FR -----------------------
# __4.1 structural equation model ------------------------------------------------------------------
model_save <- data.frame()

# ____4.1.1 Test 1: DR_total -----------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_set3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> mutate(forest_age = case_when(forest_age == "OF" ~ 3, forest_age == "LSF" ~ 2, TRUE ~ 1)) |> 
  mutate(across(c(forest_age:alti_log), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_total ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_total ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Forest stands --
model <- glmmTMB(DR_total ~ forest_age + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stan <- data$forest_age

# (3) construct composite variable Species compositiona --
model <- glmmTMB(DR_total ~ TD_SR + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$comp <- data$TD_SR

# (4) construct composite variable Structural diversity --
model <- glmmTMB(DR_total ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_total ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (5) construct composite variable Functional identity --
model <- glmmTMB(DR_total ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_total ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$stan, site = data$site, DR = data$DR_total, alti1 = data$alti1, altiq = data$altiq, comp = data$comp, stru = data$stru, func = data$func)

# full model --
model_sem_0 <- psem(glmmTMB(DR ~ altiq + stan + comp + stru + func + (1 | site), data_sem), 
                    glmmTMB(stru ~ alti1 + stan + comp + (1 | site), data_sem), 
                    glmmTMB(func ~ alti1 + stan + comp + (1 | site), data_sem), 
                    glmmTMB(comp ~ alti1 + stan + (1 | site), data_sem), 
                    func %~~% stru); summary(model_sem_0)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_0)$coefficients[, -9] |> mutate(type = "total", DR = "DR_total"), 
                                              bind_rows(summary(model_sem_0)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - nrow(summary(model_sem_0)$R2)))) |> select(!del)))

# part model within ESF --
data_sem1 <- data_sem |> filter(stan < 0) |> select(site:func)
model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_sem_1 <- psem(glmmTMB(DR ~ altiq + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ altiq + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ altiq + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_total"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan > 0 & stan < 1.5) |> select(site:func)
model_sem_2 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_total"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan > 1.5) |> select(site:func)
model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_sem_3 <- psem(glmmTMB(DR ~ altiq + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_total"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))

rm(data, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_0, model_sem_1, model_sem_2, model_sem_3)

# ____4.1.2 Test 2: DR_uproot ----------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_set3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> mutate(forest_age = case_when(forest_age == "OF" ~ 3, forest_age == "LSF" ~ 2, TRUE ~ 1)) |> 
  mutate(across(c(forest_age:alti_log), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_uproot ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_uproot ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Forest stands --
model <- glmmTMB(DR_uproot ~ forest_age + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stan <- data$forest_age

# (3) construct composite variable Species compositiona --
model <- glmmTMB(DR_uproot ~ TD_SR + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$comp <- data$TD_SR

# (4) construct composite variable Structural diversity --
model <- glmmTMB(DR_uproot ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_uproot ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (5) construct composite variable Functional identity --
model <- glmmTMB(DR_uproot ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_uproot ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$stan, site = data$site, DR = data$DR_uproot, alti1 = data$alti1, altiq = data$altiq, comp = data$comp, stru = data$stru, func = data$func)

# full model --
model_sem_0 <- psem(glmmTMB(DR ~ altiq + stan + comp + stru + func + (1 | site), data_sem), 
                    glmmTMB(stru ~ alti1 + stan + comp + (1 | site), data_sem), 
                    glmmTMB(func ~ altiq + stan + comp + (1 | site), data_sem), 
                    glmmTMB(comp ~ alti1 + stan + (1 | site), data_sem), 
                    func %~~% stru); summary(model_sem_0)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_0)$coefficients[, -9] |> mutate(type = "total", DR = "DR_uproot"), 
                                              bind_rows(summary(model_sem_0)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - nrow(summary(model_sem_0)$R2)))) |> select(!del)))

# part model within ESF --
data_sem1 <- data_sem |> filter(stan < 0) |> select(site:func)
model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_uproot"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan > 0 & stan < 1.5) |> select(site:func)
model_sem_2 <- psem(glmmTMB(DR ~ altiq + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_uproot"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan > 1.5) |> select(site:func)
model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_uproot"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))

rm(data, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_0, model_sem_1, model_sem_2, model_sem_3)

# ____4.1.3 Test 3: DR_below -----------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_set3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> mutate(forest_age = case_when(forest_age == "OF" ~ 3, forest_age == "LSF" ~ 2, TRUE ~ 1)) |> 
  mutate(across(c(forest_age:alti_log), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_below ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_below ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Forest stands --
model <- glmmTMB(DR_below ~ forest_age + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stan <- data$forest_age

# (3) construct composite variable Species compositiona --
model <- glmmTMB(DR_below ~ TD_SR + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$comp <- data$TD_SR

# (4) construct composite variable Structural diversity --
model <- glmmTMB(DR_below ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_below ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (5) construct composite variable Functional identity --
model <- glmmTMB(DR_below ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_below ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$stan, site = data$site, DR = data$DR_below, alti1 = data$alti1, altiq = data$altiq, comp = data$comp, stru = data$stru, func = data$func)

# full model --
model_sem_0 <- psem(glmmTMB(DR ~ altiq + stan + comp + stru + func + (1 | site), data_sem), 
                    glmmTMB(stru ~ altiq + stan + comp + (1 | site), data_sem), 
                    glmmTMB(func ~ alti1 + stan + comp + (1 | site), data_sem), 
                    glmmTMB(comp ~ alti1 + stan + (1 | site), data_sem), 
                    func %~~% stru); summary(model_sem_0)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_0)$coefficients[, -9] |> mutate(type = "total", DR = "DR_below"), 
                                              bind_rows(summary(model_sem_0)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - nrow(summary(model_sem_0)$R2)))) |> select(!del)))

# part model within ESF --
data_sem1 <- data_sem |> filter(stan < 0) |> select(site:func)
model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ altiq + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_below"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan > 0 & stan < 1.5) |> select(site:func)
model_sem_2 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_below"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan > 1.5) |> select(site:func)
model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_below"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))

rm(data, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_0, model_sem_1, model_sem_2, model_sem_3)

# ____4.1.4 Test 4: DR_up --------------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_set3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> mutate(forest_age = case_when(forest_age == "OF" ~ 3, forest_age == "LSF" ~ 2, TRUE ~ 1)) |> 
  mutate(across(c(forest_age:alti_log), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_up ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_up ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Forest stands --
model <- glmmTMB(DR_up ~ forest_age + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stan <- data$forest_age

# (3) construct composite variable Species compositiona --
model <- glmmTMB(DR_up ~ TD_SR + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$comp <- data$TD_SR

# (4) construct composite variable Structural diversity --
model <- glmmTMB(DR_up ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_up ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (5) construct composite variable Functional identity --
model <- glmmTMB(DR_up ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_up ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$stan, site = data$site, DR = data$DR_up, alti1 = data$alti1, altiq = data$altiq, comp = data$comp, stru = data$stru, func = data$func)

# full model --
model_sem_0 <- psem(glmmTMB(DR ~ altiq + stan + comp + stru + func + (1 | site), data_sem), 
                    glmmTMB(stru ~ alti1 + stan + comp + (1 | site), data_sem), 
                    glmmTMB(func ~ alti1 + stan + comp + (1 | site), data_sem), 
                    glmmTMB(comp ~ alti1 + stan + (1 | site), data_sem), 
                    func %~~% stru); summary(model_sem_0)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_0)$coefficients[, -9] |> mutate(type = "total", DR = "DR_up"), 
                                              bind_rows(summary(model_sem_0)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - nrow(summary(model_sem_0)$R2)))) |> select(!del)))

# part model within ESF --
data_sem1 <- data_sem |> filter(stan < 0) |> select(site:func)
model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_up"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan > 0 & stan < 1.5) |> select(site:func)
model_sem_2 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_up"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan > 1.5) |> select(site:func)
model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_up"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))

rm(data, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_0, model_sem_1, model_sem_2, model_sem_3)

# ____4.1.5 Test 5: DR_branch ----------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_set3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> mutate(forest_age = case_when(forest_age == "OF" ~ 3, forest_age == "LSF" ~ 2, TRUE ~ 1)) |> 
  mutate(across(c(forest_age:alti_log), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_branch ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_branch ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Forest stands --
model <- glmmTMB(DR_branch ~ forest_age + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stan <- data$forest_age

# (3) construct composite variable Species compositiona --
model <- glmmTMB(DR_branch ~ TD_SR + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$comp <- data$TD_SR

# (4) construct composite variable Structural diversity --
model <- glmmTMB(DR_branch ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_branch ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (5) construct composite variable Functional identity --
model <- glmmTMB(DR_branch ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_branch ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$stan, site = data$site, DR = data$DR_branch, alti1 = data$alti1, altiq = data$altiq, comp = data$comp, stru = data$stru, func = data$func)

# full model --
model_sem_0 <- psem(glmmTMB(DR ~ stan + comp + stru + func + (1 | site), data_sem), 
                    glmmTMB(stru ~ alti1 + stan + comp + (1 | site), data_sem), 
                    glmmTMB(func ~ alti1 + stan + comp + (1 | site), data_sem), 
                    glmmTMB(comp ~ alti1 + stan + (1 | site), data_sem), 
                    func %~~% stru); summary(model_sem_0)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_0)$coefficients[, -9] |> mutate(type = "total", DR = "DR_branch"), 
                                              bind_rows(summary(model_sem_0)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_0)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_0)$coefficients) - nrow(summary(model_sem_0)$R2)))) |> select(!del)))

# part model within ESF --
data_sem1 <- data_sem |> filter(stan < 0) |> select(site:func)
model_sem_1 <- psem(glmmTMB(DR ~ comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_branch"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan > 0 & stan < 1.5) |> select(site:func)
model_sem_2 <- psem(glmmTMB(DR ~ comp + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_branch"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan > 1.5) |> select(site:func)
model_sem_3 <- psem(glmmTMB(DR ~ comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_branch"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))

model_save <- model_save |> set_names("Response_1", "Predictor", "Estimate", "Std_Error", "DF", "Crit_Value", "P_Value_1", "Std_Estimate", "Type", "DR", "Chisq", "df_1", "P_Value_2", "Fisher_C", "df_2", "P_Value_3", "Response_2", "family", "link", "method", "Marginal", "Conditional")
# save(model_save, file = "save/model_Part3_Q1_data_20250703.rdata")
# write.xlsx(model_save, "save/model_Part3_Q1_data_20250703.xlsx")
rm(data, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_0, model_sem_1, model_sem_2, model_sem_3, model_save)

# 5. data visualization ----------------------------------------------------------------------------
# __5.1 Figure 01 ----------------------------------------------------------------------------------
# use PPT for drawing --

# __5.2 Figure 02 ----------------------------------------------------------------------------------
# data processing --
load(file = "data/data_set1.rdata")
data_set1 <- data_set1 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude))

# ____5.2.1 Fig. 2A --------------------------------------------------------------------------------
# TD_SR vs DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "ESF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data <- effects::allEffects(model, xlevels = list(TD_SR = seq(1, 5, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "LSF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(TD_SR = seq(1, 5, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF"))

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "OF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(TD_SR = seq(1, 5, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF"))

data <- data |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_01 <- ggplot() + 
  geom_ribbon(data = data, aes(x = TD_SR, y = fit, ymin = lower, ymax = upper, fill = type), alpha = 0.2) + 
  geom_line(data = data, aes(x = TD_SR, y = fit, color = type, linetype = type)) + 
  scale_x_continuous(breaks = c(1, 3, 5), labels = c(2, 8, 32)) +
  scale_y_continuous(limits = c(5, 55), breaks = c(10, 30, 50)) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_linetype_manual(values = c("longdash", "solid", "solid")) + 
  labs(x = "SR", y = NULL, fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.94), label = expression("OF:   "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.86), label = expression("LSF: "*italic(P)*" = 0.023"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.78), label = expression("ESF: "*italic(P)*" = 0.787"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") +
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_01

# ____5.2.2 Fig. 2B --------------------------------------------------------------------------------
# CV_TTH vs DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data <- effects::allEffects(model, xlevels = list(CV_TTH = seq(16.3, 62.3, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_TTH = seq(16.3, 62.3, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF"))

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "OF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_TTH = seq(16.3, 62.3, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF"))

data <- data |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_02 <- ggplot() + 
  geom_ribbon(data = data, aes(x = CV_TTH, y = fit, ymin = lower, ymax = upper, fill = type), alpha = 0.2) + 
  geom_line(data = data, aes(x = CV_TTH, y = fit, color = type, linetype = type)) + 
  scale_x_continuous(breaks = c(18, 40, 62)) + 
  scale_y_continuous(limits = c(5, 55), breaks = c(10, 30, 50)) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_linetype_manual(values = c("longdash", "solid", "solid")) + 
  labs(x = expression(CV[H]*" (%)"), y = NULL, fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.94), label = expression("OF:   "*italic(P)*" < 0.001"),  size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.86), label = expression("LSF: "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.78), label = expression("ESF: "*italic(P)*" = 0.999"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") +
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_02

# ____5.2.3 Fig. 2C --------------------------------------------------------------------------------
# CV_DBH vs DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data <- effects::allEffects(model, xlevels = list(CV_DBH = seq(29, 139, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_DBH = seq(29, 139, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF"))

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "OF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_DBH = seq(29, 139, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF"))

data <- data |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_03 <- ggplot() + 
  geom_ribbon(data = data, aes(x = CV_DBH, y = fit, ymin = lower, ymax = upper, fill = type), alpha = 0.2) + 
  geom_line(data = data, aes(x = CV_DBH, y = fit, color = type, linetype = type)) + 
  scale_x_continuous(breaks = c(35, 85, 135)) + 
  scale_y_continuous(limits = c(5, 55), breaks = c(10, 30, 50)) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_linetype_manual(values = c("longdash", "solid", "solid")) + 
  labs(x = expression(CV[DBH]*" (%)"), y = NULL, fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.94), label = expression("OF:   "*italic(P)*" < 0.001"),  size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.86), label = expression("LSF: "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.78), label = expression("ESF: "*italic(P)*" = 0.474"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") +
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_03

# ____5.2.4 Fig. 2D --------------------------------------------------------------------------------
# CV_NMB vs DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data <- effects::allEffects(model, xlevels = list(CV_NMB = seq(39.7, 252, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_NMB = seq(39.7, 252, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF"))

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "OF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_NMB = seq(39.7, 252, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF"))

data <- data |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_04 <- ggplot() + 
  geom_ribbon(data = data, aes(x = CV_NMB, y = fit, ymin = lower, ymax = upper, fill = type), alpha = 0.2) + 
  geom_line(data = data, aes(x = CV_NMB, y = fit, color = type)) + 
  scale_x_continuous(breaks = c(50, 150, 250)) + 
  scale_y_continuous(limits = c(5, 55), breaks = c(10, 30, 50)) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  labs(x = expression(CV[NMB]*" (%)"), y = NULL, fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.94), label = expression("OF:   "*italic(P)*" < 0.001"),  size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.86), label = expression("LSF: "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.78), label = expression("ESF: "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") +
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_04

# ____5.2.5 Fig. 2E-G ------------------------------------------------------------------------------
load(file = "save/DR_total_CV_TTH_ESF.rdata"); data_1 <- data_slope
load(file = "save/DR_total_CV_DBH_ESF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope3"))
load(file = "save/DR_total_CV_NMB_ESF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope4"))

load(file = "save/DR_total_CV_TTH_LSF.rdata"); data_1 <- bind_cols(data_1, data_slope |> set_names("slope5", "slope6"))
load(file = "save/DR_total_CV_DBH_LSF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope7"))
load(file = "save/DR_total_CV_NMB_LSF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope8"))

load(file = "save/DR_total_CV_TTH_OF.rdata"); data_1 <- bind_cols(data_1, data_slope |> set_names("slope9", "slope10"))
load(file = "save/DR_total_CV_DBH_OF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope11"))
load(file = "save/DR_total_CV_NMB_OF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope12"))

data_1 <- data_1 |> pivot_longer(cols = slope1:slope12, names_to = "slope", values_to = "value") |> mutate(n = rep(4:1, 3000)) |> 
  mutate(X = case_when(slope %in% c("slope1", "slope5", "slope9") ~ "TD_SR", slope %in% c("slope2", "slope6", "slope10") ~ "CV_TTH", slope %in% c("slope3", "slope7", "slope11") ~ "CV_DBH", TRUE ~ "CV_NMB")) |> 
  mutate(forest_age = case_when(slope %in% c("slope1", "slope2", "slope3", "slope4") ~ "young", slope %in% c("slope5", "slope6", "slope7", "slope8") ~ "middle", TRUE ~ "mature")) |> 
  mutate(x = case_when((forest_age == "young") & (n == 1) ~ 1, (forest_age == "young") & (n == 2) ~ 2, (forest_age == "young") & (n == 3) ~ 3, (forest_age == "young") & (n == 4) ~ 4, 
                       (forest_age == "middle") & (n == 1) ~ 6, (forest_age == "middle") & (n == 2) ~ 7, (forest_age == "middle") & (n == 3) ~ 8, (forest_age == "middle") & (n == 4) ~ 9, 
                       (forest_age == "mature") & (n == 1) ~ 11, (forest_age == "mature") & (n == 2) ~ 12, (forest_age == "mature") & (n == 3) ~ 13, (forest_age == "mature") & (n == 4) ~ 14))

data_2 <- read.xlsx("save/3_statistical result/model_Part1_Q3_data_20250703.xlsx") |> tibble()
data_2 <- data_2 |> filter(Y == "DR_total") |> mutate(n = case_when(Parameter == "TD_SR" ~ 4, Parameter == "CV_TTH" ~ 3, Parameter == "CV_DBH" ~ 2, TRUE ~ 1)) |> 
  mutate(forest_age = case_when(str_detect(type, "ESF") ~ "young", str_detect(type, "LSF") ~ "middle", TRUE ~ "mature")) |> 
  mutate(x = case_when((forest_age == "young") & (n == 1) ~ 1, (forest_age == "young") & (n == 2) ~ 2, (forest_age == "young") & (n == 3) ~ 3, (forest_age == "young") & (n == 4) ~ 4, 
                       (forest_age == "middle") & (n == 1) ~ 6, (forest_age == "middle") & (n == 2) ~ 7, (forest_age == "middle") & (n == 3) ~ 8, (forest_age == "middle") & (n == 4) ~ 9, 
                       (forest_age == "mature") & (n == 1) ~ 11, (forest_age == "mature") & (n == 2) ~ 12, (forest_age == "mature") & (n == 3) ~ 13, (forest_age == "mature") & (n == 4) ~ 14))

plot_05 <- ggplot() + 
  gghalves::geom_half_violin(data = data_1, aes(x = x, y = value, group = x, fill = X, color = X), side = "r", alpha = 0.3, linewidth = 0.25) + 
  geom_errorbar(data = data_2, aes(x = x, ymin = CI_low, ymax = CI_high, color = Parameter), alpha = 0.5, width = 0, linewidth = 1, show.legend = FALSE) + 
  geom_point(data = data_2, aes(x = x, y = Std_Coefficient, color = Parameter), alpha = 0.5, shape = 16, size = 1.5, show.legend = FALSE) + 
  
  geom_rect(aes(xmin = 4.7, xmax = 5.3, ymin = -0.5, ymax = 0.5), fill = "#8cc269", alpha = 0.5) + 
  geom_segment(aes(x = 0.5, xend = 4.7, y = 0, yend = 0), linetype = 5, linewidth = 0.4, alpha = 0.7) + 
  annotate("text", x = 5, y = 0, label = "ESF", color = "#000000", size = (9*0.35), family = "serif") +
  geom_segment(aes(x = 1, y = 0.40, xend = 4, yend = 0.40), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 2.5, y = 0.45, label = "***", color = "#000000", size = (12*0.35), family = "serif", angle = 90) +
  
  geom_rect(aes(xmin = 9.7, xmax = 10.3, ymin = -0.5, ymax = 0.5), fill = "#1ba784", alpha = 0.5) + 
  geom_segment(aes(x = 5.3, xend = 9.7, y = 0, yend = 0), linetype = 5, linewidth = 0.4, alpha = 0.7) + 
  annotate("text", x = 10, y = 0, label = "LSF", color = "#000000", size = (9*0.35), family = "serif") +
  geom_segment(aes(x = 6, y = 0.40, xend = 9, yend = 0.40), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 7.5, y = 0.45, label = "***", color = "#000000", size = (12*0.35), family = "serif", angle = 90) +
  geom_segment(aes(x = 7, y = 0.32, xend = 9, yend = 0.32), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 8, y = 0.37, label = "***", color = "#000000", size = (12*0.35), family = "serif", angle = 90) +
  geom_segment(aes(x = 8, y = 0.24, xend = 9, yend = 0.24), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 8.5, y = 0.29, label = "***", color = "#000000", size = (12*0.35), family = "serif", angle = 90) +
  
  geom_rect(aes(xmin = 14.7, xmax = 15.3, ymin = -0.5, ymax = 0.5), fill = "#1a6840", alpha = 0.5) + 
  geom_segment(aes(x = 10.3, xend = 14.7, y = 0, yend = 0), linetype = 5, linewidth = 0.4, alpha = 0.7) + 
  annotate("text", x = 15, y = 0, label = "OF", color = "#000000", size = (9*0.35), family = "serif") +
  geom_segment(aes(x = 11, y = 0.40, xend = 14, yend = 0.40), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 12.5, y = 0.45, label = "***", color = "#000000", size = (12*0.35), family = "serif", angle = 90) +
  geom_segment(aes(x = 12, y = 0.32, xend = 14, yend = 0.32), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 13, y = 0.37, label = "***", color = "#000000", size = (12*0.35), family = "serif", angle = 90) +
  geom_segment(aes(x = 13, y = 0.24, xend = 14, yend = 0.24), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 13.5, y = 0.29, label = "***", color = "#000000", size = (12*0.35), family = "serif", angle = 90) +
  
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, 15.3), breaks = c(14, 13, 12, 11, 9, 8, 7, 6, 4, 3, 2, 1), labels = c("SR", expression(CV[H]), expression(CV[DBH]), expression(CV[NMB]), "SR", expression(CV[H]), expression(CV[DBH]), expression(CV[NMB]), "SR", expression(CV[H]), expression(CV[DBH]), expression(CV[NMB]))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 0.5), breaks = c(-0.4, 0, 0.4), labels = c("−0.4", "0", "0.4")) + 
  scale_fill_manual(values = c("#43b244", "#43b244", "#43b244", "#1781b5")) + 
  scale_color_manual(values = c("#43b244", "#43b244", "#43b244", "#1781b5")) + 
  labs(x = NULL, y = "Standardized slopes") + 
  coord_flip() + 
  theme_bw(base_family = "serif") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        strip.text = element_blank(),  
        legend.position = "none"); plot_05

# ____5.2.6 layout ---------------------------------------------------------------------------------
plot_06 <- ggplot() + scale_x_continuous(expand = c(0, 0), limits = c(0, 0.6)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0), panel.grid = element_blank(), panel.background = element_rect(fill = "#ffffff")) + 
  annotate(geom = "text", x = 0.3, y = 5.5, label = "Total damage ratio (%)", size = (10*0.35), family = "serif", angle = 90); plot_06

plot_grid(plot_grid(plot_06, plot_06, nrow = 2), 
  plot_grid(plot_01, plot_02, nrow = 2) + 
    annotate(geom = "text", x = ggpp::as_npc(0.05), y = ggpp::as_npc(0.98), label = "a", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.05), y = ggpp::as_npc(0.48), label = "b", size = (10*0.35), family = "serif", fontface = "bold"), 
  plot_grid(plot_03, plot_04, nrow = 2) + 
    annotate(geom = "text", x = ggpp::as_npc(0.05), y = ggpp::as_npc(0.98), label = "c", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.05), y = ggpp::as_npc(0.48), label = "d", size = (10*0.35), family = "serif", fontface = "bold"), 
  plot_05, rel_widths = c(0.6, 5.7, 5.7, 6.0), nrow = 1) + 
  annotate(geom = "text", x = ggpp::as_npc(0.68), y = ggpp::as_npc(0.975), label = "e", size = (10*0.35), family = "serif", fontface = "bold") + 
  annotate(geom = "text", x = ggpp::as_npc(0.68), y = ggpp::as_npc(0.670), label = "f", size = (10*0.35), family = "serif", fontface = "bold") + 
  annotate(geom = "text", x = ggpp::as_npc(0.68), y = ggpp::as_npc(0.363), label = "g", size = (10*0.35), family = "serif", fontface = "bold")

# ggsave(file = "save/fig_2_20250703.tiff", width = 16, height = 11, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")
rm(plot_01, plot_02, plot_03, plot_04, plot_05, plot_06, data, data_1, data_2, data_slope, model)

# __5.3 Figure 03 ----------------------------------------------------------------------------------
# ____5.3.1 Fig. 3A --------------------------------------------------------------------------------
# data processing --
load(file = "data/data_set2_1.rdata")
load(file = "save/model_Part2_Q1_data_20250703.rdata")

data <- model_save_rlm |> filter(!is.na(Std_Coefficient)) |> mutate(damage_type = rep(c("DR_total", "DR_uproot", "DR_below", "DR_up", "DR_branch"), 105)) |> filter(str_detect(fix_name, "CV")) |> 
  select(index, species, damage_type, fix_estimate, Std_Coefficient) |> left_join(bind_rows(data_set2_1) |> distinct(species_CN, species_LN), by = c("species" = "species_CN"))

sp_list <- sp_list_df(sp_list = unique(data$species_LN), taxon = "plant")
set.seed(1234); tree_plant <- get_tree(sp_list = sp_list, taxon = "plant", scenario = "random_below_basal", show_grafted = FALSE)
data.frame(x = tree_plant$tip.label, y = 1:35)

data <- data |> mutate(species_LN = case_when(species_LN == "Corylus_ferox_var._thibetica" ~ "Corylus_ferox", species_LN == "Quercus_aliena_var._acutiserrata" ~ "Quercus_aliena", species_LN == "Acer_pictum_subsp._mono" ~ "Acer_pictum", species_LN == "Cornus_kousa_subsp._chinensis" ~ "Cornus_kousa", TRUE ~ species_LN)) |> 
  left_join(data.frame(species_LN = c("Carpinus turczaninowii", "Carpinus fargesiana", "Corylus ferox", "Corylus chinensis", "Betula albosinensis", "Betula luminifera", "Juglans mandshurica", "Platycarya strobilacea", "Quercus serrata", "Quercus aliena", "Quercus variabilis", "Castanea seguinii", "Castanea henryi", "Fagus engleriana", "Sorbus folgneri", "Sorbus alnifolia", "Crataegus wilsonii", "Prunus conradinae", "Prunus padus", "Populus lasiocarpa", "Populus davidiana", "Salix wallichiana", "Rhus potaninii", "Rhus chinensis", "Toxicodendron vernicifluum", "Acer davidii", "Acer pictum", "Cornus macrophylla", "Cornus kousa", "Cornus controversa", "Fraxinus chinensis", "Diospyros lotus", "Litsea ichangensis", "Lindera obtusiloba", "Pinus armandii"), 
                       ID = c(35:1), niche_width = c(0.5693252, 0.5067005, 0.4052766, 0.8033091, 0.7300763, 0.6806774, 0.7154775, 0.7134962, 0.4579634, 0.5269582, 0.5795512, 0.6247113, 0.539175, 0.4583443, 0.4789034, 0.4000346, 0.5557265, 0.5170674, 0.6401926, 0.5663812, 0.6025009, 0.6959552, 0.7151405, 0.7137618, 0.4638151, 0.7080873, 0.6139691, 0.516735, 0.3869133, 0.6158935, 0.4840941, 0.6366028, 0.5825325, 0.4570297, 0.7790024)) |> mutate(species_LN = str_replace_all(species_LN, " ", "_"))) |> arrange(ID)

# phylogenetic heatmap --
plot_01 <- ggplot(tree_plant, branch.length = "none") + 
  geom_tree(linewidth = 0.25) + 
  geom_text(data = data |> distinct(species_LN, ID) |> mutate(species_LN = str_replace_all(species_LN, "_", " ")), aes(x = 11, y = ID, label = species_LN), size = (7*0.35), hjust = "left", fontface = "italic", family = "serif") + 
  geom_point(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_TTH"), aes(x = 30, y = ID, color = Std_Coefficient), shape = 15, size = (10*0.35), stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_TTH"), aes(x = 30, y = ID, label = ifelse(Std_Coefficient > 0, "+", "−")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_DBH"), aes(x = 33, y = ID, color = Std_Coefficient), shape = 15, size = (10*0.35), stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_DBH"), aes(x = 33, y = ID, label = ifelse(Std_Coefficient > 0, "+", "−")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_NMB"), aes(x = 36, y = ID, color = Std_Coefficient), shape = 15, size = (10*0.35), stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_NMB"), aes(x = 36, y = ID, label = ifelse(Std_Coefficient > 0, "+", "−")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 30.0, y = 0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 40, hjust = 1) + 
  annotate(geom = "text", x = 33.0, y = 0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 40, hjust = 1) + 
  annotate(geom = "text", x = 36.0, y = 0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 40, hjust = 1) + 
  annotate(geom = "text", x = 3, y = -0.5, label = "Slopes", size = (9*0.35), family = "serif") + 
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, "RdBu")[10:5])(99), breaks = c(-0.15, -0.05, 0.05), labels = c("−0.15", "−0.05", "0.05")) + 
  scale_x_continuous(limits = c(0, 38)) + 
  scale_y_continuous(limits = c(-1, 36)) + 
  theme_bw(base_family = "serif") + 
  labs(x = NULL, y = NULL, color = NULL) + 
  theme(legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.2, "cm"), 
        legend.direction = "horizontal", 
        legend.position = c(0.40, 0.056),
        legend.background = element_blank(),
        legend.text = element_text(size = 9, margin = margin(t = 1, unit = "pt")),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.margin = margin(t = -22, r = 0, b = -16, l = -9)); plot_01

# ____5.3.2 Fig. 3B --------------------------------------------------------------------------------
# data processing --
data <- read.xlsx("save/model_Part2_Q2_multi_data_20250703.xlsx")
data <- data |> filter(fix_name != "(Intercept)") |> filter(index %in% c("CV_TTH", "CV_DBH", "CV_NMB")) |> filter(Y == "DR_total") |> 
  mutate(fix_name = case_when(fix_name == "niche_width" ~ "1_NB", fix_name == "woody_density" ~ "2_WD", fix_name == "P50" ~ "3_P50", fix_name == "rdmax" ~ "4_RDmax", TRUE ~ "5_Hmax")) |> arrange(index, fix_name) |> 
  mutate(ID = c(14, 16, 18, 20, 22, 26, 28, 30, 32, 34, 2, 4, 6, 8, 10)) |> set_names("fix_name", "estimate", "se", "se_adjusted", "z_palue", "p_value", "Y", "index", "ID")

# forest plot --
plot_02 <- ggplot(data, aes(x = estimate, y = ID)) + 
  geom_rect(aes(xmin = -0.044, xmax = 0.044, ymin =  0.5, ymax = 11.5), fill = "#8cc269", alpha = 0.008) + 
  geom_rect(aes(xmin = -0.044, xmax = 0.044, ymin = 12.5, ymax = 23.5), fill = "#1ba784", alpha = 0.008) + 
  geom_rect(aes(xmin = -0.044, xmax = 0.044, ymin = 24.5, ymax = 35.5), fill = "#1a6840", alpha = 0.008) + 
  
  geom_errorbarh(aes(xmax = estimate + 1.96 * se, xmin = estimate - 1.96 * se, color = index), height = 0, linewidth = 3, alpha = 0.5) + 
  geom_point(aes(color = index), size = 2.5, alpha = 0.7, shape = 16) + 
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 35.5), linetype = 2, linewidth = 0.4) + 

  annotate(geom = "text", x = -0.040, y = 10.8, label = expression(CV[H]), size = (10*0.35), family = "serif", fontface = "bold", color = "#8cc269", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 2, label = "NB", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.033, y = 2, label = expression(""^"*"), size = (12*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 4, label = "WD", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 6, label = "P50", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 8, label = "RDmax", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 10, label = "Hmax", size = (9*0.35), family = "serif", hjust = "left") + 
  
  annotate(geom = "text", x = -0.040, y = 22.8, label = expression(CV[DBH]), size = (10*0.35), family = "serif", fontface = "bold", color = "#1ba784", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 14, label = "NB", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.033, y = 14, label = expression(""^"*"), size = (12*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 16, label = "WD", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 18, label = "P50", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 20, label = "RDmax", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 22, label = "Hmax", size = (9*0.35), family = "serif", hjust = "left") + 
  
  annotate(geom = "text", x = -0.040, y = 34.8, label = expression(CV[NMB]), size = (10*0.35), family = "serif", fontface = "bold", color = "#1a6840", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 26, label = "NB", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.033, y = 26, label = expression(""^"**"), size = (12*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 28, label = "WD", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 30, label = "P50", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 32, label = "RDmax", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 34, label = "Hmax", size = (9*0.35), family = "serif", hjust = "left") + 
  scale_x_continuous(limits = c(-0.044, 0.044), breaks = c(-0.04, 0, 0.04), labels = c("−0.04", "0", "0.04")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 36)) + 
  scale_fill_manual(values = c("#1ba784", "#1a6840", "#8cc269")) + 
  scale_color_manual(values = c("#1ba784", "#1a6840", "#8cc269")) + 
  labs(x = "Standardized effect sizes", y = NULL) + 
  theme_classic(base_family = "serif") + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        axis.line.y  = element_blank(), 
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = "none", 
        plot.margin = margin(t = 0, r = 2, b = 2, l = 0)); plot_02

# ____5.3.3 layout ---------------------------------------------------------------------------------
plot_grid(plot_01, plot_02, nrow = 1, rel_widths = c(0.52, 0.48)) + 
  annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.98), label = "a", size = (10*0.35), family = "serif", fontface = "bold") + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.98), label = "b", size = (10*0.35), family = "serif", fontface = "bold")

# ggsave(file = "save/fig_3_20250703.tiff", width = 12, height = 12, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")
rm(data, sp_list, model_save_rlm, tree_plant, plot_01, plot_02)

# __5.4 Figure 04 ----------------------------------------------------------------------------------
# ____5.4.1 Fig. 4A --------------------------------------------------------------------------------
# use PPT for drawing --

# ____5.4.2 Fig. 4B --------------------------------------------------------------------------------
# use PPT for drawing --

# ____5.4.3 Fig. 4C --------------------------------------------------------------------------------
# use PPT for drawing --

# ____5.4.4 Fig. 4D --------------------------------------------------------------------------------
data.frame(n = c(1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15), 
                      type = c("young", "young", "young", "young", "middle", "middle", "middle", "middle", "mature", "mature", "mature", "mature"), 
                      effect = c(-0.45, 0, -0.19, 0, -0.45, -0.09, -0.22, 0.16, 0.01, -0.17, -0.24, 0.14)) |> 
  ggplot(aes(x = n, y = effect, fill = type)) + 
  geom_col(position = "stack", width = 0.9, alpha = 0.5) + 
  geom_hline(yintercept = 0, linewidth = 0.2) + 
  
  geom_vline(xintercept = 4, linewidth = 0.3, linetype = 5, color = "#7f7f7f") + 
  geom_vline(xintercept = 8, linewidth = 0.3, linetype = 5, color = "#7f7f7f") + 
  geom_vline(xintercept = 12, linewidth = 0.3, linetype = 5, color = "#7f7f7f") + 

  geom_rect(aes(xmin = 12.5, xmax = 13.3, ymin = -0.20, ymax = -0.24), fill = "#c5e0b4") + 
  annotate(geom = "text", x = 13.6, y = -0.22, label = "ESF", size = (10*0.35), family = "serif", hjust = "left") + 
  geom_rect(aes(xmin = 12.5, xmax = 13.3, ymin = -0.28, ymax = -0.32), fill = "#8dd3c1") + 
  annotate(geom = "text", x = 13.6, y = -0.30, label = "LSF", size = (10*0.35), family = "serif", hjust = "left") + 
  geom_rect(aes(xmin = 12.5, xmax = 13.3, ymin = -0.36, ymax = -0.40), fill = "#8cb39f") + 
  annotate(geom = "text", x = 13.6, y = -0.38, label = "OF", size = (10*0.35), family = "serif", hjust = "left") + 
  
  scale_fill_manual(values = c("#1a6840", "#1ba784", "#8cc269"), labels = c("OF", "LSF", "ESF")) + 
  scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15), labels = c("", "Elevation", "", "", "Species\nrichness", "", "", "Structural\ndiversity", "", "", "Functional\nidentity", "")) + 
  scale_y_continuous(limits = c(-0.45, 0.2), breaks = c(-0.4, -0.2, 0, 0.2), labels = c("−0.4", "−0.2", "0", "0.2")) + 
  labs(x = NULL, y = "Total effects", fill = NULL) + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 9, color = "#000000"), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2))

# ggsave(file = "save/fig_4_20250703.tiff", width = 7.6, height = 6.4, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")

# __5.5 Figure S1 ----------------------------------------------------------------------------------
png_1 <- "data/icon/Uprooting.png"
png_2 <- "data/icon/Clear-bole broken.png"
png_3 <- "data/icon/Crown broken.png"
png_4 <- "data/icon/Branch broken.png"
load(file = "data/data_set1.rdata")

data_1 <- data_set1 |> filter(forest_age == "ESF") |> mutate(n_uproot = sum(n_uproot), n_below = sum(n_below), n_up = sum(n_up), n_branch = sum(n_branch), n_well = sum(n_total - n_damaged)) |> slice_head(n = 1) |> 
  select(n_uproot:n_branch, n_well) |> t() |> data.frame() |> rownames_to_column() |> set_names("damage_type", "n") |> mutate(ratio = round(n/sum(n)*100, 2)) |> select(!n)
data_1$damage_type <- c("1_Uprooting", "2_Clear-bole broken", "3_Crown broken", "4_Branch broken", "5_Undamaged"); data_1

data_2 <- data_set1 |> filter(forest_age == "LSF") |> mutate(n_uproot = sum(n_uproot), n_below = sum(n_below), n_up = sum(n_up), n_branch = sum(n_branch), n_well = sum(n_total - n_damaged)) |> slice_head(n = 1) |> 
  select(n_uproot:n_branch, n_well) |> t() |> data.frame() |> rownames_to_column() |> set_names("damage_type", "n") |> mutate(ratio = round(n/sum(n)*100, 2)) |> select(!n)
data_2$damage_type <- c("1_Uprooting", "2_Clear-bole broken", "3_Crown broken", "4_Branch broken", "5_Undamaged"); data_2

data_3 <- data_set1 |> filter(forest_age == "OF") |> mutate(n_uproot = sum(n_uproot), n_below = sum(n_below), n_up = sum(n_up), n_branch = sum(n_branch), n_well = sum(n_total - n_damaged)) |> slice_head(n = 1) |> 
  select(n_uproot:n_branch, n_well) |> t() |> data.frame() |> rownames_to_column() |> set_names("damage_type", "n") |> mutate(ratio = round(n/sum(n)*100, 2)) |> select(!n)
data_3$damage_type <- c("1_Uprooting", "2_Clear-bole broken", "3_Crown broken", "4_Branch broken", "5_Undamaged"); data_3

data <- bind_rows(data_1 |> mutate(type = "Early secondary-growth forest"), data_2 |> mutate(type = "Late secondary-growth forest"), data_3 |> mutate(type = "Old-growth forest")) |> mutate(ratio = ifelse(damage_type == "5_Undamaged", 100 - ratio, ratio))

ggplot(data, aes(x = damage_type, y = ratio, fill = type)) + 
  geom_col(position = "dodge", width = 0.8, alpha = 0.5) + 
  scale_x_discrete(labels = c("Uprooting", "Clear-bole\nbroken", "Crown\nbroken", "Branch\nbroken", "Total")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 27), breaks = c(5, 15, 25)) +
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  geom_image(aes(x = 1, y = 24), image = png_1, size = 0.19) + 
  geom_image(aes(x = 2, y = 24), image = png_2, size = 0.19) + 
  geom_image(aes(x = 3, y = 24), image = png_3, size = 0.19) + 
  geom_image(aes(x = 4, y = 24), image = png_4, size = 0.19) + 
  labs(x = NULL, y = "Proportion (%)", fill = "Forest stand categories") + 
  theme_bw(base_family = "serif") + 
  theme(panel.grid = element_blank(), 
        legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.5, "cm"), 
        legend.position = c(0.297, 0.6), 
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10))

# ggsave(file = "save/fig_S1_20250703.tiff", width = 10, height = 8, units = "cm", dpi = 300, limitsize = FALSE)
rm(data, data_1, data_2, data_3, png_1, png_2, png_3, png_4)

# __5.6 Figure S2 ----------------------------------------------------------------------------------
load(file = "data/data_set1.rdata")
data_set1 <- data_set1 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude))

# OF --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "OF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data1 <- effects::allEffects(model, xlevels = list(TD_SR = seq(1, 4.25, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF")
data1 <- data1 |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "OF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data2 <- effects::allEffects(model, xlevels = list(CV_TTH = seq(16.3, 62.3, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF")
data2 <- data2 |> mutate(CV_TTH = CV_TTH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "OF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data3 <- effects::allEffects(model, xlevels = list(CV_DBH = seq(29, 139, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF")
data3 <- data3 |> mutate(CV_DBH = CV_DBH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "OF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data4 <- effects::allEffects(model, xlevels = list(CV_NMB = seq(41.9, 184, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF")
data4 <- data4 |> mutate(CV_NMB = CV_NMB/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_01 <- ggplot() + 
  geom_line(data = data1, aes(x = TD_SR, y = fit, color = type, linetype = type), color = "#f43e06") + 
  geom_line(data = data2, aes(x = CV_TTH, y = fit, color = type, linetype = type), color = "#8cc269") + 
  geom_line(data = data3, aes(x = CV_DBH, y = fit, color = type, linetype = type), color = "#1ba784") + 
  geom_line(data = data4, aes(x = CV_NMB, y = fit, color = type, linetype = type), color = "#1a6840") + 
  scale_x_continuous(limits = c(0, 4.4), breaks = c(0, 2.2, 4.4)) +
  scale_y_continuous(limits = c(14, 40), breaks = c(14, 27, 40)) + 
  labs(x = "Diversity metrics", y = "Total damage ratio (%)", fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.94), label = "Old-growth forest", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.73), y = ggpp::as_npc(0.94), label = "Standardized slopes", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.87), label = expression(SR*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.87), label = expression("−0.166"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.80), label = expression(CV[H]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.80), label = expression("−0.175"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.73), label = expression(CV[DBH]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.73), label = expression("−0.195"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.66), label = expression(CV[NMB]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.66), label = expression("−0.201"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_01

# LSF --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "LSF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data1 <- effects::allEffects(model, xlevels = list(TD_SR = seq(1, 4.39, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF")
data1 <- data1 |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data2 <- effects::allEffects(model, xlevels = list(CV_TTH = seq(17.2, 56.8, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF")
data2 <- data2 |> mutate(CV_TTH = CV_TTH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data3 <- effects::allEffects(model, xlevels = list(CV_DBH = seq(31.5, 119, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF")
data3 <- data3 |> mutate(CV_DBH = CV_DBH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data4 <- effects::allEffects(model, xlevels = list(CV_NMB = seq(39.7, 252, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF")
data4 <- data4 |> mutate(CV_NMB = CV_NMB/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_02 <- ggplot() + 
  geom_line(data = data1, aes(x = TD_SR, y = fit, color = type, linetype = type), color = "#f43e06") + 
  geom_line(data = data2, aes(x = CV_TTH, y = fit, color = type, linetype = type), color = "#8cc269") + 
  geom_line(data = data3, aes(x = CV_DBH, y = fit, color = type, linetype = type), color = "#1ba784") + 
  geom_line(data = data4, aes(x = CV_NMB, y = fit, color = type, linetype = type), color = "#1a6840") + 
  scale_x_continuous(limits = c(0, 4.4), breaks = c(0, 2.2, 4.4)) +
  scale_y_continuous(limits = c(10, 40), breaks = c(10, 25, 40)) + 
  labs(x = "Diversity metrics", y = "Total damage ratio (%)", fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.94), label = "Late secondary-growth forest", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.73), y = ggpp::as_npc(0.94), label = "Standardized slopes", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.87), label = expression(SR*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.87), label = expression("−0.076"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.80), label = expression(CV[H]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.80), label = expression("−0.141"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.73), label = expression(CV[DBH]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.73), label = expression("−0.179"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.66), label = expression(CV[NMB]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.66), label = expression("−0.201"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_02

# ESF --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "ESF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data1 <- effects::allEffects(model, xlevels = list(TD_SR = seq(1.58, 4.75, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")
data1 <- data1 |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data2 <- effects::allEffects(model, xlevels = list(CV_TTH = seq(19.1, 55.8, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")
data2 <- data2 |> mutate(CV_TTH = CV_TTH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data3 <- effects::allEffects(model, xlevels = list(CV_DBH = seq(29.4, 126, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")
data3 <- data3 |> mutate(CV_DBH = CV_DBH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_set1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_set1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data4 <- effects::allEffects(model, xlevels = list(CV_NMB = seq(46.7, 248, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")
data4 <- data4 |> mutate(CV_NMB = CV_NMB/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_03 <- ggplot() + 
  geom_line(data = data1, aes(x = TD_SR, y = fit, color = type, linetype = type), color = "#f43e06") + 
  geom_line(data = data2, aes(x = CV_TTH, y = fit, color = type, linetype = type), color = "#8cc269") + 
  geom_line(data = data3, aes(x = CV_DBH, y = fit, color = type, linetype = type), color = "#1ba784") + 
  geom_line(data = data4, aes(x = CV_NMB, y = fit, color = type, linetype = type), color = "#1a6840") + 
  scale_x_continuous(limits = c(0, 5), breaks = c(0, 2.5, 5)) +
  scale_y_continuous(limits = c(9, 25), breaks = c(9, 17, 25)) + 
  labs(x = "Diversity metrics", y = "Total damage ratio (%)", fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.94), label = "Early secondary-growth forest", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.73), y = ggpp::as_npc(0.94), label = "Standardized slopes", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.87), label = expression(SR*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.87), label = expression(" 0.013"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.80), label = expression(CV[H]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.80), label = expression(" 0.000"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.73), label = expression(CV[DBH]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.73), label = expression(" 0.032"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.66), label = expression(CV[NMB]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.66), label = expression("−0.193"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_03

plot_grid(plot_01, plot_02, plot_03, nrow = 3) +  
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.990), label = "a", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.655), label = "b", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.322), label = "c", size = (10*0.35), family = "serif", fontface = "bold")

# ggsave(file = "save/fig_S2_20250703.tiff", width = 12, height = 18, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")
rm(plot_01, plot_02, plot_03, data1, data2, data3, data4, model)

# __5.7 Figure S3 ----------------------------------------------------------------------------------
# data processing --
load(file = "data/data_set2_1.rdata")
load(file = "save/model_Part2_Q1_data_20250703.rdata")

data <- model_save_rlm |> filter(!is.na(Std_Coefficient)) |> mutate(damage_type = rep(c("DR_total", "DR_uproot", "DR_below", "DR_up", "DR_branch"), 105)) |> filter(str_detect(fix_name, "CV")) |> 
  select(index, species, damage_type, fix_estimate, Std_Coefficient) |> left_join(bind_rows(data_set2_1) |> distinct(species_CN, species_LN), by = c("species" = "species_CN"))

sp_list <- sp_list_df(sp_list = unique(data$species_LN), taxon = "plant")
set.seed(1234); tree_plant <- get_tree(sp_list = sp_list, taxon = "plant", scenario = "random_below_basal", show_grafted = FALSE)
data.frame(x = tree_plant$tip.label, y = 1:35)

data <- data |> mutate(species_LN = case_when(species_LN == "Corylus_ferox_var._thibetica" ~ "Corylus_ferox", species_LN == "Quercus_aliena_var._acutiserrata" ~ "Quercus_aliena", species_LN == "Acer_pictum_subsp._mono" ~ "Acer_pictum", species_LN == "Cornus_kousa_subsp._chinensis" ~ "Cornus_kousa", TRUE ~ species_LN)) |> 
  left_join(data.frame(species_LN = c("Carpinus turczaninowii", "Carpinus fargesiana", "Corylus ferox", "Corylus chinensis", "Betula albosinensis", "Betula luminifera", "Juglans mandshurica", "Platycarya strobilacea", "Quercus serrata", "Quercus aliena", "Quercus variabilis", "Castanea seguinii", "Castanea henryi", "Fagus engleriana", "Sorbus folgneri", "Sorbus alnifolia", "Crataegus wilsonii", "Prunus conradinae", "Prunus padus", "Populus lasiocarpa", "Populus davidiana", "Salix wallichiana", "Rhus potaninii", "Rhus chinensis", "Toxicodendron vernicifluum", "Acer davidii", "Acer pictum", "Cornus macrophylla", "Cornus kousa", "Cornus controversa", "Fraxinus chinensis", "Diospyros lotus", "Litsea ichangensis", "Lindera obtusiloba", "Pinus armandii"), 
                       ID = c(35:1), niche_width = c(0.5693252, 0.5067005, 0.4052766, 0.8033091, 0.7300763, 0.6806774, 0.7154775, 0.7134962, 0.4579634, 0.5269582, 0.5795512, 0.6247113, 0.539175, 0.4583443, 0.4789034, 0.4000346, 0.5557265, 0.5170674, 0.6401926, 0.5663812, 0.6025009, 0.6959552, 0.7151405, 0.7137618, 0.4638151, 0.7080873, 0.6139691, 0.516735, 0.3869133, 0.6158935, 0.4840941, 0.6366028, 0.5825325, 0.4570297, 0.7790024)) |> mutate(species_LN = str_replace_all(species_LN, " ", "_"))) |> arrange(ID)
range(data$Std_Coefficient)

# phylogenetic heatmap --
ggplot(tree_plant, branch.length = "none") + 
  geom_tree() + 
  geom_text(data = data |> distinct(species_LN, ID) |> mutate(species_LN = str_replace_all(species_LN, "_", " ")), aes(x = 11, y = ID, label = species_LN), size = (7*0.35), hjust = "left", fontface = "italic", family = "serif") + 
  annotate(geom = "text", x = 1.7, y = 36.4, label = "Slopes", size = (9*0.35), family = "serif", fontface = "bold") + 
  
  geom_point(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_TTH"), aes(x = 22, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_TTH"), aes(x = 22, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_DBH"), aes(x = 24, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_DBH"), aes(x = 24, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_NMB"), aes(x = 26, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_NMB"), aes(x = 26, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 24, y = 36.8, label = "Uprooting", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 21, y = -0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 23, y = -0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 25, y = -0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 45) + 
  
  geom_point(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_TTH"), aes(x = 29, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_TTH"), aes(x = 29, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_DBH"), aes(x = 31, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_DBH"), aes(x = 31, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_NMB"), aes(x = 33, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_NMB"), aes(x = 33, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 31, y = 36.8, label = "Clear-bole", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 31, y = 36.0, label = "broken", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 28, y = -0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 30, y = -0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 32, y = -0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 45) + 
  
  geom_point(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_TTH"), aes(x = 36, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_TTH"), aes(x = 36, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_DBH"), aes(x = 38, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_DBH"), aes(x = 38, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_NMB"), aes(x = 40, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_NMB"), aes(x = 40, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 38, y = 36.8, label = "Crown", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 38, y = 36.0, label = "broken", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 35, y = -0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 37, y = -0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 39, y = -0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 45) + 
  
  geom_point(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_TTH"), aes(x = 43, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_TTH"), aes(x = 43, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_DBH"), aes(x = 45, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_DBH"), aes(x = 45, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_NMB"), aes(x = 47, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_NMB"), aes(x = 47, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 45, y = 36.8, label = "Branch", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 45, y = 36.0, label = "broken", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 42, y = -0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 44, y = -0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 46, y = -0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 45) + 
  
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, "RdBu")[10:2])(99), limits = c(-0.2, 0.42), breaks = c(-0.2, 0, 0.2, 0.4), labels = c("−0.2", " 0", " 0.2", " 0.4")) + 
  scale_y_continuous(limits = c(-2, 38)) + 
  theme_bw(base_family = "serif") + 
  labs(x = NULL, y = NULL, color = NULL) + 
  theme(legend.key.width = unit(0.2, "cm"), 
        legend.key.height = unit(0.5, "cm"), 
        legend.position = c(0.08, 0.80), 
        legend.background = element_blank(),
        legend.text = element_text(size = 9, margin = margin(l = 2, unit = "pt")),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.margin = margin(t = -25, r = -10, b = -25, l = -15))

# ggsave(file = "save/fig_S3_20250703.tiff", width = 14, height = 14, units = "cm", dpi = 300, limitsize = FALSE)
rm(data, model_save_rlm, sp_list, tree_plant)

# __5.8 Figure S4 ----------------------------------------------------------------------------------
# data processing --
load(file = "data/data_set1.rdata")

# DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ forest_age + (1 | site), data_set1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_01 <- ggplot(data_set1, aes(x = forest_age, y = DR_total*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Total damage", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.45) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "***", size = (12*0.35), family = "serif") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.81), xend = 2, yend = ggpp::as_npc(0.81)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 1.5, y = ggpp::as_npc(0.82), label = "***", size = (12*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_01

# DR_uproot --
model <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ forest_age + (1 | site), data_set1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_02 <- ggplot(data_set1, aes(x = forest_age, y = DR_uproot*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Uprooting", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.45) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "***", size = (12*0.35), family = "serif") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.81), xend = 2, yend = ggpp::as_npc(0.81)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 1.5, y = ggpp::as_npc(0.82), label = "**", size = (12*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_02

# DR_below --
model <- glmmTMB(cbind(n_below, n_total - n_below) ~ forest_age + (1 | site), data_set1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_03 <- ggplot(data_set1, aes(x = forest_age, y = DR_below*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Clear-bole broken", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "*", size = (12*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_03

# DR_up --
model <- glmmTMB(cbind(n_up, n_total - n_up) ~ forest_age + (1 | site), data_set1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_04 <- ggplot(data_set1, aes(x = forest_age, y = DR_up*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Crown broken", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "***", size = (12*0.35), family = "serif") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.81), xend = 2, yend = ggpp::as_npc(0.81)), color = "#4d4d4d", linewidth = 0.45) + 
  annotate(geom = "text", x = 1.5, y = ggpp::as_npc(0.82), label = "***", size = (12*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_04

# DR_branch --
model <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ forest_age + (1 | site), data_set1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_05 <- ggplot(data_set1, aes(x = forest_age, y = DR_branch*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Branch broken", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "**", size = (12*0.35), family = "serif") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.81), xend = 2, yend = ggpp::as_npc(0.81)), color = "#4d4d4d", linewidth = 0.45) + 
  annotate(geom = "text", x = 1.5, y = ggpp::as_npc(0.82), label = "***", size = (12*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_05

# TD_SR --
model <- glmmTMB(TD_SR ~ forest_age + (1 | site), data_set1)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_06 <- ggplot(data_set1, aes(x = forest_age, y = TD_SR, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.7, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  scale_y_continuous(limits = c(2, 32), breaks = c(2, 17, 32), labels = c("2", "17", "32")) + 
  labs(x = NULL, y = "SR") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "**", size = (12*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_06

plot_grid(plot_01, plot_02, plot_03, plot_04, plot_05, plot_06, nrow = 3) + 
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.990), label = "a", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.656), label = "c", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.323), label = "e", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.52), y = ggpp::as_npc(0.990), label = "b", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.52), y = ggpp::as_npc(0.656), label = "d", size = (10*0.35), family = "serif", fontface = "bold") + 
    annotate(geom = "text", x = ggpp::as_npc(0.52), y = ggpp::as_npc(0.323), label = "f", size = (10*0.35), family = "serif", fontface = "bold")

# ggsave(file = "save/fig_S4_20250703.tiff", width = 12, height = 18, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")
rm(model, plot_01, plot_02, plot_03, plot_04, plot_05, plot_06)

# __5.9 Figure S5 ----------------------------------------------------------------------------------
# use PPT for drawing --

# __5.0 Figure S6 ----------------------------------------------------------------------------------
ggplot() + scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) + scale_y_continuous(expand = c(0, 0), limits = c(0.9, 4.9)) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0), panel.grid = element_blank(), panel.background = element_rect(fill = "#ffffff")) + 
  
  geom_rect(aes(xmin = 0.5, xmax = 1.9, ymin = 4.1, ymax = 4.7), fill = "#7f7f7f", alpha = 0.5) + 
  geom_text(aes(x = 1.2, y = 4.4, label = "Elevation", family = "serif"), size = (9*0.35)) + 
  
  geom_rect(aes(xmin = 3.1, xmax = 4.5, ymin = 4.1, ymax = 4.7), fill = "#7f7f7f", alpha = 0.5) + 
  geom_text(aes(x = 3.8, y = 4.4, label = "Species\nrichness", family = "serif"), size = (9*0.35)) + 
  
  geom_rect(aes(xmin = 0.5, xmax = 1.9, ymin = 2.5, ymax = 3.1), fill = "#7f7f7f", alpha = 0.5) + 
  geom_polygon(aes(x = c(0.5, 0.7, 0.5), y = c(2.5, 2.5, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(0.5, 0.7, 0.5), y = c(2.8, 3.1, 3.1)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(1.7, 1.9, 1.9), y = c(2.5, 2.5, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(1.7, 1.9, 1.9), y = c(3.1, 3.1, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_text(aes(x = 1.2, y = 2.8, label = "Structural\ndiversity", family = "serif"), size = (9*0.35)) + 
  
  geom_rect(aes(xmin = 3.1, xmax = 4.5, ymin = 2.5, ymax = 3.1), fill = "#7f7f7f", alpha = 0.5) + 
  geom_polygon(aes(x = c(3.1, 3.3, 3.1), y = c(2.5, 2.5, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(3.1, 3.3, 3.1), y = c(2.8, 3.1, 3.1)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(4.3, 4.5, 4.5), y = c(2.5, 2.5, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(4.3, 4.5, 4.5), y = c(3.1, 3.1, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_text(aes(x = 3.8, y = 2.8, label = "Functional\nidentity", family = "serif"), size = (9*0.35)) + 
  
  geom_rect(aes(xmin = 1.8, xmax = 3.2, ymin = 0.9, ymax = 1.5), fill = "#7f7f7f", alpha = 0.5) + 
  geom_text(aes(x = 2.5, y = 1.2, label = "Damage\nratio", family = "serif"), size = (9*0.35)) + 
  
  # 01 line
  geom_segment(aes(x = 0.2, y = 4.4, xend = 0.5, yend = 4.4), linewidth = 0.2) + 
  geom_segment(aes(x = 0.2, y = 4.4, xend = 0.2, yend = 1.2), linewidth = 0.2) + 
  geom_segment(aes(x = 0.2, y = 1.2, xend = 1.8, yend = 1.2), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 0.2, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 0.2, y = 3.75, label = expression(italic("1")), size = (9*0.35), family = "serif") + 
  
  # 02 line
  geom_segment(aes(x = 1.9, y = 4.4, xend = 3.1, yend = 4.4), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 2.5, y0 = 4.4, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 2.5, y = 4.4, label = expression(italic("2")), size = (9*0.35), family = "serif") + 
  
  # 03 line
  geom_segment(aes(x = 1.2, y = 4.1, xend = 1.2, yend = 3.1), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 1.2, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 1.2, y = 3.75, label = expression(italic("3")), size = (9*0.35), family = "serif") + 
  
  # 04 line
  geom_segment(aes(x = 1.2, y = 4.1, xend = 3.8, yend = 3.1), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 2.11, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 2.11, y = 3.75, label = expression(italic("4")), size = (9*0.35), family = "serif") + 
  
  # 05 line
  geom_segment(aes(x = 3.8, y = 4.1, xend = 1.2, yend = 3.1), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 2.89, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 2.89, y = 3.75, label = expression(italic("5")), size = (9*0.35), family = "serif") + 
  
  # 06 line
  geom_segment(aes(x = 3.8, y = 4.1, xend = 3.8, yend = 3.1), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 3.8, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 3.8, y = 3.75, label = expression(italic("6")), size = (9*0.35), family = "serif") + 
  
  # 07 line
  geom_segment(aes(x = 4.8, y = 4.4, xend = 4.5, yend = 4.4), linewidth = 0.2) + 
  geom_segment(aes(x = 4.8, y = 4.4, xend = 4.8, yend = 1.2), linewidth = 0.2) + 
  geom_segment(aes(x = 4.8, y = 1.2, xend = 3.2, yend = 1.2), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 4.8, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 4.8, y = 3.75, label = expression(italic("7")), size = (9*0.35), family = "serif") + 
  
  # 08 line
  geom_segment(aes(x = 1.2, y = 2.5, xend = 2.5, yend = 1.5), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 1.85, y0 = 2.0, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 1.85, y = 2.0, label = expression(italic("8")), size = (9*0.35), family = "serif") + 
  
  # 09 line
  geom_segment(aes(x = 3.8, y = 2.5, xend = 2.5, yend = 1.5), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 3.15, y0 = 2.0, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 3.15, y = 2.0, label = expression(italic("9")), size = (9*0.35), family = "serif") + 
  
  # 10 line
  geom_curve(aes(x = 1.9, y = 2.8, xend = 3.1, yend = 2.8), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.25, curvature = 0.3) + 
  geom_curve(aes(x = 3.1, y = 2.8, xend = 1.9, yend = 2.8), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.25, curvature = -0.3) + 
  ggforce::geom_circle(aes(x0 = 2.5, y0 = 2.6, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 2.5, y = 2.6, label = expression(italic("10")), size = (9*0.35), family = "serif")

# ggsave(file = "save/fig_S6_20250703.tiff", width = 8, height = 6, units = "cm", dpi = 300, limitsize = FALSE)
rm(list = ls())
