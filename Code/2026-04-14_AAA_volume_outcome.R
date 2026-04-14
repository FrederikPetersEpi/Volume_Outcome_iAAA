# =============================================================================
# BAA Volume-Outcome Analysis 
# =============================================================================
# Purpose:   Analyze hospital volume vs. outcomes for open AAA repair (DIGG data).
# Author:    Frederik Peters
# Date:      2026-04-06
# =============================================================================

rm(list = ls())
gc()

# =============================================================================
# 1. LIBRARIES 
# =============================================================================
library(data.table)
library(mgcv)
library(marginaleffects)
library(ggplot2)
library(tableone)
library(patchwork)
library(flextable)
library(scales)
library(forester)
library(pROC)
library(caret)

# =============================================================================
# 2. LOAD & PREPARE DATA
# =============================================================================
setwd("C:/DIGG/baa/Data/")

# --- Facility-level info (exported once) ------------------------------------
dat.export = data.table(readRDS("data_gesamt_annotated.rds"))
dat.export = dat.export[jahr > 2017 & geschlecht != 3 & patientverstorben != 1]

# Postcode → Kreis-ID (two-step for intra-city regions)
dat.export[plz == "00000", plz := NA]
dat.export[, plz := na.omit(plz)[1], by = facilityid]
dat.export = dat.export[!is.na(plz)]

zuordnungkreis = fread("plz_to_kreis_id.csv", encoding = "Latin-1")
zuordnungkreis[,n_char_plz:=nchar(Postleitzahl)]
zuordnungkreis[n_char_plz==5,plz:=as.character(Postleitzahl)]
zuordnungkreis[n_char_plz==4,plz:=as.character(paste0(0,Postleitzahl))]
zuordnungkreis[,Postleitzahl:=NULL]
zuordnungkreis[,n_char_plz:=NULL]
names(zuordnungkreis)=c("kreis_id","plz")

setkey(dat.export, plz)
setkey(zuordnungkreis, plz)
dat.export = zuordnungkreis[dat.export]

# Second pass for remaining postcodes
remaining = dat.export[is.na(kreis_id)]
zuordnungkreis2 = fread("plz_to_kreis_id_2.csv", encoding = "Latin-1")
zuordnungkreis2[,n_char_plz:=nchar(Postleitzahl)]
zuordnungkreis2[n_char_plz==5,plz:=as.character(Postleitzahl)]
zuordnungkreis2[n_char_plz==4,plz:=as.character(paste0(0,Postleitzahl))]
zuordnungkreis2[,Postleitzahl:=NULL]
zuordnungkreis2[,n_char_plz:=NULL]
names(zuordnungkreis2)=c("kreis_id","plz")
setkey(remaining, plz)
setkey(zuordnungkreis2, plz)
remaining = zuordnungkreis2[remaining]
dat.export = rbind(dat.export[!is.na(kreis_id)], remaining, fill = TRUE)

dat.export[kreis_id == 5118, kreis_id := 5119]
dat.export[kreis_id == 5912, kreis_id := 5911]

dat.export[, uhc := ifelse(krankenhaus_art == "Unive", 1L, 0L)]
dat.export[is.na(uhc),uhc:=0L]
dat.export=unique(dat.export[,.(gefaess_leistung_vc71=max(gefaess_leistung_vc71,na.rm=T),
                           gefaess_leistung_vc67=max(gefaess_leistung_vc67,na.rm=T),
                           anzahl_betten=max(anzahl_betten,na.rm=T),
                           uhc=max(uhc,na.rm=T),
                           open_n=max(open_n,na.rm=T),
                           endo_n=max(endo_n,na.rm=T)),.(facilityid,kreis_id,ik_nummer)])

# --- Main patient-level data ------------------------------------------------
dat = data.table(readRDS("data_gesamt.rds")) # 47821
dat = dat[geschlecht == 3, sex:=NA] # 47821
dat = dat[jahr > 2017] # 47821-30906
dat = dat[patientverstorben != 1] # 30906-30780
dat[,ruptured:=1,]
dat[is.na(klinik_ruptur)|klinik_ruptur==2,ruptured:=0] # if NA or 1 iAAA
dat = dat[ruptured==0] # 30780-28199
dat = dat[opverfahren != 3] # 28199-27770

setkey(dat.export, facilityid)
setkey(dat, facilityid)
dat = dat.export[dat]
dat = dat[!is.na(kreis_id)] # 27770-27095

# Certification
keyzert = fread("zertifizierung.csv")
names(keyzert)[1] = "ik_nummer"
keyzert = unique(keyzert)[,.SD[1], by = ik_nummer][, c(1, 2, 8:14)]
dat[,ik_nummer:=as.integer(ik_nummer)]
setkey(dat, ik_nummer)
setkey(keyzert, ik_nummer)
dat = keyzert[dat]
dat = dat[FALSCH != "x"] # 27095-24599

dat[, zertifizierung := 0L]
for (y in 2018:2024) {
  col = paste0("zert_", y)
  if (col %in% names(dat)) dat[jahr == y & get(col) == "x", zertifizierung := 1L]
}


# =============================================================================
# 3. DERIVE VOLUME & PATIENT VARIABLES 
# =============================================================================
dat[, `:=`(open_real = sum(opverfahren == 1L, na.rm = TRUE),
           endo_real = sum(opverfahren == 2L, na.rm = TRUE)),
    by = .(facilityid, jahr)]

dat[, `:=`(open_real_max = max(open_real, na.rm = TRUE),
           endo_real_max = max(endo_real, na.rm = TRUE)),
    by = facilityid]

# BMI
dat[groesse < 150 | groesse > 210, groesse := NA_real_]
dat[gewicht < 45 | gewicht > 160, gewicht := NA_real_]
dat[, bmi := gewicht / (0.01 * groesse)^2]


# Shorten variable names
setnames(dat, gsub("risikokonstellation_comm_", "", names(dat)))

# Categorical variables 
dat[, `:=`(
  koronareherzerkrankung = fcase(koronareherzerkrankung == 0L, "no",
                                 koronareherzerkrankung == 1L, "yes",
                                 default = NA_character_),
  nyha_cat = fcase(herzinsuffizienznachnyha %in% 1:2, "NYHA I–II",
                   herzinsuffizienznachnyha %in% 3:4, "NYHA III–IV",
                   default = "none"),
  copd_cat = fcase(copdcold %in% 1:2, "GOLD I–II",
                   copdcold %in% 3:4, "GOLD III–IV",
                   default = "none"),
  asa_cat = fcase(asaklassifikation %in% 1:2, "ASA I–II",
                  asaklassifikation == 3L, "ASA III",
                  asaklassifikation %in% 4:5, "ASA IV–V",
                  default = NA_character_),
  asa_3_to5 = fifelse(asaklassifikation %in% 3:5, "no", "yes"),
  prior_aorta_surgery = fifelse(vorhergehendeaortenchirurgieoderaortenstent==1,"yes","no"),
  diabetes = fifelse(diabetes==1,"yes","no"),
  hypertonus = fifelse(hypertonus==1,"yes","no"),
  dyslipoprot = fifelse(dyslipoprot==1,"yes","no"),
  nikotinakt = fifelse(nikotinakt==1,"yes","no"),
  aneurysmaand = fifelse(aneurysmaand==1,"yes","no"),
  renalebegleiterkrankungen = fifelse(renalebegleiterkrankungen==1,"yes","no"),
  pulmonalebegleiterkrankungen = fifelse(pulmonalebegleiterkrankungen==1,"yes","no"),
  kardialebegleiterkrankungen = fifelse(kardialebegleiterkrankungen==1,"yes","no"),
  lokalisation = fifelse(lokalisation==1,"juxtarenal","infrarenal"),
  iliacalebeteiligung = fifelse(iliacalebeteiligung==1,"yes","no"),
  age_cat = cut(alterbeiop, breaks = c(-Inf, 64, 74, 84, Inf),
                labels = c("<65", "65–74", "75–84", ">=85"), right = TRUE),
  bmi_cat = cut(bmi, breaks = c(-Inf, 18.5, 25, 30, 35, Inf),
                labels = c("<18.5", "18.5–24.9", "25–29.9", "30–34.9", ">=35"),
                right = FALSE),
  diam_cat = cut(maximalerdurchmesserdesaneurysmas,
                 breaks = c(-Inf, 49, 54, 59, 69, Inf),
                 labels = c("<50", "50–54", "55–59", "60–69", ">=70"),
                 right = TRUE),
  fontaine_cat = fcase(stadium == 0L, "0 (keine AVK)",
                       stadium %in% 1:2, "I–IIa",
                       stadium == 3L, "IIb",
                       stadium %in% 4:5, "III–IV",
                       default = NA_character_)
)]

factor_cols = c("aneurysmaand", "nikotinakt", "dyslipoprot", "hypertonus",
                "diabetes", "lokalisation", "iliacalebeteiligung",
                "renalebegleiterkrankungen", "pulmonalebegleiterkrankungen",
                "kardialebegleiterkrankungen", "asa_cat", "fontaine_cat",
                "asa_3_to5", "prior_aorta_surgery",
                "koronareherzerkrankung", "nyha_cat", "copd_cat")
dat[, (factor_cols) := lapply(.SD, as.factor), .SDcols = factor_cols]

# Assign unknown to Fontaine stage and prior Surgery
dat[is.na(fontaine_cat), fontaine_cat := "unknown"]
dat[is.na(prior_aorta_surgery), prior_aorta_surgery := "unknown"]
dat[is.na(bmi_cat), bmi_cat := "unknown"]
#dat[is.na(asa_cat), asa_cat := "unknown"]
dat[is.na(renalebegleiterkrankungen), renalebegleiterkrankungen := "no"]
dat[is.na(pulmonalebegleiterkrankungen), pulmonalebegleiterkrankungen := "no"]
dat[is.na(kardialebegleiterkrankungen), kardialebegleiterkrankungen := "no"]

# Binary outcomes (
dat[, `:=`(
  tod                = fifelse(patientpostoperativverstorben == 1L | patientverstorben == 1L | !is.na(sterbedatum), 1L, 0L, na = 0L),
  major_complication = fifelse(neuaufgetreteneherzinsuffizienz == 1L | neuerherzinfarkt == 1L | neuepneumonie == 1L | kompembol == 1L | neuerespiratorischeinsuffizienz == 1L | neuerapoplex == 1L | sepsis == 1L | sirs == 1L | neuerenaleinsuffizienz == 1L | darmischaemie == 1L | abdominell_nachblutung == 1L | platzbauch == 1L | dialyse == 1L | reintubat == 1L | langzeitbeatmung == 1L | unterstuetzverf == 1L | zweiteingriff == 1L | vasculaererevisionzf == 1L | laparatomie == 1L | haematomausraeumung == 1L, 1L, 0L, na = 0L),
  revision           = fifelse(mitrevisionsbeduerftigenlokalenkomplikationenabdominell == 1L | mitrevisionsbeduerftigenlokalenkomplikationeninguinal == 1L, 1L, 0L, na = 0L)
)]

# FTR if major complication and death
dat[major_complication == 1L & tod == 0L, ftr := 0L]
dat[major_complication == 1L & tod == 1L, ftr := 1L]
table(dat$ftr,useNA="always") # 2212 and 343

# Composite endpoint of major complications and death for in-depth analysis of factors
dat[ , compo := 0L]
dat[major_complication == 1L | tod == 1L, compo := 1L]
table(dat$compo,useNA="always") # 2212 and 343

# =============================================================================
# 4. PROPENSITY SCORE MODEL (Open vs EVAR)
# =============================================================================
dat = dat[opverfahren %in% c(1L, 2L)]
dat[, treatment := fifelse(opverfahren == 1L, 1L, 0L)]

ps_vars = c("jahr", "age_cat", "diam_cat", "asa_cat", "geschlecht", "bmi_cat",
            "fontaine_cat", "copd_cat", "nyha_cat", "koronareherzerkrankung",
            "prior_aorta_surgery", "diabetes", "hypertonus", "dyslipoprot",
            "nikotinakt", "aneurysmaand", "renalebegleiterkrankungen",
            "pulmonalebegleiterkrankungen", "kardialebegleiterkrankungen",
            "lokalisation", "iliacalebeteiligung","geschlecht")
dat[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = ps_vars]
dat = na.omit(dat, cols = ps_vars) # 24599-23671

ps_model = glm(treatment ~ jahr + age_cat + diam_cat + asa_cat + geschlecht +
                 bmi_cat + fontaine_cat + copd_cat + nyha_cat +
                 koronareherzerkrankung + prior_aorta_surgery +
                 diabetes + hypertonus + dyslipoprot + nikotinakt + aneurysmaand +
                 renalebegleiterkrankungen + pulmonalebegleiterkrankungen +
                 kardialebegleiterkrankungen + lokalisation + iliacalebeteiligung,
                 data = dat, family = binomial(link = "logit"))

# AUC value
library(pROC)
dat$probs = predict(ps_model, type = "response")
roc_score = roc(dat$treatment, dat$probs)
print(auc(roc_score))
nrow(dat)

# PS table 
sum_ps = summary(ps_model)
est = sum_ps$coefficients[-1, 1]
se = sum_ps$coefficients[-1, 2]
or_tab = data.table(
  Variable = names(est),
  OR_95CI = sprintf("%.2f (%.2f-%.2f)", round(exp(est), 2), round(exp(est - 1.96 * se), 2), round(exp(est + 1.96 * se), 2)),
  Importance = round(varImp(ps_model, scale = FALSE)$Overall, 2)
)
or_tab$Variable = c(
  "Year of surgery", "Age 65–74 years (ref: <65)", "Age 75–84 years", "Age ≥85 years",
  "Aneurysm diameter 50–54 mm (ref: <50)", "Aneurysm diameter 55–59 mm",
  "Aneurysm diameter 60–69 mm", "Aneurysm diameter ≥70 mm",
  "ASA I–II (ref: III)", "ASA IV–V", "Female sex (ref: male)",
  "BMI 18.5–24.9 kg/m² (ref: <18.5)", "BMI 25–29.9 kg/m²",
  "BMI 30–34.9 kg/m²", "BMI ≥35 kg/m²", "BMI unknown" ,
  "Fontaine I–IIa (ref: none)", "Fontaine IIb", "Fontaine III–IV", "Fontaine unknown",
  "COPD GOLD III–IV (ref: I–II)", "No COPD",
  "NYHA I–II (ref: none)", "NYHA III–IV",
  "Coronary heart disease", "Prior aortic surgery",
  "Prior aortic surgery unknown", "Diabetes mellitus", "Arterial hypertension",
  "Dyslipidaemia", "Current smoking", "Other aneurysm locations",
  "Renal comorbidities)", "Pulmonary comorbidities",
  "Cardiac comorbidities", "Location (ref.: juxtarenal)",
  "Iliac involvement"
)
ft_ps = flextable(or_tab[, .(Variable, OR_95CI, Importance)])
ft_ps = font(ft_ps, fontname = "Arial", part = "all")
ft_ps = fontsize(ft_ps, size = 8, part = "all")
ft_ps = autofit(ft_ps)
save_as_docx(ft_ps, path = paste0(Sys.Date(), "_table_ps_model.docx"))

# Patients
uniqueN(dat)

# Centres
uniqueN(dat$facilityid)

dat[, propensity_score := predict(ps_model, type = "response")]

# PS density plot
dat[, group := fifelse(opverfahren == 1L, "Open", "EVAR")]

label_dt = dat[, .N, by = group]
label_dt[, label := paste0(group, "\n(n = ", formatC(N, big.mark = ","), ")")]
label_map = label_dt[, setNames(label, group)]
dat[, group_label := label_map[group]]
dat[, group_label := factor(group_label, levels = label_dt$label)]

p_ps = ggplot(dat, aes(x = propensity_score, fill = group_label, colour = group_label)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  scale_fill_manual(values = c("#E63946", "#457B9D")) +
  scale_colour_manual(values = c("#E63946", "#457B9D")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), labels = percent_format()) +
  labs(x = "Estimated Probability of Open Surgical Repair (Propensity Score)",
       y = "Density", fill = NULL, colour = NULL) +
  theme_minimal() + 
  theme(legend.position = c(0.88, 0.88))
print(p_ps)

ggsave("ps_distribution.png", plot = p_ps, width = 10, height = 6, dpi = 300, bg = "white")

# =============================================================================
# 5. OPEN SURGERY CASES ONLY
# =============================================================================
dt.open = dat[opverfahren == 1L] # 23670-4796
rm(dat)

dt.open[open_real_max > 40L, open_real_max := 40L]
dt.open[open_real > 40L, open_real := 40L]
dt.open[open_n > 40L, open_n := 40L]

summary(dt.open$endo_real_max)
summary(dt.open$endo_real)
summary(dt.open$endo_n)

dt.open[endo_real_max > 40L, endo_real_max := 80L]
dt.open[endo_real > 40L, endo_real := 80L]
dt.open[endo_n > 40L, endo_n := 80L]

dt.open[, share_open := 100 * (open_n / (open_n + endo_n))]
dt.open[, `:=`(high_10 = fifelse(open_real_max >= 10L, 1L, 0L),
               high_15 = fifelse(open_real_max >= 15L, 1L, 0L))]

ps_vars = c("open_real_max", "open_real", "open_n", "endo_real_max", "endo_real", "endo_n",
            "share_open", "high_10", "high_15", "tod","ftr", "major_complication", "fontaine_cat",
            "bmi_cat", "asa_cat","nyha_cat","propensity_score",
            "age_cat", "renalebegleiterkrankungen")
dt.open[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = ps_vars]

# =============================================================================
# 6. MAIN VOLUME-OUTCOME MODELS (GAM)
# =============================================================================
outcome_labels = c(tod = "Death", ftr = "Failure to Rescue",
                   major_complication = "Major Complication", compo = "Major Complication or death")
exposure_labels = c(open_real = "Annual volume, DGG registry",
                    open_real_max = "Highest annual volume, DGG registry",
                    open_n = "Annual volume, GBA quality report",
                    high_10 = "Highest annual volume ≥10, DGG registry",
                    high_15 = "Highest annual volume ≥15, DGG registry")

exposures = c("open_real", "open_real_max", "open_n", "high_10", "high_15")
outcomes = c("tod", "ftr", "major_complication", "compo")

results_long = rbindlist(lapply(outcomes, function(out) {
  rbindlist(lapply(exposures, function(exp) {
    term = ifelse(exp %in% c("open_real", "open_real_max", "open_n"),
                  paste0("scale(", exp, ")"), exp)
    
    frm = as.formula(paste0(
      out, " ~ ", term, " + fontaine_cat + bmi_cat + asa_cat + age_cat + ",
      "renalebegleiterkrankungen + nyha_cat + scale(propensity_score) + ",
      "s(facilityid, bs = 're')"))
    
    fit = gam(frm, data = dt.open, family = binomial(link = "logit"))
    summ = summary(fit)
    coef_name = ifelse(exp %in% c("open_real", "open_real_max", "open_n"),
                       paste0("scale(", exp, ")"), exp)
    idx = which(names(summ$p.coeff) == coef_name)
    
    est = summ$p.coeff[idx]
    se = summ$se[idx]
    OR_str = sprintf("%.2f (%.2f-%.2f)", exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))
    
    data.table(outcome = out, exposure = exp, OR_str = OR_str)
  }))
}))

results_long[, `:=`(outcome_label = outcome_labels[outcome],
                    exposure_label = exposure_labels[exposure])]
results_wide = dcast(results_long, exposure_label ~ outcome_label, value.var = "OR_str")
setnames(results_wide, "exposure_label", "Exposure")

ft_results = flextable(results_wide)
ft_results = font(ft_results, fontname = "Arial", part = "all")
ft_results = fontsize(ft_results, size = 8, part = "all")
ft_results = autofit(ft_results);ft_results

save_as_docx(ft_results, path = paste0(Sys.Date(), "_table2_outcomes.docx"))

# =============================================================================
# 7. VOLUME-OUTCOME PLOTS
# =============================================================================
plots = lapply(outcomes, function(out) {
  dt.open.plot = dt.open[!is.na(get(out))]
  mean_out = mean(dt.open.plot[[out]], na.rm = TRUE)
  ylim_up = ceiling(mean_out * 3 * 10) / 10
  
  frm = as.formula(paste0(
    out, " ~ s(open_real_max, bs = 'ps') + fontaine_cat + bmi_cat + asa_cat + ",
    "age_cat + renalebegleiterkrankungen + nyha_cat + scale(propensity_score) + ",
    "s(facilityid, bs = 're')"))
  
  fit = gam(frm, data = dt.open.plot, family = binomial(link = "logit"))
  print(summary(fit))
  center_data = dt.open.plot[!is.na(get(out)) & !is.na(open_real_max),
                        .(n = .N, events = sum(get(out)), rate = mean(get(out))),
                        by = .(facilityid, open_real_max)]
  center_data = center_data[rate <= ylim_up]
  center_data[, vol_group := fcase(open_real_max <= 9L, "0-9",
                                   open_real_max <= 14L, "10-14",
                                   default = "15+")]
  
  pred_plot = plot_predictions(fit, condition = "open_real_max", type = "response") +
    labs(y = "Probability", x = "Annual Volume", title = outcome_labels[out], caption = "") +
    theme_minimal(base_size = 9) +
    scale_y_continuous(limits = c(0, ylim_up), labels = number_format(accuracy = 0.01),
                       breaks = pretty_breaks(n = )) +
    theme(plot.title = element_text(face = "bold", size = 9))
  
  pred_plot +
    geom_point(data = center_data, aes(x = open_real_max, y = rate, size = n, color = vol_group),
               alpha = 0.3, stroke = 0) +
    scale_size_continuous(range = c(1.5, 6), guide = "none") +
    scale_color_manual(values = c("0-9" = "#E31A1C", "10-14" = "#FDAF6C", "15+" = "#1F78B4"),
                       name = "Annual volume", guide = guide_legend(override.aes = list(size = 3))) +
    theme(legend.position = "bottom")
})

wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(paste0(Sys.Date(), "_volume_outcome_plots.png"),
       width = 11, height = 10, dpi = 600, bg = "white", units = "in", device = "png")

# =============================================================================
# 8a. FOREST PLOT – Covariate effects on Composite Endpoint
# =============================================================================
dt.open[, (c("renalebegleiterkrankungen", "pulmonalebegleiterkrankungen",
             "kardialebegleiterkrankungen", "iliacalebeteiligung")) :=
          lapply(.SD, function(x) relevel(factor(x), ref = "no")),
        .SDcols = c("renalebegleiterkrankungen", "pulmonalebegleiterkrankungen",
                    "kardialebegleiterkrankungen", "iliacalebeteiligung")]

dt.open[, high_15 := factor(high_15, levels = c(0L, 1L), labels = c("0-14 cases", "≥15 cases"))]
dt.open[, sex := factor(geschlecht, levels = c(1L, 2L), labels = c("men", "women"))]

dt.open[, age_cat := relevel(factor(age_cat), ref = "65–74")]
dt.open[, asa_cat := relevel(factor(asa_cat), ref = "ASA III")]
dt.open[, bmi_cat := relevel(factor(bmi_cat), ref = "18.5–24.9")]
dt.open[, renalebegleiterkrankungen := relevel(factor(renalebegleiterkrankungen), ref = "no")]
dt.open[, propensity_score := scale(propensity_score)]

model_ftr = gam(compo ~ high_15 + fontaine_cat + bmi_cat + asa_cat + age_cat +
                  nyha_cat + renalebegleiterkrankungen + propensity_score +
                  s(facilityid, bs = "re"),
                data = dt.open, family = binomial(link = "logit"))

sum_ftr = summary(model_ftr)
est_ftr = sum_ftr$p.coeff[-1]
se_ftr = sum_ftr$se[c(-1,-21)]

forest_tab = data.table(
  variable = names(est_ftr),
  Estimates = round(exp(est_ftr), 3),
  LL = round(exp(est_ftr - 1.96 * se_ftr), 3),
  UL = round(exp(est_ftr + 1.96 * se_ftr), 3)
)

cat_vars = c("high_15", "fontaine_cat", "bmi_cat", "asa_cat", "age_cat",
             "renalebegleiterkrankungen", "nyha_cat")
res_cat = rbindlist(lapply(cat_vars, function(v) {
  dt.open[!is.na(get(v)), .(total = .N, level = get(v)), by = v][
    order(match(level, levels(dt.open[[v]])))
  ][, .(variable = paste0(v, level), variable2 = v, variable3 = level, total)]
}))

res_num = data.table(variable = "propensity_score", variable2 = "propensity_score",
                     variable3 = "propensity_score", total = nrow(dt.open))

res_display = rbind(res_cat, res_num)

res_display[, row_nr := .I]

setkey(res_display,variable)
setkey(forest_tab,variable)
forest_tab = forest_tab[res_display,][order(row_nr)]

forest_tab[is.na(Estimates), variable := variable2]
forest_tab[!is.na(Estimates), variable := paste0("   ", variable3)]
forest_tab[is.na(Estimates),variable:=variable2,]
forest_tab[!is.na(Estimates),variable:=paste0("   ",variable3),]


forest_tab[variable == "high_15", variable := "Highest annual volume, ref.: 0-14 cases"]
forest_tab[variable == "fontaine_cat", variable := "Fontaine stage, ref.: none"]
forest_tab[variable == "bmi_cat", variable := "BMI category, kg/m², ref.: 18.5–24.9"]
forest_tab[variable == "asa_cat", variable := "ASA classification, ref.: ASA III"]
forest_tab[variable == "age_cat", variable := "Age group, years, ref.: 65–74"]
forest_tab[variable == "nyha_cat", variable := "NYHA class, ref.: none"]
forest_tab[variable == "renalebegleiterkrankungen", variable := "Renal comorbidities, ref.: none"]
forest_tab[variable == "   propensity_score", variable := "Propensity open repair, change 1 SD"]

forest_tab[, `:=`(Variable = variable, N = total)]

forest_plot_os = forester(
  left_side_data = forest_tab[, .(Variable, N)],
  estimate = forest_tab$Estimates,
  ci_low = forest_tab$LL,
  ci_high = forest_tab$UL,
  display = TRUE, estimate_precision = 2, file_path = "figure_2.png", dpi = 600,
  xlim = c(0.1, 10), null_line_at = 1, font_family = "Fira Sans",
  x_scale_linear = FALSE, arrows = TRUE,
  arrow_labels = c("lower", "higher"),
  estimate_col_name = "Odds Ratio (95% CI)", point_sizes = 3
)

# ==================================================================================
# 8b. FOREST PLOT – Covariate effects and structural variables on Composite Endpoint
# ==================================================================================
dt.open.b = copy(dt.open)

# Convert 0/1 binary variables to "no"/"yes" factors immediately after copying
binary_vars = c("uhc", "gefaess_leistung_vc71", "gefaess_leistung_vc67", "zertifizierung")
for(v in binary_vars) {
  if(v %in% names(dt.open.b)) {
    dt.open.b[, (v) := factor(get(v), levels = c(0, 1), labels = c("no", "yes"))]
  }
}

dt.open.b[, share_open := as.numeric(scale(share_open))]
dt.open.b[, endo_real_max := as.numeric(scale(endo_real_max))]
dt.open.b[, anzahl_betten := as.numeric(scale(anzahl_betten))]

model_ftr = gam(
  compo ~ high_15 + fontaine_cat + bmi_cat + asa_cat + age_cat +
    nyha_cat + renalebegleiterkrankungen + propensity_score +
    share_open + endo_real_max + anzahl_betten + uhc +
    gefaess_leistung_vc71 + gefaess_leistung_vc67 + zertifizierung +
    s(facilityid, bs = "re"),
  data = dt.open.b, family = binomial(link = "logit")
)

sum_ftr = summary(model_ftr)
est_ftr = sum_ftr$p.coeff[-1]
se_ftr = sum_ftr$se[c(-1,-28)]

forest_tab = data.table(
  variable = names(est_ftr),
  Estimates = round(exp(est_ftr), 3),
  LL = round(exp(est_ftr - 1.96 * se_ftr), 3),
  UL = round(exp(est_ftr + 1.96 * se_ftr), 3)
)
forest_tab

cat_vars = c(
  "high_15", "fontaine_cat", "bmi_cat", "asa_cat", "age_cat",
  "renalebegleiterkrankungen", "nyha_cat", "uhc",
  "gefaess_leistung_vc71", "gefaess_leistung_vc67","zertifizierung"
)

res_cat = rbindlist(lapply(cat_vars, function(v) {
  x = dt.open.b[!is.na(get(v)), .(total = .N, level = get(v)), by = v]
  if (is.factor(dt.open.b[[v]]) || is.character(dt.open.b[[v]])) {
    x[, level := factor(level, levels = levels(factor(dt.open.b[[v]])))]
    x = x[order(level)]
  } else {
    x = x[order(level)]
  }
  x[, .(variable = paste0(v, level), variable2 = v, variable3 = as.character(level), total)]
}), fill = TRUE)

res_num = data.table(
  variable = c("propensity_score",  "endo_real_max", "anzahl_betten"),
  variable2 = c("propensity_score", "endo_real_max", "anzahl_betten"),
  variable3 = c("propensity_score", "endo_real_max", "anzahl_betten"),
  total = nrow(dt.open.b)
)
res_num

res_display = rbind(res_cat, res_num, fill = TRUE)
res_display[, row_nr := .I]

setkey(res_display, variable)
setkey(forest_tab, variable)
forest_tab = forest_tab[res_display, ][order(row_nr)]

forest_tab[is.na(Estimates), variable := variable2]
forest_tab[!is.na(Estimates), variable := paste0("   ", variable3)]
forest_tab[is.na(Estimates), variable := variable2]
forest_tab[!is.na(Estimates), variable := paste0("   ", variable3)]

forest_tab[variable == "high_15", variable := "Highest annual volume, ref.: 0–14 cases"]
forest_tab[variable == "fontaine_cat", variable := "Fontaine stage, ref.: none"]
forest_tab[variable == "bmi_cat", variable := "BMI category, kg/m², ref.: 18.5–24.9"]
forest_tab[variable == "asa_cat", variable := "ASA classification, ref.: ASA III"]
forest_tab[variable == "age_cat", variable := "Age group, years, ref.: 65–74"]
forest_tab[variable == "nyha_cat", variable := "NYHA class, ref.: none"]
forest_tab[variable == "renalebegleiterkrankungen", variable := "Renal comorbidities, ref.: none"]
forest_tab[variable == "uhc", variable := "University hospital center"]
forest_tab[variable == "gefaess_leistung_vc71", variable := "Emergency surgical care available"]
forest_tab[variable == "gefaess_leistung_vc67", variable := "Vascular intensive care unit available"]
forest_tab[variable == "peripheral", variable := "Center located in peripheral district"]
forest_tab[variable == "zertifizierung", variable := "Certified center (DGG)"]
#forest_tab[variable == "   share_open", variable := "Maximum share of open repairs, change 1 SD"]
forest_tab[variable == "   endo_real_max", variable := "Maximum EVAR volume, change 1 SD"]
forest_tab[variable == "   anzahl_betten", variable := "Total number of beds, change 1 SD"]
forest_tab[variable == "   propensity_score", variable := "Propensity open repair, change 1 SD"]

forest_tab[, `:=`(Variable = variable, N = total)]

forest_plot_os = forester(
  left_side_data = forest_tab[, .(Variable, N)],
  estimate = forest_tab$Estimates,
  ci_low = forest_tab$LL,
  ci_high = forest_tab$UL,
  display = TRUE,
  estimate_precision = 2,
  file_path = "figure_2_b.png",
  dpi = 600,
  xlim = c(0.1, 10),
  null_line_at = 1,
  font_family = "Fira Sans",
  x_scale_linear = FALSE,
  arrows = TRUE,
  arrow_labels = c("lower", "higher"),
  estimate_col_name = "Odds Ratio (95% CI)",
  point_sizes = 3
)
# =============================================================================
# 9. DESCRIPTIVE TABLES
# =============================================================================
# Table 1 – Patient characteristics by high-volume (≥15)
patient_vars = c("sex", "fontaine_cat", "age_cat", "bmi_cat", "asa_cat", "diam_cat",
                 "diabetes", "hypertonus", "dyslipoprot", "nikotinakt", "aneurysmaand",
                 "renalebegleiterkrankungen", "pulmonalebegleiterkrankungen",
                 "kardialebegleiterkrankungen", "iliacalebeteiligung", "lokalisation",
                 "nyha_cat", "copd_cat", "koronareherzerkrankung", "prior_aorta_surgery")
cat_vars = patient_vars[-1]

tab1 = CreateTableOne(vars = patient_vars, strata = "high_15", data = dt.open,
                      factorVars = cat_vars, addOverall = TRUE, includeNA = TRUE, test = FALSE)
tab1_df = print(tab1, noSpaces = TRUE, smd = TRUE, includeNA = TRUE,
                catDigits = 1, contDigits = 1, pDigits = 2, showAllLevels = FALSE, explain = FALSE)
tab1_df = as.data.frame(cbind(tab1_df, rownames(tab1_df)))
names(tab1_df)[5] = "Variable"
rownames(tab1_df) = NULL

var_labels = c(
  "fontaine_cat" = "Fontaine stage", "age_cat" = "Age group, years",
  "bmi_cat" = "BMI category, kg/m²", "asa_cat" = "ASA classification",
  "diam_cat" = "Aneurysm diameter, mm", "diabetes" = "Diabetes mellitus",
  "hypertonus" = "Arterial hypertension", "dyslipoprot" = "Dyslipidemia",
  "nikotinakt" = "Current smoking", "aneurysmaand" = "Other aneurysm",
  "renalebegleiterkrankungen" = "Renal comorbidities",
  "pulmonalebegleiterkrankungen" = "Pulmonary comorbidities",
  "kardialebegleiterkrankungen" = "Cardiac comorbidities",
  "iliacalebeteiligung" = "Iliac involvement",
  "lokalisation" = "Juxtarenal location", "nyha_cat" = "NYHA class",
  "copd_cat" = "COPD GOLD stage", "koronareherzerkrankung" = "Coronary artery disease",
  "prior_aorta_surgery" = "Prior aortic surgery or stent"
)
tab1_df$Variable = sapply(tab1_df$Variable, function(x) {
  varname = trimws(sub("\\s*=.*$", "", x))
  if (varname %in% names(var_labels)) sub(varname, var_labels[varname], x, fixed = TRUE) else x
})
tab1_df$Variable = gsub(" = yes", "", tab1_df$Variable, fixed = TRUE)
tab1_df = cbind(tab1_df[, 5, drop = FALSE], tab1_df[, 1:4])

ft_tab1 = flextable(tab1_df)
ft_tab1 = font(ft_tab1, fontname = "Arial", part = "all")
ft_tab1 = fontsize(ft_tab1, size = 8, part = "all")
ft_tab1 = autofit(ft_tab1)
save_as_docx(ft_tab1, path = paste0(Sys.Date(), "_table1_patients_high_15.docx"))

# Table 2 – Centre characteristics by high-volume
dt.open[,high_10_num:=ifelse(high_10=="≥10 cases",1L,0L)]
dt.open[,high_15_num:=ifelse(high_15=="≥15 cases",1L,0L)]
dt.open[,sex_num:=ifelse(sex=="women",1,0)]

pct_vars = c("tod", "ftr", "major_complication", "compo", "sex_num")
dt.open[, (pct_vars) := lapply(.SD, `*`, 100), .SDcols = pct_vars]

outcomes     = c("tod", "ftr", "major_complication", "compo")
patient_vars = c("alterbeiop","bmi","maximalerdurchmesserdesaneurysmas","sex_num")
center_vars = c("high_15_num", "open_real_max", "open_n", "share_open",
                "open_real_max","endo_n","endo_real_max", "uhc", "anzahl_betten", "gefaess_leistung_vc71",
                "gefaess_leistung_vc67","zertifizierung")

dt.center = dt.open[, c(
  .SD[, lapply(.SD, mean, na.rm = TRUE), .SDcols = c(outcomes, patient_vars)],
  .SD[, lapply(.SD, function(x) max(x, na.rm = TRUE)), .SDcols = center_vars]
), by = facilityid]

patient_vars = setdiff(c(outcomes,patient_vars,center_vars),"high_15_num")
cat_vars2 = c("uhc", "gefaess_leistung_vc71", "gefaess_leistung_vc67",
               "zertifizierung")
tab2 = CreateTableOne(vars = patient_vars,strata = "high_15_num", data = dt.center,
                      factorVars = cat_vars2, addOverall = TRUE, includeNA = TRUE, test = FALSE)

tab2_df = print(tab2, noSpaces = TRUE, smd = TRUE, includeNA = TRUE,
                catDigits = 1, contDigits = 1, pDigits = 2, showAllLevels = FALSE)
tab2_df = as.data.frame(tab2_df)
tab2_df = cbind(Variable = rownames(tab2_df), tab2_df)
rownames(tab2_df) = NULL

var_labels2 = c(
  NA, "Death rate", "Failure-to-rescue rate", "Major complication rate",
  "Major complication or death rate", "Mean patient age, years", "Mean BMI, kg/m²",
  "Mean aneurysm diameter, mm", "Share of women", "Maximum annual volume",
  "Maximum annual volume (quality report)", "Maximum share of open repairs",
  "Maximum EVAR volume", "Maximum EVAR volume  (quality report)","University hospital center", "Total number of beds",
  "Emergency surgical care available", "Vascular intensive care unit available", "Certified center (DGG)"
)
tab2_df$Variable = var_labels2
colnames(tab2_df)[3:4] = c("1-14 cases", "≥15 cases")

ft_tab2 = flextable(tab2_df)
ft_tab2 = font(ft_tab2, fontname = "Arial", part = "all")
ft_tab2 = fontsize(ft_tab2, size = 8, part = "all")
ft_tab2 = autofit(ft_tab2);ft_tab2
save_as_docx(ft_tab2, path = paste0(Sys.Date(), "_table2_centers.docx"))

# Patients
uniqueN(dt.open)
# Centres
uniqueN(dt.open$facilityid)

table(dt.open$jahr)

# =============================================================================
# END OF SCRIPT
# =============================================================================