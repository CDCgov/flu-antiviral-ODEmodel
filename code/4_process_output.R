## Set up ---------

xfun::pkg_attach(c("patchwork", "ggrepel", "scico", "RColorBrewer", "here", "tidyverse", "janitor", "lhs", "gtools"))

filepath <- ""

nsim <- 100

# load functions and data for calibration plots
source(here(filepath, 'code/0_functions.R'))
source(here(filepath, 'code/0_params.R'))
source(here(filepath, 'code/0_calibration_data.R'))

# plotting parameters
age_labels <- c("0-4 yrs", "5-17 yrs", "18-49 yrs", "50-64 yrs", "65+ yrs")
age_labels_nice <- c("0-4 yrs", "5-17 yrs", "18-49 yrs", "50-64 yrs", "≥65 yrs")

cols_severity <- brewer.pal(6, 'Spectral')[c(1, 2, 6)]
cols_age <- scico(7, palette = 'lapaz')[c(2:6)]

figpath <- "figs"

# scenario information
scenario_tab <- read.csv(here(filepath,'inputs/scenarios.csv'), row.names = 1)
epi_scenarios <- names(scenario_tab)

scenario_short <- unique(str_split_i(epi_scenarios, pattern = "_", i = 3))
scenario_labels <- c("Baseline", "No transmission impact", "Small transmission impact", "Large transmission impact")


# Process output files ------------------

# main scenario analysis (from 1_run_model_original.R)
load(here(filepath, paste0("out/original_100.RData")) )

out_all <- combine_out(out) %>% 
  select(time, scenario, sim, contains("all")) %>% 
  separate(scenario, c("X", "severity", "scenario"), sep = "_") %>%
  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels))

# sensitivity scenario analysis (from 2_run_model_sensitivity.R)
load(here(filepath, paste0("out/sensitivity_100.RData")) )

out_all_sens <- combine_out(out_sens) %>% 
  select(time, scenario, sim, contains("all")) %>% 
  separate(scenario, c("X", "severity", "scenario"), sep = "_") %>%
  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels))

# extended spectrum analyses (from 3_run_model_spectrum.R)
load(here(filepath, paste0("out/spectrum_0.025.RData")))


# Plot contact matrix input --------

contact_mat %>% mutate(.row = colnames(.)) %>%
  pivot_longer(-.row, names_to = ".col", values_to = "val") %>%
  mutate(
    .row = factor(paste(.row, "yrs"), levels = paste(colnames(contact_mat), "yrs")),     
    .col = factor(paste(.col, "yrs"), levels = paste(colnames(contact_mat), "yrs"))     
  ) %>%
  ggplot(aes(x = .row, y = .col, fill = val)) + 
  geom_tile(color = "white") +
  scale_fill_scico(palette = "lapaz", name = "Average\ncontacts\nper day") +
  coord_fixed() + get_theme(txt = 14) +
  scale_x_discrete(labels = age_labels_nice) + scale_y_discrete(labels = age_labels_nice) +
  labs(x = "Respondent age group", y = "Contact age group")

# save_pdf function to save plot with special character '≥' (see 0_functions.R)
save_pdf(here(figpath, "input_contacts.pdf"), width = 8, height = 6) 


# Calibration plot: incidence of symptomatic illnesses --------------

origin <- as.Date("2010-08-01")  # choose any non-leap year origin

# set earliest / latest observed peaks from past seasons
# dates are approximate as peaks are reported in weeks
# doesn't include early peak of 2022/23 
rect_df <- data.frame(
  xmin = as.Date("2010-12-20"),  # late-December
  xmax = as.Date("2011-03-17"),  # mid-March
  ymin = -Inf,
  ymax = Inf
)

p_timing <- 
  out_all %>% 
  group_by(sim, severity, scenario, scenario_names) %>%
  mutate(I_inc = c(0, diff(count_Is_all + count_Ia1_all + count_Ia2_all))) %>%
  ungroup() %>%
  # approximate time entering symptomatic infectious rather than leaving
  mutate(time = time - 3) %>% filter( time >=0) %>%
  mutate(
    date  = origin + (time - 1),
    month = month(date, label = TRUE, abbr = TRUE),  
    month = factor(month, levels = month.abb[c(7:12, 1:6)], ordered = TRUE)
  ) %>%
  ggplot(aes(x = date, y = I_inc, 
             color = severity, group = interaction(severity, sim))) + 
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey60", alpha = 0.25) +
  geom_line(linewidth = 0.5) + get_theme(txt = 11) + facet_wrap(~ scenario_names, ncol = 2) +
  scale_y_continuous(labels = scales::comma) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual("Season\nseverity", values = cols_severity[c(1, 3)]) +
  labs(x = "Month", y = "New symptomatic illnesses")

ggsave(p_timing, file = here(figpath, "calibration_curve.pdf"), device = cairo_pdf, width = 10, height = 6)


# Calibration plot: total symptomatic illnesses ----------------

Is_summary <- out_all %>% filter(time == max(time)) %>%
  mutate(total_Is = count_Is_all + count_Ia1_all + count_Ia2_all) %>%
  group_by(scenario, scenario_names, severity) %>%
  # get mean and 2.5th-97.5th percentiles (see 0.functions.R)
  get_summary(contains("total")) %>%
  gather(var, val, contains("total")) %>%
  separate(var, c("X", "var", "metric"), "_") %>%
  spread(metric, val)

p_Is_base <- 
  Is_summary %>% filter(var == "Is", scenario %in% c("baseline")) %>%
  ggplot() +
  geom_point(data = db, aes(x = season, y = symptomatic_illnesses, color = severity), size = 3) +
  geom_errorbar(data = db, inherit.aes = F,
                aes(x = season, color = severity,
                    ymin = symptomatic_illnesses_95_ui_lower_bound, 
                    ymax = symptomatic_illnesses_95_ui_upper_bound), 
                width = 0.1, position = position_dodge(width = 0.2)) + 
  geom_hline(aes(yintercept = mean, color = severity), alpha = 0.35, linewidth = 1.5) +
  get_theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_manual("Season\nseverity", breaks = c("high", "moderate", "low"), values = cols_severity) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Season", y = "Symptomatic illnesses", title = "Baseline scenario")

p_Is_bigtrans <- 
  Is_summary %>% filter(var == "Is", scenario %in% c("bigtrans")) %>%
  ggplot() +
  geom_point(data = db, aes(x = season, y = symptomatic_illnesses, color = severity), size = 3) +
  geom_errorbar(data = db, inherit.aes = F,
                aes(x = season, color = severity,
                    ymin = symptomatic_illnesses_95_ui_lower_bound, 
                    ymax = symptomatic_illnesses_95_ui_upper_bound), 
                width = 0.1, position = position_dodge(width = 0.2)) + 
  geom_hline(aes(yintercept = mean, color = severity), alpha = 0.35, linewidth = 1.5) +
  get_theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_manual("Season\nseverity", breaks = c("high", "moderate", "low"), values = cols_severity) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Season", y = "Symptomatic illnesses", title = "Large transmission impact scenario")

p_Is_base / p_Is_bigtrans + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")
ggsave(file = here(figpath, "calibration_ill.pdf"), width = 10, height = 11)


# Calibration plot: total hospitalizations -----------------

hosp_summary <- 
  out_all %>% filter(time == max(time)) %>% 
  group_by(scenario, severity) %>%
  get_summary(H_all)

p_hosp <- 
  hosp_summary %>% filter(scenario %in% c("baseline")) %>%
  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels)) %>%
  ggplot() +
  geom_point(data = db, aes(x = season, y = hospitalizations, color = severity), size = 3) +
  geom_errorbar(data = db, inherit.aes = F,
                aes(x = season, color = severity,
                    ymin = hospitalizations_95_ui_lower_bound, 
                    ymax = hospitalizations_95_ui_upper_bound), 
                width = 0.1, position = position_dodge(width = 0.2)) + 
  geom_hline(aes(yintercept = H_all_mean, color = severity), alpha = 0.35, linewidth = 1.5) +
  get_theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_manual("Season\nseverity", breaks = c("high", "moderate", "low"), values = cols_severity) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Season", y = "Hospitalizations")
ggsave(p_hosp, file = here(figpath, "calibration_hosp.png"), width = 10, height = 6, dpi = 600)


hosp_summary %>% filter(scenario %in% c("bigtrans")) %>%
  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels)) %>%
  ggplot() +
  geom_point(data = db, aes(x = season, y = hospitalizations, color = severity), size = 3) +
  geom_errorbar(data = db, inherit.aes = F,
                aes(x = season, color = severity,
                    ymin = hospitalizations_95_ui_lower_bound, 
                    ymax = hospitalizations_95_ui_upper_bound), 
                width = 0.1, position = position_dodge(width = 0.2)) + 
  geom_hline(aes(yintercept = H_all_mean, color = severity), alpha = 0.35, linewidth = 1.5) +
  get_theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_manual("Season\nseverity", breaks = c("high", "moderate", "low"), values = cols_severity) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Season", y = "Hospitalizations")
ggsave(file = here(figpath, "calibration_hosp_large.pdf"), width = 10, height = 6)


# Calibration plot: age distributions -----------

tmp <- out %>% 
  select(sim, scenario, time, contains("count_I"), matches("H[0-9]")) %>% 
  filter(time == max(time)) %>% 
  separate(scenario, c("X", "severity", "scenario"), sep = "_")

tmp_ill <- 
  tmp %>% gather(var, val, contains("count")) %>%
  mutate(agegrp = substr(var, nchar(var), nchar(var))) %>% 
  group_by(sim, severity, scenario, agegrp) %>%
  summarize(illnesses = sum(val)) 

tmp_hosp <- 
  tmp %>% gather(var, val, matches("H[0-9]")) %>%
  mutate(agegrp = substr(var, nchar(var), nchar(var))) %>% 
  group_by(sim, severity, scenario, agegrp) %>%
  summarize(hosp = sum(val)) 

# calibration data
dat <- db_age %>% 
  select(season, severity, age_txt, illnesses = symptomatic_illnesses, hosp = hospitalizations) %>%
  group_by(severity, season) %>%
  mutate(across(c(illnesses, hosp), function(x) x/sum(x))) %>%
  filter(severity %in% c("high", "low"))

p_age <-
  left_join(tmp_ill, tmp_hosp) %>%
  mutate(across(c(illnesses, hosp), function(x) x/sum(x))) %>%
  # average across simulations
  group_by(severity, scenario, agegrp) %>%
  summarize(across(c(illnesses, hosp), mean)) %>% 
  ungroup() %>%
  filter(scenario == "baseline") %>% 
  select(-scenario) %>%
  mutate(age_txt = factor(age_labels[as.numeric(agegrp)], levels = age_labels),
         season = "Model") %>%
  bind_rows(dat) %>%
  gather(var, prop, illnesses, hosp) %>%
  mutate(var = ifelse(var == "illnesses", "Symptomatic illnesses", "Hospitalizations"),
         var = factor(var, levels = c("Symptomatic illnesses", "Hospitalizations")),
         severity = ifelse(severity == "high", "High severity", "Low severity")) %>%
  ggplot(aes(x = season, y = prop, fill = age_txt, colour = season == "Model")) + 
  geom_bar(stat = "identity", position = "stack", linewidth = 0.7) + 
  facet_grid(var ~ severity, scales = "free") +
  scale_y_continuous("Outcomes (% of total)", labels = scales::percent) +
  labs(x = "Season") + 
  scale_fill_manual("Age group", values = cols_age, labels = age_labels_nice) +
  get_theme(txt = 11, axis.text.x = element_text(angle = 90)) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NULL), na.value  = NA, guide = "none")

# save_pdf function to save plot with special character '≥' (see 0_functions.R)
save_pdf(p_age, file = here(figpath, "calibration_age.pdf"), width = 8, height = 6)


# Plot hospitalizations prevented (all ages) ------------------------

tmp_H <- out_all %>% 
  filter(time == max(time)) %>%
  select(sim, severity, scenario, H_all) 

prevented_H <- 
  tmp_H %>% spread(scenario, H_all) %>%
  # calc absolute difference and percentage difference  
  mutate(across(.cols = contains("trans"), 
                .fns = list(directperc   = function(x) (baseline - notrans)/baseline * 100,
                            indirectperc = function(x) (notrans - x)/baseline * 100,
                            totalperc    = function(x) (baseline - x)/baseline * 100,
                            
                            directabs   = function(x) baseline - notrans,
                            indirectabs = function(x) notrans - x,
                            totalabs    = function(x) baseline - x),
                .names = "{.col}_{.fn}")) %>%
  # rename the original columns to follow the same pattern as _abs/_perc
  rename_with(~ paste0(.x, "_orig"), .cols = all_of(scenario_short[scenario_short != "baseline"])) %>% 
  # then gather
  pivot_longer(cols = contains("trans"), names_to   = c("scenario", ".value"), names_sep  = "_") %>%
  # then summarize
  group_by(severity, scenario) %>%
  get_summary(contains(c("orig", "direct", "indirect", "total"))) %>%
  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels))


p_h_total <- 
  prevented_H %>% filter(scenario != "baseline") %>%
  ggplot(aes(x = scenario_names, y = totalabs_mean, color = severity, group = severity)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) + 
  geom_errorbar(aes(ymin = totalabs_lcl, ymax = totalabs_ucl), position = position_dodge(width = 0.5), width = 0.3) +
  geom_text_repel(aes(label = paste0(sprintf("%.1f", totalperc_mean), "%")), size = 3.5,
                  position = position_dodge(width = 0.5), hjust= -0.3, color = "black") +
  get_theme(txt = 11) + 
  scale_x_discrete(labels = function(l) str_wrap(l, width = 10)) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual("Season severity", values = cols_severity[c(1, 3)]) +
  labs(title = "Total effects", x = "Scenario", y = "Number of hospitalizations prevented")

p_h_ind <- 
  prevented_H %>% filter(scenario != "baseline") %>%
  ggplot(aes(x = scenario_names, y = indirectabs_mean, color = severity, group = severity)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) + 
  geom_errorbar(aes(ymin = indirectabs_lcl, ymax = indirectabs_ucl), position = position_dodge(width = 0.5), width = 0.3) +
  geom_text_repel(aes(label = paste0(sprintf("%.1f", indirectperc_mean), "%")), size = 3.5,
            position = position_dodge(width = 0.5), hjust= -0.3, color = "black") +
  get_theme(txt = 11) + 
  scale_x_discrete(labels = function(l) str_wrap(l, width = 10)) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual("Season severity", values = cols_severity[c(1, 3)]) +
  labs(title = "Indirect effects", x = "Scenario", y = "Number of hospitalizations prevented")

p_h <- p_h_total + p_h_ind + plot_layout(guides = "collect") + 
        plot_annotation(tag_levels = 'A') & theme(legend.position = "bottom")
ggsave(p_h, file = here(figpath, "prevented_h_all.png"), width = 8.5, height = 4.5, dpi = 600)


# Percentage prevented
prevented_H %>% select(severity, scenario_names, contains("perc_mean")) %>% 
  arrange(scenario_names)

# Absolute number prevented
prevented_H %>% select(severity, scenario_names, contains("totalabs")) %>% 
  arrange(scenario_names)


# Plot extended analyses: spectrum of av effect & % undergoing prompt tx  ------------------

# For visualization purposes, plot p_av marker at mean of previously sampled ranges
if(!exists("pars_fixed")) {
  pars_fixed <- 
    data.frame(vals = (maxvals + minvals)/2,
               pars = c(paste0(parbyage, 1:m), parnames[parnames != parbyage])) %>% spread(pars, vals) %>% 
    pivot_longer(cols = contains(parbyage), values_to = parbyage) %>%
    mutate(agegrp = as.numeric(str_match(name, pattern = "[0-9]"))) %>% select(-name) %>% 
    mutate(p_prescribed1 = p_seekcare * p_d1 * p_prescribed,
           p_prescribed2 = p_seekcare * p_d2 * p_prescribed,
           
           p_a1 = p_prescribed1 * p_adherence,
           p_a2 = p_prescribed2 * p_adherence
    )
}

plt_pts <- 
  data.frame(scenario = c("No transmission impact", "Small transmission impact", "Large transmission impact"),
           av_hosp = c(0.8, 0.8, 0.8),
           av_trans = c(1, 0.9, 0.8),
           p_av = mean(pars_fixed$p_a1 + pars_fixed$p_a2),
           prevented = 0,
           tx = 0)


# AV effectiveness against hospitalization vs transmission
# comparison contour lines
comp_lines <- c(1, 5, 10, 15)

tmp_hosp_trans <-  
  out_spectrum %>% filter(combination == "av_trans x av_hosp") %>%
  left_join(
    out_spectrum %>% filter(av_trans == 1, av_hosp == 1, combination == "av_trans x av_hosp") %>% 
      rename(baseline_H = H_all) %>% select(severity, scenario, baseline_H)
  ) %>%
  mutate(severity = paste(str_to_title(severity), "severity season")) %>%
  mutate(prevented = (baseline_H - H_all)/baseline_H * 100) 

# quoted numbers (effectiveness = 100%)
tmp_hosp_trans %>% filter(av_trans == 0, av_hosp == 0)


p_trans_hosp <- tmp_hosp_trans %>%
  ggplot(aes(x = (1 - av_hosp)*100, y = (1 - av_trans)*100, fill = prevented, color = prevented)) +
  geom_tile() + facet_wrap(~ severity) +
  geom_contour(aes(z = prevented), breaks = comp_lines, 
               linetype = "dashed", colour = "white", linewidth = 0.4, alpha = 0.5) +
  stat_contour(aes(z = prevented, label = after_stat(level)), breaks = comp_lines,
               geom = "text", color = "white", size = 2.5, check_overlap = TRUE, alpha = 0.5) +
  geom_point(data = plt_pts, col = "white") +
  geom_text(data = plt_pts, aes(label = scenario), col = "white", hjust = -0.05, vjust = -0.5) +
  scale_fill_scico("Prevented\nhospitalizations\n(%)", palette = "lapaz") +
  scale_color_scico(guide = "none", palette = "lapaz") +
  get_theme(txt = 11.5) +
  labs(x = "Reduction in hospitalization risk (%)", y = "Reduction in transmission risk (%)")


# % undergoing prompt Tx vs AV effectiveness against transmission
# comparison contour lines
comp_lines2 <- c(10, 30, 50, 70, 90)

tmp_pav_trans <- 
  out_spectrum %>% filter(combination == "av_trans x p_av") %>%
  left_join(
    out_spectrum %>% filter(av_trans == 1, p_av == 0, combination == "av_trans x p_av") %>% 
      rename(baseline_H = H_all) %>% select(severity, scenario, baseline_H)
  ) %>%
  mutate(severity = paste(str_to_title(severity), "severity season")) %>%
  mutate(prevented = (baseline_H - H_all)/baseline_H * 100)


# quoted numbers - both 100%
tmp_pav_trans %>% filter(av_trans == 0, p_av == 1)

# quoted numbers - 20% effectiveness and 25% on treatment
tmp_pav_trans %>% filter(av_trans == 0.8, p_av == 0.25)

# quoted numbers - 20% effectiveness and 100% on treatment
tmp_pav_trans %>% filter(av_trans == 0.8, p_av == 1)


p_trans_pav <- 
  tmp_pav_trans %>%
  ggplot(aes(x = p_av*100, y = (1 - av_trans)*100, fill = prevented, color = prevented)) +
  geom_tile() + facet_wrap(~ severity) +
  geom_contour(aes(z = prevented), breaks = comp_lines2, 
               linetype = "dashed", colour = "white", linewidth = 0.4, alpha = 0.5) +
  stat_contour(aes(z = prevented, label = after_stat(level)), breaks = comp_lines2,
               geom = "text", color = "white", size = 2.5, check_overlap = TRUE, alpha = 0.5) +
  geom_point(data = plt_pts, col = "white") +
  geom_text(data = plt_pts, aes(label = str_wrap(scenario, width = 30)), 
            col = "white", hjust = -0.05, vjust = -0.2) +
  scale_fill_scico("Prevented\nhospitalizations\n(%)", palette = "lapaz") +
  scale_color_scico(guide = "none", palette = "lapaz") +
  get_theme(txt = 11.5) +
  labs(x = "Symptomatic individuals undergoing prompt treatment (%)", y = "Reduction in transmission risk (%)")

p_trans_hosp / p_trans_pav + plot_annotation(tag_levels = 'A')
ggsave(file = here(figpath, "spectrum_both.png"), width = 8.5, height = 8.5, dpi = 800)


# Increasing each component of p_av --------------------------

# get parameters comprising the % people completing prompt Tx
pars <- c("pars_fixed$p_seekcare * (pars_fixed$p_d1 + pars_fixed$p_d2)", "pars_fixed$p_prescribed", "pars_fixed$p_adherence")

# find all permutations 
perms_id <- permutations(3, 3)
perms <- matrix(pars[perms_id], nrow = nrow(perms_id))
colnames(perms) <- c("par1", "par2", "par3")

perms_out <- perms %>% data.frame() %>% 
  mutate(iteration = 1:nrow(perms)) %>% 
  rowwise() %>%
  mutate(out1 = mean( Reduce("*", lapply(pars[pars != par1], function(expr) eval(parse(text = expr)))) ),
         out2 = mean( Reduce("*", lapply(pars[!(pars %in% c(par1, par2))], function(expr) eval(parse(text = expr)))) ),
         out3 = 1 )  %>% ungroup() %>% 
  mutate(across(par1:par3, function(x) gsub(x = x, pattern = "pars_fixed\\$", replacement = "")),
         across(par1:par3, function(x) ifelse(x == "p_seekcare * (p_d1 + p_d2)", "p_seekcare_2d", x)))


# find prevented hospitalizations at each new value
tmp <- data.frame(iteration = rep(1:nrow(perms_out), each = 2*nrow(plt_pts)),
                  severity  = rep(rep(c("Low severity season", "High severity season"), each = nrow(plt_pts)), times = nrow(perms_out)),
                  av_trans = rep(rep(plt_pts$av_trans, times = 2), times = nrow(perms_out)) )

dat_plt <- 
  tmp_pav_trans %>% filter(av_trans %in% plt_pts$av_trans) %>%
  select(av_trans, p_av, severity, prevented, contains("count")) %>% 
  full_join(
    perms_out %>% left_join(tmp) 
  ) %>%
  group_by(av_trans, severity, iteration) %>%
  mutate(prevented1 = prevented[which.min(abs(p_av - out1))],
         prevented2 = prevented[which.min(abs(p_av - out2))],
         prevented3 = prevented[which.min(abs(p_av - out3))], 
         
         tx = count_Ia1_all + count_Ia2_all, 
         tx1 = tx[which.min(abs(p_av - out1))],
         tx2 = tx[which.min(abs(p_av - out2))],
         tx3 = tx[which.min(abs(p_av - out3))]
  ) %>% ungroup() %>%
  distinct(av_trans, severity, iteration, 
           par1, par2, par3, 
           prevented1, prevented2, prevented3, 
           tx1, tx2, tx3) 

# Max treatment #s
dat_plt %>% distinct(tx3, severity, av_trans)

cols_params <-  scico(7, palette = 'lapaz')[c(3, 4, 6)]


# Function to plot different numbers of sequential increases

create_barplot <- function(dat){
  tmp <- dat %>% 
    group_by(av_trans, severity, iteration) %>%
    mutate(diff = c(prevented[1], diff(prevented)),
           parlabel = case_when(par == "p_seekcare_2d" ~ "Sought prompt care",
                                par == "p_prescribed" ~ "Prescribed antivirals",
                                par == "p_adherence" ~ "Adhered whilst infectious"),
           parlabel = factor(parlabel, 
                             levels = c("Sought prompt care", 
                                        "Prescribed antivirals", 
                                        "Adhered whilst infectious"))) %>%
    left_join(plt_pts %>% select(av_trans, scenario)) %>% ungroup() %>% 
    mutate(block = ceiling(iteration/factorial(length(pars)-1))) 
  
  keep <- tmp %>%
    group_by(block, iteration, av_trans, scenario, severity) %>%
    summarise(input_string = paste(par, collapse = "_")) %>%
    group_by(block, av_trans, scenario, severity) %>%
    filter(!duplicated(input_string)) %>%
    select(-input_string)
  
  tmp %>% 
    semi_join(keep, by = c("block", "iteration", "av_trans", "scenario", "severity")) %>%
    mutate(scenario = factor(scenario, levels = rev(scenario_labels))) %>%
    ggplot(aes(x = iteration, y = diff, fill = parlabel, group = order)) +
    geom_bar(stat = "identity", position = position_stack(reverse = T), color = "black", linewidth = 0.5) + 
    facet_grid(scenario ~ severity) +
    scale_fill_manual("Component", values = cols_params) + 
    labs(x = "Permutation", y = "Prevented hospitalizations (%)") +
    get_theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}


# Actual plotting

# Add lines to show where you started 
start_pts <- prevented_H %>% filter(scenario != "baseline") %>%
  select(scenario = scenario_names, severity, totalperc_mean) %>%
  mutate(severity = paste(str_to_title(severity), "severity season"))


# Stacked bar plot (first increase)
dat_plt %>% pivot_longer(cols = c(par1, prevented1),
                         names_to = c(".value", "order"),
                         names_pattern = "([a-z]+)(\\d+)") %>%
  create_barplot() + 
  geom_hline(data = start_pts, aes(yintercept = totalperc_mean), linetype = "dashed") +
  xlab("Component")
ggsave(file = here(figpath, "spectrum_extension_bar1.png"), width = 8.5, height = 7, dpi = 600)


# Stacked bar plot (first & second increases)
dat_plt %>% pivot_longer(cols = c(par1:par2, prevented1:prevented2),
                         names_to = c(".value", "order"),
                         names_pattern = "([a-z]+)(\\d+)") %>%
  create_barplot() + 
  geom_hline(data = start_pts, aes(yintercept = totalperc_mean), linetype = "dashed")
ggsave(file = here(figpath, "spectrum_extension_bar2.pdf"), width = 8.5, height = 7)


# Stacked bar plot (all sequential increases)
dat_plt %>% pivot_longer(cols = c(par1:par3, prevented1:prevented3),
                         names_to = c(".value", "order"),
                         names_pattern = "([a-z]+)(\\d+)") %>%
  create_barplot()  + 
  geom_hline(data = start_pts, aes(yintercept = totalperc_mean), linetype = "dashed")
ggsave(file = here(figpath, "spectrum_extension_bar3.pdf"), width = 8.5, height = 7)


# quoted max impacts
dat_plt %>% pivot_longer(cols = c(par1:par3, prevented1:prevented3),
                         names_to = c(".value", "order"),
                         names_pattern = "([a-z]+)(\\d+)") %>% 
  filter(order == 3) %>% distinct(prevented, av_trans, severity)


# Plot number undergoing prompt treatment (spectrum) ------------------

# comparison contour lines
comp_lines_tx <- c(1e6, 1e7, 2e7, 3e7)

p_tx_pav <- 
  tmp_pav_trans %>% mutate(tx = count_Ia1_all + count_Ia2_all) %>%
  ggplot(aes(x = p_av*100, y = (1 - av_trans)*100, fill = tx, color = tx)) +
  geom_tile() + facet_wrap(~ severity) +
 geom_contour(aes(z = tx), breaks = comp_lines_tx, 
               linetype = "dashed", colour = "white", linewidth = 0.4, alpha = 0.5) +
 stat_contour(aes(z = tx, label = after_stat(paste0(level / 1e6, "M"))), breaks = comp_lines_tx,
                geom = "text", color = "white", size = 2.5, check_overlap = TRUE, alpha = 0.5) +
 geom_point(data = plt_pts, col = "white") +
 geom_text(data = plt_pts, aes(label = str_wrap(scenario, width = 30)), 
            col = "white", hjust = -0.05, vjust = -0.2) +
  scale_fill_scico("Number\ntreated", palette = "lapaz", labels = function(x) ifelse(x == 0, "0", paste0(x / 1e6, "M"))) +
  scale_color_scico(guide = "none", palette = "lapaz") +
  get_theme(txt = 11.5) +
  labs(x = "Symptomatic individuals undergoing prompt treatment (%)", y = "Reduction in transmission risk (%)")
ggsave(p_tx_pav, file = here(figpath, "spectrum_tx.pdf"), width = 10, height = 5)

# quoted numbers
tmp_pav_trans %>% mutate(tx = count_Ia1_all + count_Ia2_all) %>%
  filter(p_av == 0.25, av_trans == 0.8)
  
tmp_pav_trans %>% mutate(tx = count_Ia1_all + count_Ia2_all) %>%
  filter(p_av == 1, av_trans %in% c(0, 1))

  
# Plot sensitivity analyses (scenario & extended): hosp risk reduction = 10% --------------

tmp_H_sens <- out_all_sens %>% 
  filter(time == max(time)) %>%
  select(sim, severity, scenario, H_all) 

sens_H <- 
  tmp_H_sens %>% spread(scenario, H_all) %>%
  # calc absolute difference and percentage difference  
  mutate(across(.cols = contains("trans"), 
                .fns = list(directperc   = function(x) (baseline - notrans)/baseline * 100,
                            indirectperc = function(x) (notrans - x)/baseline * 100,
                            totalperc    = function(x) (baseline - x)/baseline * 100,
                            
                            directabs   = function(x) baseline - notrans,
                            indirectabs = function(x) notrans - x,
                            totalabs    = function(x) baseline - x),
                .names = "{.col}_{.fn}")) %>%
  # rename the original columns to follow the same pattern as _abs/_perc
  rename_with(~ paste0(.x, "_orig"), .cols = all_of(scenario_short[scenario_short != "baseline"])) %>% 
  # then gather
  pivot_longer(cols = contains("trans"), names_to   = c("scenario", ".value"), names_sep  = "_") %>%
  # then summarize
  group_by(severity, scenario) %>%
  get_summary(contains(c("orig", "direct", "indirect", "total"))) %>%
  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels))

# total effects: scenario analyses
p_sens_total <- 
  sens_H %>% filter(scenario != "baseline") %>%
  ggplot(aes(x = scenario_names, y = totalabs_mean, color = severity, group = severity)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) + 
  geom_errorbar(aes(ymin = totalabs_lcl, ymax = totalabs_ucl), position = position_dodge(width = 0.5), width = 0.3) +
  geom_text(aes(label = paste0(sprintf("%.1f", totalperc_mean), "%")), size = 3.5,
            position = position_dodge(width = 0.5), hjust= -0.3, color = "black") +
  get_theme(txt = 13) + 
  scale_x_discrete(labels = function(l) str_wrap(l, width = 10)) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual("Season\nseverity", values = cols_severity[c(1, 3)]) +
  labs(x = "Scenario", y = "Number of hospitalizations prevented")
ggsave(p_sens_total, file = here(figpath, "prevented_h_sens_total.pdf"), width = 7, height = 5)


# total effects: extended spectrum analyses
sens_pts <- 
  data.frame(scenario = c("No transmission impact", "Small transmission impact", "Large transmission impact"),
             av_hosp = c(0.9, 0.9, 0.9),
             av_trans = c(1, 0.9, 0.8),
             p_av = mean(pars_fixed$p_a1 + pars_fixed$p_a2),
             prevented = 0)

tmp_sens_plt <-
  out_spectrum %>% filter(combination == "av_trans x p_av sensitivity") %>%
  left_join(
    out_spectrum %>% filter(av_trans == 1, p_av == 0, combination == "av_trans x p_av sensitivity") %>% 
      rename(baseline_H = H_all) %>% select(severity, scenario, baseline_H)
  ) %>%
  mutate(severity = paste(str_to_title(severity), "severity")) %>%
  mutate(prevented = (baseline_H - H_all)/baseline_H * 100) 

p_sens_spec <- tmp_sens_plt %>%
  mutate(severity = paste(severity, "season")) %>%
  ggplot(aes(x = p_av*100, y = (1 - av_trans)*100, fill = prevented, color = prevented)) +
  geom_tile() + facet_wrap(~ severity) +
  geom_contour(aes(z = prevented), breaks = comp_lines2, 
               linetype = "dashed", colour = "white", linewidth = 0.4, alpha = 0.5) +
  stat_contour(aes(z = prevented, label = after_stat(level)), breaks = comp_lines2,
               geom = "text", color = "white", size = 2.5, check_overlap = TRUE, alpha = 0.5) +
  geom_point(data = sens_pts, col = "white") +
  geom_text(data = sens_pts, aes(label = str_wrap(scenario, width = 30)), 
            col = "white", hjust = -0.05, vjust = -0.2) +
  scale_fill_scico("Prevented\nhospitalizations\n(%)", palette = "lapaz") +
  scale_color_scico(guide = "none", palette = "lapaz") +
  get_theme(txt = 12.5) +
  labs(x = "Symptomatic individuals undergoing prompt treatment (%)", y = "Reduction in transmission risk (%)")
ggsave(p_sens_spec, file = here(figpath, "prevented_h_sens_spec.pdf"), width = 10, height = 5)


# quoted numbers - both 100%
tmp_sens_plt %>% filter(av_trans == 0, p_av == 1)


## End ---------------------


