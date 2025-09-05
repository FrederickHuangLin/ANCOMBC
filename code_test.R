data(atlas1006, package = "microbiome")

# Subset to baseline
pseq = phyloseq::subset_samples(atlas1006, time == 0)

# Re-code the bmi group
meta_data = microbiome::meta(pseq)
meta_data$bmi = recode(meta_data$bmi_group,
                       obese = "obese",
                       severeobese = "obese",
                       morbidobese = "obese")

# Note that by default, levels of a categorical variable in R are sorted
# alphabetically. In this case, the reference level for `bmi` will be
# `lean`. To manually change the reference level, for instance, setting `obese`
# as the reference level, use:
meta_data$bmi = factor(meta_data$bmi, levels = c("obese", "overweight", "lean"))
# You can verify the change by checking:
# levels(meta_data$bmi)

# Create the region variable
meta_data$region = recode(as.character(meta_data$nationality),
                          Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE",
                          CentralEurope = "CE", EasternEurope = "EE",
                          .missing = "unknown")

phyloseq::sample_data(pseq) = meta_data

# Subset to lean, overweight, and obese subjects
pseq = phyloseq::subset_samples(pseq, bmi %in% c("lean", "overweight", "obese"))
# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
pseq = phyloseq::subset_samples(pseq, ! region %in% c("EE", "unknown"))

print(pseq)

set.seed(123)
# It should be noted that we have set the number of bootstrap samples (B) equal
# to 10 in the 'trend_control' function for computational expediency.
# However, it is recommended that users utilize the default value of B,
# which is 100, or larger values for optimal performance.
output = ancombc3_daa(data = pseq, tax_level = "Family",
                      fix_formula = "age + region + bmi", rand_formula = NULL,
                      p_adj_method = "holm", pseudo_sens = TRUE,
                      prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                      group = "bmi", alpha = 0.05, n_cl = 2, verbose = TRUE,
                      global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                      iter_control = list(tol = 1e-2, max_iter = 20,
                                          verbose = TRUE),
                      em_control = list(tol = 1e-5, max_iter = 100),
                      lme_control = lme4::lmerControl(),
                      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                      trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                  nrow = 2,
                                                                  byrow = TRUE),
                                                           matrix(c(-1, 0, 1, -1),
                                                                  nrow = 2,
                                                                  byrow = TRUE),
                                                           matrix(c(1, 0, 1, -1),
                                                                  nrow = 2,
                                                                  byrow = TRUE)),
                                           node = list(2, 2, 1),
                                           solver = "ECOS",
                                           B = 10))

res_prim = output$res

df_age = res_prim %>%
    dplyr::select(taxon, ends_with("age"))
df_fig_age = df_age %>%
    dplyr::filter(diff_age == 1) %>%
    dplyr::arrange(desc(lfc_age)) %>%
    dplyr::mutate(direct = ifelse(lfc_age > 0, "Positive LFC", "Negative LFC"),
                  color = ifelse(diff_robust_age, "aquamarine3", "black"))
df_fig_age$taxon = factor(df_fig_age$taxon, levels = df_fig_age$taxon)
df_fig_age$direct = factor(df_fig_age$direct,
                           levels = c("Positive LFC", "Negative LFC"))

fig_age = df_fig_age %>%
    ggplot(aes(x = taxon, y = lfc_age, fill = direct)) +
    geom_bar(stat = "identity", width = 0.7, color = "black",
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lfc_age - se_age, ymax = lfc_age + se_age),
                  width = 0.2, position = position_dodge(0.05), color = "black") +
    labs(x = NULL, y = "Log fold change",
         title = "Log fold changes as one unit increase of age") +
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1,
                                     color = df_fig_age$color))
fig_age

res_trend = output$res_trend

df_fig_trend = res_trend %>%
    dplyr::filter(diff_abn == 1) %>%
    dplyr::mutate(lfc1 = round(lfc_bmioverweight, 2),
                  lfc2 = round(lfc_bmilean, 2),
                  color = ifelse(diff_robust_abn, "aquamarine3", "black")) %>%
    tidyr::pivot_longer(cols = lfc1:lfc2,
                        names_to = "group", values_to = "value") %>%
    dplyr::arrange(taxon)

df_fig_trend$group = recode(df_fig_trend$group,
                            `lfc1` = "Overweight - Obese",
                            `lfc2` = "Lean - Obese")
df_fig_trend$group = factor(df_fig_trend$group,
                            levels = c("Overweight - Obese",
                                       "Lean - Obese"))

lo = floor(min(df_fig_trend$value))
up = ceiling(max(df_fig_trend$value))
mid = (lo + up)/2
fig_trend = df_fig_trend %>%
    ggplot(aes(x = group, y = taxon, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         na.value = "white", midpoint = mid, limit = c(lo, up),
                         name = NULL) +
    geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
    labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(color = df_fig_trend %>%
                                         dplyr::distinct(taxon, color) %>%
                                         .$color))
fig_trend
