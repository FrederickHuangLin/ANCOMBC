# =============================================================================
# ANCOM-BC3 omnibus inference and component labels
# =============================================================================
#
# omnibus_pvalue() combines the DA and DP P values of the taxa of the omnibus
# universe U, the taxa eligible for DA testing, DP testing, or both: for a
# taxon eligible for both, the Bonferroni-adjusted minimum
# min{1, 2 min(P_DA, P_DP)}; for a taxon eligible for one component, its
# component P value. ancombc3_component_labels() labels each omnibus discovery
# by the component P values that support it.

# Per-taxon vector v in the order of the taxa taxa, the names of has_da: by
# name when both v and has_da are named, by position otherwise (v must then
# have p elements). The result is unnamed.
.ac3_omnibus_match <- function(v, taxa, p, name) {
  if (!is.null(taxa) && !is.null(names(v))) {
    absent <- setdiff(taxa, names(v))
    if (length(absent) > 0L) {
      stop(name, " has no value for ", length(absent), " taxon/taxa of has_da: ",
           paste(utils::head(absent, 5L), collapse = ", "), ".")
    }
    return(unname(v[taxa]))
  }
  if (length(v) != p) stop(name, " must have one value per taxon of has_da.")
  unname(v)
}

# Omnibus P values of one tested covariate. has_da and has_dp flag the taxa
# with a DA and a DP P value; a taxon with neither receives NA. The taxa are
# those of has_da: when has_da is named, every named argument is matched to
# its names; unnamed arguments are matched by position. The result is unnamed,
# in the order of has_da.
omnibus_pvalue <- function(p_da, p_dp, has_da, has_dp) {
  p <- length(has_da)
  taxa <- names(has_da)
  p_da <- .ac3_omnibus_match(p_da, taxa, p, "p_da")
  p_dp <- .ac3_omnibus_match(p_dp, taxa, p, "p_dp")
  has_dp <- .ac3_omnibus_match(has_dp, taxa, p, "has_dp")
  has_da <- unname(has_da)
  out <- rep(NA_real_, p)
  both <- has_da & has_dp
  only_da <- has_da & !has_dp
  only_dp <- !has_da & has_dp
  if (any(both)) out[both] <- pmin(1, 2 * pmin(p_da[both], p_dp[both]))
  out[only_da] <- p_da[only_da]
  out[only_dp] <- p_dp[only_dp]
  out
}

# Component labels of the taxa of one tested covariate. With R omnibus
# discoveries at level alpha among the |U| taxa of the universe, the component
# threshold tau is half the rejection cutoff of the omnibus P values under the
# multiplicity adjustment p_adj_method:
#   "BH" (and NULL, the default): R alpha / |U|;
#   "BY": R alpha / (|U| c(|U|)), with c(m) = sum_{i <= m} 1 / i;
#   "bonferroni": alpha / |U|;
#   "holm" and "hochberg": alpha / (|U| - R + 1), the critical value of the
#     last rejected hypothesis;
#   "none": alpha.
# A discovery eligible for both components is labeled
# "DA_driven" (P_DA <= tau < P_DP), "DP_driven" (P_DP <= tau < P_DA), "DA_DP"
# (both at most tau), or "omnibus_only" (both above tau); a discovery eligible
# for one component is labeled by that component. Taxa of the universe that are
# not discovered are labeled "not_discovered", and taxa outside the universe
# "not_in_universe".
# The taxa are matched as in omnibus_pvalue(): by the names of has_da when it
# is named, by position otherwise; the result is unnamed, in the order of
# has_da.
ancombc3_component_labels <- function(discovered, has_da, has_dp, p_da, p_dp,
                                      alpha = 0.1, p_adj_method = NULL) {
  p <- length(has_da)
  taxa <- names(has_da)
  discovered <- .ac3_omnibus_match(discovered, taxa, p, "discovered")
  has_dp <- .ac3_omnibus_match(has_dp, taxa, p, "has_dp")
  p_da <- .ac3_omnibus_match(p_da, taxa, p, "p_da")
  p_dp <- .ac3_omnibus_match(p_dp, taxa, p, "p_dp")
  has_da <- unname(has_da)
  in_universe <- has_da | has_dp
  label <- rep("not_in_universe", length(discovered))
  label[in_universe & !discovered] <- "not_discovered"

  only_da <- discovered & has_da & !has_dp
  only_dp <- discovered & has_dp & !has_da
  both <- discovered & has_da & has_dp
  label[only_da] <- "DA_driven"
  label[only_dp] <- "DP_driven"

  n_universe <- sum(in_universe)
  n_discoveries <- sum(discovered)
  tau <- if (n_universe > 0L) n_discoveries * alpha / (2 * n_universe) else 0
  if (!is.null(p_adj_method) && !identical(p_adj_method, "BH") && n_universe > 0L) {
    tau <- switch(
      p_adj_method,
      BY = n_discoveries * alpha / (n_universe * sum(1 / seq_len(n_universe))) / 2,
      bonferroni = alpha / n_universe / 2,
      holm = ,
      hochberg = alpha / (n_universe - n_discoveries + 1) / 2,
      none = alpha / 2,
      stop("p_adj_method must be one of: BH, BY, holm, hochberg, bonferroni, none.")
    )
  }
  da_pass <- has_da & is.finite(p_da) & p_da <= tau
  dp_pass <- has_dp & is.finite(p_dp) & p_dp <= tau
  label[both] <- "omnibus_only"
  label[both & da_pass & !dp_pass] <- "DA_driven"
  label[both & dp_pass & !da_pass] <- "DP_driven"
  label[both & da_pass & dp_pass] <- "DA_DP"
  label
}
