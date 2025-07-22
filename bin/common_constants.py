PVAL_DSET = "p_value"
NEG_LOG_PVAL_DSET = "neg_log_10_p_value"
CHR_DSET = "chromosome"
BP_DSET = "base_pair_location"
OR_DSET = "odds_ratio"
RANGE_U_DSET = "ci_upper"
RANGE_L_DSET = "ci_lower"
BETA_DSET = "beta"
ZSCORE_DSET = "z_score"
SE_DSET = "standard_error"
RSID = "rsid"
EFFECT_DSET = "effect_allele"
OTHER_DSET = "other_allele"
FREQ_DSET = "effect_allele_frequency"
HM_CC_DSET = "hm_coordinate_conversion"
HM_OR_DSET = "hm_odds_ratio"
HM_RANGE_U_DSET = "hm_ci_upper"
HM_RANGE_L_DSET = "hm_ci_lower"
HM_BETA_DSET = "hm_beta"
HM_EFFECT_DSET = "hm_effect_allele"
HM_OTHER_DSET = "hm_other_allele"
HM_FREQ_DSET = "hm_effect_allele_frequency"
HM_CODE = "hm_code"

# Map each canonical column name to a list of sensible alternative names
COLUMN_ALIASES = {
    PVAL_DSET: [
        "p",
        "pval",
        "pvalue",
        "p-value",
        "P",
        "p_val",
        "p.value",
        "pv",
        "pval_nominal",
        "p_val_nominal",
        "p-val",
    ],
    CHR_DSET: [
        "chr",
        "chrom",
        "chromosome_name",
        "CHR",
        "CHROM",
        "chromosome_number",
        "chrom_no",
        "chromosome",
        "chromosome_no",
        "chromosomeid",
    ],
    BP_DSET: [
        "bp",
        "pos",
        "position",
        "BP",
        "BP_POS",
        "bp_position",
        "basepair",
        "base_pair",
        "coordinate",
        "bp_coord",
        "base_pair_pos",
        "bp_start",
    ],
    OR_DSET: [
        "or",
        "oddsratio",
        "OddsRatio",
        "OR",
        "odds_ratio",
        "OR_wald",
        "odds_ratio_wald",
    ],
    RANGE_U_DSET: [
        "or_upper",
        "ci_high",
        "upper_ci",
        "OR_U95",
        "ci_upper",
        "ci_upp",
        "or_upper_ci",
        "upper_ci95",
        "ci_u95",
        "CI_upper",
    ],
    RANGE_L_DSET: [
        "or_lower",
        "ci_low",
        "lower_ci",
        "OR_L95",
        "ci_lower",
        "ci_low95",
        "or_lower_ci",
        "lower_ci95",
        "ci_l95",
        "CI_lower",
    ],
    BETA_DSET: [
        "effect_size",
        "beta_estimate",
        "BETA",
        "beta_coef",
        "beta_coefficient",
        "beta",
        "slope",
        "effect",
        "beta_effect",
    ],
    ZSCORE_DSET: [
        "z",
        "zscore",
        "ZSCORE",
        "z_stat",
        "zstat",
        "z-statistic",
        "z_score",
        "zvalue",
        "z_statistic",
    ],
    SE_DSET: [
        "se",
        "stderr",
        "std_err",
        "se_beta",
        "sebeta",
        "se_beta_est",
    ],
    RSID: [
        "snp",
        "snpid",
        "markername",
        "variant_id",
        "rs_id",
        "variantid",
        "varid",
        "snp_id",
    ],
    EFFECT_DSET: [
        "ea",
        "alt",
        "effect_allele",
        "alt_allele",
        "ALT_ALLELE",
        "risk_allele",
    ],
    OTHER_DSET: [
        "non_effect_allele",
        "ref",
        "non_effect",
        "other",
        "other_allele",
        "reference_allele",
        "oa",
        "ref_allele",
    ],
    FREQ_DSET: ["eaf", "freq", "ALT_AF"],
    NEG_LOG_PVAL_DSET: ["neglog10p", "neg_log_p", "neg_log10_p", "neg_log_pvalue"],
}


# Reverse lookup table: alias -> canonical column name.
# is quite a long dictonary, i.e. {p: p_value, pval: p_value, pvalue: p_value,...}
ALIAS_LOOKUP = {
    alias.lower(): canonical
    for canonical, aliases in COLUMN_ALIASES.items()
    for alias in aliases + [canonical]
}


HARMONISER_ARG_MAP = {
    CHR_DSET: "--chrom_col",
    BP_DSET: "--pos_col",
    EFFECT_DSET: "--effAl_col",
    OTHER_DSET: "--otherAl_col",
    RSID: "--rsid_col",
    BETA_DSET: "--beta_col",
    OR_DSET: "--or_col",
    ZSCORE_DSET: "--zscore_col",
    RANGE_L_DSET: "--or_col_lower",
    RANGE_U_DSET: "--or_col_upper",
    FREQ_DSET: "--eaf_col",
    HM_CC_DSET: "--hm_coordinate_conversion",
}

STRAND_COUNT_ARG_MAP = {
    CHR_DSET: "--chrom_col",
    BP_DSET: "--pos_col",
    EFFECT_DSET: "--effAl_col",
    OTHER_DSET: "--otherAl_col",
    RSID: "--rsid_col",
    HM_CC_DSET: "--hm_coordinate_conversion",
}

DEFAULT_CHROMS = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "X",
    "Y",
]
