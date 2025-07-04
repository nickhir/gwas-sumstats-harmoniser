PVAL_DSET = "p_value"
CHR_DSET = "chromosome"
BP_DSET = "base_pair_location"
OR_DSET = "odds_ratio"
RANGE_U_DSET = "ci_upper"
RANGE_L_DSET = "ci_lower"
BETA_DSET = "beta"
ZSCORE_DSET = "z_score"
RSID = "variant_id"
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
HM_VAR_ID = "hm_variant_id"
HM_CODE = "hm_code"


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
    "MT",
]
