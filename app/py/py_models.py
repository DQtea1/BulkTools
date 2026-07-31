# contains : 
# evaluate_signature_threshold_cv() : Evaluates a signature for a specific model
# sample_confidence_report() : Evaluates the confidence we can have in a specific prediction


import os
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, confusion_matrix
from joblib import Parallel, delayed


# Ceiling on replicates for ci_method="pipeline". Each replicate re-runs the
# whole CV, so an uncapped n_bootstrap=1000 means ~3.5 minutes for an opt-in
# mode; 300 costs ~65 s on 8 threads (measured, ~0.21 s per replicate wall).
#
# Measured on nested prefixes of 400 replicates (17R/23NR), largest shift in a
# 95% bound between n=300 and n=400:
#   sensitivity 0.000   balanced_accuracy 0.001   auc 0.008   coverage 0.000
#   specificity 0.021   ppv 0.026
# The two larger ones are NOT slow Monte Carlo convergence: they are the upper
# bound flipping between the top two attainable values (0.974 vs 1.000) because
# with n=40 these metrics are discrete and their distribution piles up at 1.0.
# More replicates will not resolve that -- only more patients would. So 300 is
# enough for anything you would quote.
#
# Override deliberately with  py_models._PIPELINE_MAX_BOOTSTRAP = N.
_PIPELINE_MAX_BOOTSTRAP = 300


def _n_jobs():
    """Parallel worker budget. Honors SHINY_N_CORES (set by the launchers to
    ~3/4 of the host cores); otherwise falls back to 3/4 of os.cpu_count().
    Always >= 1. Bootstrap loops below use a thread backend, so this never
    pessimizes: worst case it behaves like the serial loop."""
    env = os.environ.get("SHINY_N_CORES", "")
    try:
        n = int(env)
        if n >= 1:
            return n
    except (ValueError, TypeError):
        pass
    cpu = os.cpu_count() or 1
    return max(1, int(cpu * 0.75))


# -----------------------------------------------------------------------------
# Shared helpers
#
# Used by BOTH evaluate_signature_threshold_cv() and sample_confidence_report().
# They used to exist as one private copy inside each function, which let the two
# drift apart (different kwarg names, different defaults, different tie handling).
# Keeping a single definition means the two entry points can only ever apply
# the same rule.
# -----------------------------------------------------------------------------

def _safe_div(num, den):
    return np.nan if den == 0 else num / den


def _percentile_ci(values, ci_level=0.95):
    # Drop non-finite, not just NaN. roc_curve's thresholds[0] is +inf in
    # sklearn >= 1.3, and a degenerate bootstrap replicate (no finite threshold
    # reaching the target specificity catches any responder) selects it, so
    # t_high_boot really does contain +inf on imbalanced cohorts. np.quantile
    # then does inf - inf and returns nan with a RuntimeWarning, turning the
    # reported CI into garbage. The consumers of these arrays already mask with
    # np.isfinite; this makes the CI agree with them. n_bootstrap_valid_* in the
    # summary reports how many replicates survived.
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (np.nan, np.nan)
    alpha = 1 - ci_level
    lo = np.quantile(arr, alpha / 2)
    hi = np.quantile(arr, 1 - alpha / 2)
    return float(lo), float(hi)


def _pick_single_threshold_from_roc(y_true, scores, criterion="youden"):
    """
    scores oriented so that higher = more likely responder.
    Returns chosen threshold + ROC arrays.
    """
    fpr, tpr, thresholds = roc_curve(y_true, scores)

    # roc_curve prepends thresholds[0] = +inf (the "call nothing positive"
    # operating point). It is not a usable cut, and on a degenerate fold where
    # tpr - fpr <= 0 everywhere nanargmax returns index 0 and hands back +inf,
    # which then poisons median(fold_thresholds). Score only finite candidates.
    finite = np.isfinite(thresholds)
    if not finite.any():
        raise ValueError("roc_curve returned no finite threshold.")

    if criterion == "youden":
        # J = sensitivity + specificity -1 = tpr - fpr
        j = np.where(finite, tpr - fpr, -np.inf)
        idx = int(np.nanargmax(j))
    elif criterion == "closest_topleft":
        d = np.where(finite, np.sqrt((fpr - 0.0) ** 2 + (tpr - 1.0) ** 2), np.inf)
        idx = int(np.nanargmin(d))
    else:
        raise ValueError("threshold_criterion must be 'youden' or 'closest_topleft'.")

    return float(thresholds[idx]), {
        "fpr": fpr,
        "tpr": tpr,
        "thresholds": thresholds,
        "best_idx": idx,
    }


def _pick_gray_thresholds_from_roc(y_true, scores, target_sens=0.90, target_spec=0.90):
    """
    Returns (t_low, t_high) based on ROC on oriented scores (higher => responder):
      - t_low: highest threshold achieving target sensitivity ON THE TRAINING
               ROC -> scores < t_low are called "confident negatives"
      - t_high: lowest threshold achieving target specificity ON THE TRAINING
               ROC -> scores >= t_high are called "confident positives"

    Both are boundary choices: the most aggressive threshold that still meets
    the constraint, hence the most optimistic. The target is an in-sample
    constraint and NOT a guarantee about new samples -- realised out-of-fold
    sensitivity below t_low runs ~0.02-0.07 short of it at cohort sizes of
    40-50. gray_zone["realised"] reports what was actually achieved; use that,
    not the target, when describing the zone.

    Labels are assigned by _gray_zone_labels(), which is the single place the
    rule lives; see it for the tie convention and the inverted-zone case.
    """
    fpr, tpr, thresholds = roc_curve(y_true, scores)
    spec = 1 - fpr

    # Finite candidates only -- see _pick_single_threshold_from_roc. thresholds[0]
    # is +inf and always satisfies "spec >= target", so without this filter a
    # degenerate replicate (no finite high-specificity threshold catches any
    # responder) returns t_high = +inf. Observed live on an 8R/32NR cohort.
    finite = np.isfinite(thresholds)
    if not finite.any():
        raise ValueError("roc_curve returned no finite threshold.")
    idx_finite = np.where(finite)[0]

    # t_low from target sensitivity
    idx_sens_valid = np.where((tpr >= target_sens) & finite)[0]
    if len(idx_sens_valid) == 0:
        # Fallback: pick threshold with max sensitivity (usually very low threshold)
        idx_low = idx_finite[int(np.nanargmax(tpr[idx_finite]))]
    else:
        # Among thresholds with target sensitivity, maximize specificity
        idx_low = idx_sens_valid[int(np.nanargmax(spec[idx_sens_valid]))]

    # t_high from target specificity
    idx_spec_valid = np.where((spec >= target_spec) & finite)[0]
    if len(idx_spec_valid) == 0:
        # Fallback: pick threshold with max specificity (usually very high threshold)
        idx_high = idx_finite[int(np.nanargmax(spec[idx_finite]))]
    else:
        # Among thresholds with target specificity, maximize sensitivity
        idx_high = idx_spec_valid[int(np.nanargmax(tpr[idx_spec_valid]))]

    t_low = float(thresholds[idx_low])
    t_high = float(thresholds[idx_high])

    # In rare/noisy cases t_low > t_high (gray zone "negative width")
    # We keep values + flag upstream, but can also collapse if needed.
    return t_low, t_high, {
        "fpr": fpr,
        "tpr": tpr,
        "spec": spec,
        "thresholds": thresholds,
        "idx_low": idx_low,
        "idx_high": idx_high,
    }


def _gray_zone_labels(scores, t_low, t_high):
    """
    Canonical gray-zone rule. Returns labels in {-1, 0, 1}:
        1  confident responder      : score >= t_high
        0  confident non-responder  : score <  t_low
       -1  uncertain                : everything else

    The asymmetry (< for t_low, >= for t_high) is deliberate. roc_curve computes
    tpr[i] under "score >= thresholds[i]", so the target-sensitivity guarantee
    attached to t_low counts a sample sitting exactly on t_low as a caught
    responder. The old rule ("score <= t_low" => confident negative) swallowed
    those samples and so broke the very guarantee t_low was selected to satisfy.
    t_high needs no such correction: specificity counts negatives with
    score < t_high, so "score >= t_high" is already the matching rule.

    Inverted zone (t_low > t_high): the two criteria overlap on [t_high, t_low)
    and both fire. They disagree, so those samples are labelled uncertain. The
    two former copies of this rule resolved that overlap in OPPOSITE directions
    (one let positives overwrite negatives, the other returned negative first),
    which meant the same score could come back 0 or 1 depending on which
    function you asked. Marking the conflict uncertain also keeps a degenerate
    zone visible in the coverage figure instead of silently pushing it to 100%.

    Broadcasts: scores and thresholds may each be scalar or array.
    """
    s = np.asarray(scores)
    is_neg = s < t_low
    is_pos = s >= t_high
    return np.where(is_pos & ~is_neg, 1, np.where(is_neg & ~is_pos, 0, -1)).astype(int)


def _recentre_prob(p, p_operating):
    """
    Rescale a calibrated probability so that the DECISION BOUNDARY sits at 0.5.

    The binary call is "score >= thr", which under a monotone calibration is
    exactly "p >= p_operating", where p_operating = P(responder | score = thr).
    Youden's threshold does not generally land on p = 0.5: it drifts with class
    prevalence. Measured on synthetic cohorts, p_operating is ~0.47 at 25R/25NR,
    ~0.38-0.54 at 17R/23NR, ~0.20 at 8R/32NR and ~0.84 at 30R/10NR. Anything
    calling itself "confidence" while measuring distance from 0.5 is therefore
    measuring from the wrong point, and on an imbalanced cohort it will report a
    sample sitting right on the decision boundary as highly confident.

    Maps [0, p_operating, 1] -> [0, 0.5, 1], linear on each side. Monotone, so
    the ranking of samples is unchanged; only the zero point moves.

    Broadcasts: p and p_operating may each be scalar or array.
    """
    p = np.asarray(p, dtype=float)
    p0 = np.asarray(p_operating, dtype=float)
    eps = 1e-12
    below = 0.5 * p / np.maximum(p0, eps)
    above = 0.5 + 0.5 * (p - p0) / np.maximum(1.0 - p0, eps)
    return np.clip(np.where(p >= p0, above, below), 0.0, 1.0)


def evaluate_signature_threshold_cv(
    df_scores: pd.DataFrame,
    sample_ID_responders,
    sample_ID_non_responders,
    score_col: str = None,
    # CV
    n_folds: int = 5,
    n_subsamples: int = 5,
    shuffle: bool = True,
    random_state: int = 42,
    # Threshold selection (single threshold)
    threshold_criterion: str = "youden",  # "youden" | "closest_topleft"
    higher_score_is_responder: bool = None,
    # Gray zone
    use_gray_zone: bool = True,
    gray_target_sensitivity: float = 0.90,   # for t_low
    gray_target_specificity: float = 0.90,   # for t_high
    # Bootstrap CIs
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    bootstrap_seed: int = 123,
    bootstrap_stratified: bool = True,
    ci_method: str = "conditional",   # "conditional" (fast) | "pipeline" (honest, slow,
                                      #   replicates capped at _PIPELINE_MAX_BOOTSTRAP)
    # Outputs
    return_oof_predictions: bool = True,
):
    """
    Cross-validated decision threshold for a fixed per-sample score (1D):
      - single-threshold classification, the threshold being ROC-derived on each
        fold's training data and bagged over n_subsamples subsamples of it
      - optional gray zone (double threshold)
      - confusion matrices, Se/Sp/PPV/NPV, ROC/AUC
      - bootstrap confidence intervals, CONDITIONAL on the fitted threshold
        (see results["bootstrap_ci_interpretation"]); the gray-zone coverage
        and rejection_rate intervals are the exception and refit the thresholds
        on every replicate

    No model is fitted on features: the score arrives precomputed and the only
    thing estimated here is where to cut it.

    Formerly called nested_cv_signature, which was wrong twice over. There is no
    nesting: nothing is ever scored on held-out inner data and no model or
    hyper-parameter is selected, so there is no inner validation to nest -- the
    inner split only generates subsamples for bagging the threshold. And the
    cross-validation protects less than the name suggested: the score is computed
    upstream without using the labels, so the ROC/AUC reported here is a
    resubstitution ROC over the whole cohort. Only the threshold-dependent
    metrics (Se/Sp/PPV/NPV) are genuinely cross-validated.

    Inputs
    ------
    df_scores : pd.DataFrame
        Index = sample IDs. Contains one score column (or specify score_col).
    sample_ID_responders : list-like
        IDs for positive class (responder = 1).
    sample_ID_non_responders : list-like
        IDs for negative class (non-responder = 0).

    Returns
    -------
    results : dict
        Contains ROC, confusion matrices, metrics, bootstrap CIs, gray zone outputs,
        thresholds per fold, and OOF predictions.
    """

    # -------------------------------------------------------------------------
    # Helpers
    # -------------------------------------------------------------------------
    def _compute_binary_metrics(y_true, y_pred):
        """
        y_true, y_pred in {0,1}
        """
        cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()

        sens = _safe_div(tp, tp + fn)  # recall positive
        spec = _safe_div(tn, tn + fp)
        ppv = _safe_div(tp, tp + fp)
        npv = _safe_div(tn, tn + fn)

        acc = _safe_div(tp + tn, tp + tn + fp + fn)
        bal_acc = np.nan if (np.isnan(sens) or np.isnan(spec)) else (sens + spec) / 2
        fpr = np.nan if np.isnan(spec) else 1 - spec
        fnr = np.nan if np.isnan(sens) else 1 - sens
        prevalence = _safe_div(np.sum(y_true == 1), len(y_true))

        # Risks conditional on predicted label
        risk_fp_given_pred_pos = np.nan if np.isnan(ppv) else 1 - ppv
        risk_fn_given_pred_neg = np.nan if np.isnan(npv) else 1 - npv

        out = {
            "n": int(len(y_true)),
            "tn": int(tn),
            "fp": int(fp),
            "fn": int(fn),
            "tp": int(tp),
            "sensitivity": sens,
            "specificity": spec,
            "ppv": ppv,
            "npv": npv,
            "fpr": fpr,
            "fnr": fnr,
            "accuracy": acc,
            "balanced_accuracy": bal_acc,
            "prevalence": prevalence,
            "risk_false_positive_given_pred_positive": risk_fp_given_pred_pos,
            "risk_false_negative_given_pred_negative": risk_fn_given_pred_neg,
        }
        return out, cm

    def _bootstrap_ci_binary_metrics(
        y_true,
        y_pred,
        y_score=None,
        n_boot=1000,
        ci_level=0.95,
        seed=123,
        stratified=True
    ):
        """
        Bootstrap percentile CIs for binary metrics (and AUC if y_score provided).
        """
        y_true = np.asarray(y_true, dtype=int)
        y_pred = np.asarray(y_pred, dtype=int)
        y_score = None if y_score is None else np.asarray(y_score, dtype=float)

        n = len(y_true)
        idx_all = np.arange(n)
        idx_pos = np.where(y_true == 1)[0]
        idx_neg = np.where(y_true == 0)[0]

        metric_names = [
            "sensitivity", "specificity", "ppv", "npv",
            "fpr", "fnr", "accuracy", "balanced_accuracy",
            "risk_false_positive_given_pred_positive",
            "risk_false_negative_given_pred_negative",
        ]
        samples = {m: [] for m in metric_names}
        if y_score is not None:
            samples["auc"] = []

        # One bootstrap replicate -> dict of metric values (+ auc). Per-replicate
        # seeds keep results independent of execution order / worker count.
        def _one_metric_boot(child_seed):
            rng_b = np.random.default_rng(child_seed)
            if stratified:
                if len(idx_pos) == 0 or len(idx_neg) == 0:
                    boot_idx = rng_b.choice(idx_all, size=n, replace=True)
                else:
                    boot_pos = rng_b.choice(idx_pos, size=len(idx_pos), replace=True)
                    boot_neg = rng_b.choice(idx_neg, size=len(idx_neg), replace=True)
                    boot_idx = np.concatenate([boot_pos, boot_neg])
                    rng_b.shuffle(boot_idx)
            else:
                boot_idx = rng_b.choice(idx_all, size=n, replace=True)

            yb = y_true[boot_idx]
            pb = y_pred[boot_idx]
            mets, _ = _compute_binary_metrics(yb, pb)
            row = {m: mets[m] for m in metric_names}
            if y_score is not None:
                sb = y_score[boot_idx]
                if len(np.unique(yb)) < 2:
                    row["auc"] = np.nan
                else:
                    fpr_b, tpr_b, _ = roc_curve(yb, sb)
                    row["auc"] = auc(fpr_b, tpr_b)
            return row

        child_seeds = np.random.SeedSequence(seed).spawn(n_boot)
        _n = _n_jobs()
        if _n == 1:
            rows = [_one_metric_boot(s) for s in child_seeds]
        else:
            try:
                rows = Parallel(n_jobs=_n, prefer="threads")(
                    delayed(_one_metric_boot)(s) for s in child_seeds
                )
            except Exception:
                rows = [_one_metric_boot(s) for s in child_seeds]

        for row in rows:
            for m in metric_names:
                samples[m].append(row[m])
            if y_score is not None:
                samples["auc"].append(row["auc"])

        out = {}
        for m, vals in samples.items():
            lo, hi = _percentile_ci(vals, ci_level=ci_level)
            out[m] = {"ci_low": lo, "ci_high": hi}
        return out

    # Which metrics carry cross-validation protection and which do not. auc/roc
    # are computed from the raw score over the whole cohort (see the note in
    # section 3); everything derived from the confusion matrix comes from
    # out-of-fold predictions. Tabulating them together without this distinction
    # is what made the AUC look cross-validated.
    _VALIDATION = {
        "auc": "resubstitution (whole cohort, no CV protection)",
        "prevalence": "cohort property (not a performance metric)",
        "coverage": "out-of-fold",
        "rejection_rate": "out-of-fold",
    }

    def _metrics_ci_table(point_estimates: dict, ci_dict: dict, metric_order=None,
                          population=None, validation=True):
        """
        population: string naming WHICH samples these metrics describe. The
        confident-subset table used to be schema-identical to the whole-cohort
        one, so two tables describing different populations could be read side
        by side as a like-for-like comparison -- and the confident subset always
        looks better, because the hard cases were removed from it.
        """
        rows = []
        if metric_order is None:
            metric_order = list(point_estimates.keys())
        for m in metric_order:
            if m not in point_estimates:
                continue
            pe = point_estimates[m]
            ci = ci_dict.get(m, {"ci_low": np.nan, "ci_high": np.nan})
            row = {
                "metric": m,
                "estimate": pe,
                "ci_low": ci["ci_low"],
                "ci_high": ci["ci_high"],
            }
            if validation:
                row["validation"] = _VALIDATION.get(m, "out-of-fold")
            if population is not None:
                row["population"] = population
            rows.append(row)
        cols = ["metric", "estimate", "ci_low", "ci_high"]
        if validation:
            cols.append("validation")
        if population is not None:
            cols.append("population")
        return pd.DataFrame(rows, columns=cols)

    # -------------------------------------------------------------------------
    # 1) Validate inputs, build labels, orient the score
    # -------------------------------------------------------------------------
    if not isinstance(df_scores, pd.DataFrame):
        raise TypeError("df_scores must be a pandas DataFrame.")
    if df_scores.shape[1] < 1:
        raise ValueError("df_scores must contain at least one column.")
    if score_col is None:
        score_col = df_scores.columns[0]
    if score_col not in df_scores.columns:
        raise ValueError(f"score_col '{score_col}' not found in df_scores.")
    if not df_scores.index.is_unique:
        raise ValueError("df_scores.index must contain unique sample IDs.")

    responders = set(sample_ID_responders)
    non_responders = set(sample_ID_non_responders)

    overlap = responders & non_responders
    if len(overlap) > 0:
        ex = sorted(overlap, key=str)[:10]
        raise ValueError(f"Some IDs are in both classes (e.g., {ex}).")

    # Canonical (sorted) row order. Iteration order of a set depends on
    # PYTHONHASHSEED, which Python randomises per process, so building all_ids
    # straight from the union reshuffled df on every launch; since
    # StratifiedKFold(shuffle=True) permutes positions, the folds, thresholds
    # and every downstream metric silently changed between runs despite
    # random_state being fixed. Sorting also makes the result depend only on
    # WHICH samples are passed, not on the order the caller happened to pass
    # them in (upstream merges/filters in R can reorder the ID vectors).
    all_ids = sorted(responders | non_responders, key=str)
    missing = [sid for sid in all_ids if sid not in df_scores.index]
    if len(missing) > 0:
        raise ValueError(f"{len(missing)} IDs not found in df_scores.index (e.g., {missing[:10]}).")

    df = df_scores.loc[all_ids, [score_col]].copy()
    df["y"] = 0
    df.loc[df.index.isin(list(responders)), "y"] = 1
    df["score_raw"] = pd.to_numeric(df[score_col], errors="coerce")

    if df["score_raw"].isna().any():
        bad = df.index[df["score_raw"].isna()].tolist()[:10]
        raise ValueError(f"NaN/non-numeric scores detected (e.g., {bad}).")

    # infer orientation if needed
    if higher_score_is_responder is None:
        mean_pos = df.loc[df["y"] == 1, "score_raw"].mean()
        mean_neg = df.loc[df["y"] == 0, "score_raw"].mean()
        higher_score_is_responder = bool(mean_pos >= mean_neg)

    if higher_score_is_responder:
        df["score"] = df["score_raw"].values
        orientation_txt = "higher score => responder"
    else:
        df["score"] = -df["score_raw"].values
        orientation_txt = "lower raw score => responder (internally inverted score)"

    score_oriented = df["score"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=int)
    sample_ids = df.index.to_numpy()

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos < 2 or n_neg < 2:
        raise ValueError("Need at least 2 samples per class.")

    n_folds_eff = min(n_folds, n_pos, n_neg)
    if n_folds_eff < 2:
        raise ValueError("Effective n_folds < 2 (not enough samples per class).")

    # -------------------------------------------------------------------------
    # 2) Cross-validation: one threshold per fold, applied to its held-out samples
    # -------------------------------------------------------------------------
    fold_cv = StratifiedKFold(
        n_splits=int(n_folds_eff), shuffle=shuffle, random_state=random_state
    )

    # OOF containers
    # Every sample must be labelled by exactly one fold. Tracked explicitly
    # rather than by refilling a copy of the score array: that copy ended up
    # bit-identical to the input score, which is precisely why the ROC computed
    # from it is resubstitution and not out-of-fold (see section 3).
    fold_covered = np.zeros(len(y), dtype=bool)
    oof_pred = np.full(len(y), -1, dtype=int)         # standard binary pred
    oof_thr = np.full(len(y), np.nan, dtype=float)    # threshold used in that fold

    # Gray zone OOF containers
    if use_gray_zone:
        oof_pred_gray = np.full(len(y), -999, dtype=int)  # -1 uncertain, 0,1 valid
        oof_t_low = np.full(len(y), np.nan, dtype=float)
        oof_t_high = np.full(len(y), np.nan, dtype=float)
    else:
        oof_pred_gray = None
        oof_t_low = None
        oof_t_high = None

    fold_details = []
    fold_thresholds = []
    fold_gray_thresholds = []

    for fold_id, (train_idx, test_idx) in enumerate(fold_cv.split(score_oriented, y), start=1):
        score_train, y_train = score_oriented[train_idx], y[train_idx]
        score_test, y_test = score_oriented[test_idx], y[test_idx]

        # THRESHOLD BAGGING (deliberately not nested CV).
        #
        # The subsample split is used only to draw subsamples of this fold's
        # training set; a threshold is fitted on each subsample and the median
        # is taken. No validation happens on those subsamples and none is
        # intended -- there is no model and no hyper-parameter to select, so
        # there is nothing to validate. Hence the validation index returned by
        # the splitter is discarded below.
        #
        # Measured on 9 synthetic cohorts x 10 CV partitions, against fitting a
        # single threshold on the full fold training set:
        #   mean balanced accuracy   0.7001 (bagged) vs 0.6995 (single) - a wash
        #   run-to-run SD of bacc    0.01345 vs 0.01432  -> bagging ~6% tighter
        #   run-to-run SD of spec    0.02483 vs 0.02687  -> bagging ~8% tighter
        # So the bagging earns a small stability gain and is kept; what was
        # wrong was calling it nesting.
        n_pos_train = int((y_train == 1).sum())
        n_neg_train = int((y_train == 0).sum())
        n_subsamples_eff = min(n_subsamples, n_pos_train, n_neg_train)

        subsample_thresholds_single = []
        subsample_thresholds_low = []
        subsample_thresholds_high = []

        threshold_fit_mode = "bagged_over_subsamples"
        if n_subsamples_eff < 2:
            # Too few samples per class to subsample: fit once on fold train.
            threshold_fit_mode = "single_fit_fold_train"
            thr_single, _ = _pick_single_threshold_from_roc(y_train, score_train, criterion=threshold_criterion)
            subsample_thresholds_single = [thr_single]

            if use_gray_zone:
                t_low, t_high, _ = _pick_gray_thresholds_from_roc(
                    y_train, score_train,
                    target_sens=gray_target_sensitivity,
                    target_spec=gray_target_specificity
                )
                subsample_thresholds_low = [t_low]
                subsample_thresholds_high = [t_high]
        else:
            subsample_cv = StratifiedKFold(
                n_splits=int(n_subsamples_eff),
                shuffle=shuffle,
                random_state=random_state + fold_id
            )
            # The validation index is intentionally discarded: the split is a
            # subsampler here, not a validation scheme (see the note above).
            for subsample_idx, _unused_val_idx in subsample_cv.split(score_train, y_train):
                score_subsample = score_train[subsample_idx]
                yi_train = y_train[subsample_idx]

                thr_i, _ = _pick_single_threshold_from_roc(
                    yi_train, score_subsample, criterion=threshold_criterion
                )
                subsample_thresholds_single.append(thr_i)

                if use_gray_zone:
                    t_low_i, t_high_i, _ = _pick_gray_thresholds_from_roc(
                        yi_train, score_subsample,
                        target_sens=gray_target_sensitivity,
                        target_spec=gray_target_specificity
                    )
                    subsample_thresholds_low.append(t_low_i)
                    subsample_thresholds_high.append(t_high_i)

        # Bag the subsample thresholds -> this fold's threshold
        thr_fold = float(np.median(subsample_thresholds_single))
        fold_thresholds.append(thr_fold)

        # Apply standard binary threshold on the held-out fold
        y_pred_test = (score_test >= thr_fold).astype(int)

        # Store OOF
        fold_covered[test_idx] = True
        oof_pred[test_idx] = y_pred_test
        oof_thr[test_idx] = thr_fold

        # Gray zone
        gray_info = {}
        if use_gray_zone:
            t_low_fold = float(np.median(subsample_thresholds_low))
            t_high_fold = float(np.median(subsample_thresholds_high))

            pred_gray_test = _gray_zone_labels(score_test, t_low_fold, t_high_fold)

            oof_pred_gray[test_idx] = pred_gray_test
            oof_t_low[test_idx] = t_low_fold
            oof_t_high[test_idx] = t_high_fold

            fold_gray_thresholds.append((t_low_fold, t_high_fold))
            gray_info = {
                "t_low_fold": t_low_fold,
                "t_high_fold": t_high_fold,
                "gray_zone_inverted_order": bool(t_low_fold > t_high_fold),
                # The binary threshold and the two gray thresholds are picked by
                # independent criteria, so nothing forces thr_fold to lie inside
                # the gray band. When it does not, a sample can be called
                # "responder" by the binary rule and "non_responder_confident" by
                # the gray rule. Recorded per fold so the contradiction is
                # visible instead of silent.
                "threshold_inside_gray_zone": bool(t_low_fold <= thr_fold <= t_high_fold),
                "subsample_thresholds_low": [float(t) for t in subsample_thresholds_low],
                "subsample_thresholds_high": [float(t) for t in subsample_thresholds_high],
            }

        fold_details.append({
            "fold": fold_id,
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "n_pos_test": int(np.sum(y_test == 1)),
            "n_neg_test": int(np.sum(y_test == 0)),
            "threshold_fold": thr_fold,
            "subsample_thresholds_single": [float(t) for t in subsample_thresholds_single],
            "threshold_fit_mode": threshold_fit_mode,
            **gray_info
        })

    # Basic check
    if not fold_covered.all() or np.any(oof_pred < 0):
        raise RuntimeError("OOF predictions are incomplete.")

    if use_gray_zone and np.any(oof_pred_gray == -999):
        raise RuntimeError("OOF gray-zone predictions are incomplete.")

    # -------------------------------------------------------------------------
    # 3) ROC/AUC (RESUBSTITUTION) + out-of-fold confusion matrix and metrics
    #
    # The ROC below is NOT out-of-fold, despite what this block used to be
    # called. score_oriented is just the oriented score for every sample --
    # cross-validation never touched it, because the score is computed upstream
    # without using the labels. So this is the plain resubstitution ROC over the
    # whole cohort and the AUC carries no cross-validation protection. It is
    # still a fair summary of the score's separation IF the signature was fixed
    # in advance; it is optimistic exactly when it matters, i.e. if the signature
    # was chosen by looking at these same patients.
    #
    # The confusion matrix and everything derived from it ARE cross-validated:
    # they come from oof_pred, where each sample was classified by a threshold
    # fitted on folds that excluded it.
    # -------------------------------------------------------------------------
    fpr, tpr, roc_thresholds = roc_curve(y, score_oriented)
    roc_auc = float(auc(fpr, tpr))

    standard_metrics, cm = _compute_binary_metrics(y, oof_pred)
    cm_df = pd.DataFrame(
        cm,
        index=["True_NonResponder(0)", "True_Responder(1)"],
        columns=["Pred_NonResponder(0)", "Pred_Responder(1)"]
    )

    # -------------------------------------------------------------------------
    # Bootstrap CI for standard metrics (+ AUC)
    #
    # "conditional" (default): resample (y_true, y_pred) with the out-of-fold
    #   predictions FROZEN. Fast (~3 s at 1000 replicates) but the interval is
    #   conditional on the threshold that was fitted -- it excludes the error in
    #   estimating that threshold.
    #
    # "pipeline": resample the cohort and re-run the ENTIRE procedure on each
    #   resample (fold partition, subsample bagging, threshold selection, gray
    #   zone), so threshold-estimation and partition variance are included.
    #   ~160 ms per replicate, i.e. ~4 min at 1000 -- opt in for a final run.
    #
    # Measured difference on a 17R/23NR cohort (pipeline / conditional widths):
    #   sensitivity 1.00  specificity 0.87  ppv 0.80  npv 1.12  auc 1.08
    # so the conditional interval is at worst ~20% too narrow and sometimes
    # wider. Do not expect the slow mode to change conclusions; it exists so the
    # number quoted in a paper does not rest on a conditional interval.
    #
    # CAVEAT on "pipeline": a bootstrap resample contains duplicates, and a
    # duplicate can land in both a training and a test fold, which leaks and
    # makes this interval somewhat too NARROW as well. It is closer to honest
    # than the conditional one, not exactly honest. A leak-free version would
    # need out-of-bag evaluation (.632+ family), which changes the estimator.
    # -------------------------------------------------------------------------
    if ci_method not in ("conditional", "pipeline"):
        raise ValueError("ci_method must be 'conditional' or 'pipeline'.")

    _PIPE_METRICS = ["sensitivity", "specificity", "ppv", "npv", "fpr", "fnr",
                     "accuracy", "balanced_accuracy",
                     "risk_false_positive_given_pred_positive",
                     "risk_false_negative_given_pred_negative", "auc"]
    pipeline_rows = None

    n_pipeline = min(int(n_bootstrap), _PIPELINE_MAX_BOOTSTRAP)

    if ci_method == "pipeline" and n_pipeline > 0:
        _idx_pos_p = np.where(y == 1)[0]
        _idx_neg_p = np.where(y == 0)[0]
        _idx_all_p = np.arange(len(y))

        def _one_pipeline_boot(child_seed):
            rng_p = np.random.default_rng(child_seed)
            if bootstrap_stratified and len(_idx_pos_p) > 0 and len(_idx_neg_p) > 0:
                bi = np.concatenate([
                    rng_p.choice(_idx_pos_p, size=len(_idx_pos_p), replace=True),
                    rng_p.choice(_idx_neg_p, size=len(_idx_neg_p), replace=True),
                ])
            else:
                bi = rng_p.choice(_idx_all_p, size=len(_idx_all_p), replace=True)
            if len(np.unique(y[bi])) < 2:
                return None
            # Duplicated draws need distinct IDs: the CV requires a unique index
            # and builds its label vector from set()s, which would silently
            # collapse the resample back to fewer samples.
            ids_b = [f"b{k:05d}" for k in range(len(bi))]
            df_b = pd.DataFrame({score_col: score_oriented[bi]}, index=ids_b)
            resp_b = [ids_b[k] for k in range(len(bi)) if y[bi[k]] == 1]
            non_b = [ids_b[k] for k in range(len(bi)) if y[bi[k]] == 0]
            try:
                out = evaluate_signature_threshold_cv(
                    df_b, resp_b, non_b,
                    score_col=score_col,
                    n_folds=n_folds,
                    n_subsamples=n_subsamples,
                    shuffle=shuffle,
                    # vary the partition too, so fold-partition variance is
                    # sampled rather than held fixed across replicates
                    random_state=int(rng_p.integers(0, 2**31 - 1)),
                    threshold_criterion=threshold_criterion,
                    higher_score_is_responder=True,   # score_oriented is already oriented
                    use_gray_zone=use_gray_zone,
                    gray_target_sensitivity=gray_target_sensitivity,
                    gray_target_specificity=gray_target_specificity,
                    n_bootstrap=0,                    # no nested bootstrap
                    ci_method="conditional",          # and no recursion
                    return_oof_predictions=False,
                )
            except Exception:
                return None
            row = {k: out["metrics"].get(k, np.nan) for k in _PIPE_METRICS
                   if k != "auc"}
            row["auc"] = out["auc"]
            if use_gray_zone and out.get("gray_zone") is not None:
                row["coverage"] = out["gray_zone"]["coverage"]
                row["rejection_rate"] = out["gray_zone"]["rejection_rate"]
            return row

        _pipe_seeds = np.random.SeedSequence(bootstrap_seed + 7).spawn(n_pipeline)
        _njp = _n_jobs()
        if _njp == 1:
            _pipe_out = [_one_pipeline_boot(sd) for sd in _pipe_seeds]
        else:
            try:
                _pipe_out = Parallel(n_jobs=_njp, prefer="threads")(
                    delayed(_one_pipeline_boot)(sd) for sd in _pipe_seeds
                )
            except Exception:
                _pipe_out = [_one_pipeline_boot(sd) for sd in _pipe_seeds]
        pipeline_rows = [r_ for r_ in _pipe_out if r_ is not None]

    if pipeline_rows:
        _pipe_df = pd.DataFrame(pipeline_rows)
        standard_ci = {}
        for _m in _pipe_df.columns:
            _lo, _hi = _percentile_ci(_pipe_df[_m].to_numpy(), ci_level=ci_level)
            standard_ci[_m] = {"ci_low": _lo, "ci_high": _hi}
    else:
        standard_ci = _bootstrap_ci_binary_metrics(
            y_true=y,
            y_pred=oof_pred,
            y_score=score_oriented,
            n_boot=int(n_bootstrap),
            ci_level=ci_level,
            seed=bootstrap_seed,
            stratified=bootstrap_stratified
        )
    standard_ci["auc"] = standard_ci.get("auc", {"ci_low": np.nan, "ci_high": np.nan})

    standard_metric_estimates = {**standard_metrics, "auc": roc_auc}
    standard_metrics_table = _metrics_ci_table(
        standard_metric_estimates,
        standard_ci,
        metric_order=[
            "auc", "sensitivity", "specificity", "ppv", "npv",
            "fpr", "fnr", "accuracy", "balanced_accuracy",
            "risk_false_positive_given_pred_positive",
            "risk_false_negative_given_pred_negative",
            "prevalence"
        ],
        population=f"full cohort (n={len(y)})"
    )

    # prevalence is deterministic here (bootstrap CI omitted unless you want to add)
    # fill prevalence CI with NaN (already)
    # -------------------------------------------------------------------------
    # 4) Gray zone: coverage, realised rates, confident-subset metrics
    # -------------------------------------------------------------------------
    gray_zone_results = None
    if use_gray_zone:
        # 3-state confusion-like table: true class x predicted class {0, uncertain, 1}
        pred_labels_gray = oof_pred_gray.copy()
        # map uncertain -1 -> column label
        three_way_table = pd.crosstab(
            pd.Series(y, name="true"),
            pd.Series(pred_labels_gray, name="pred_gray"),
            dropna=False
        ).reindex(index=[0, 1], columns=[0, -1, 1], fill_value=0)

        three_way_table.index = ["True_NonResponder(0)", "True_Responder(1)"]
        three_way_table.columns = ["Pred_0_confident", "Pred_uncertain", "Pred_1_confident"]

        # ---------------------------------------------------------------
        # What the gray zone ACTUALLY delivered, out-of-fold.
        #
        # t_low is picked as the highest threshold whose TRAINING-ROC tpr still
        # reaches the target -- the boundary of the constraint, so the maximum-
        # optimism choice. The realised rate is therefore systematically short.
        # Measured over 30 cohorts per configuration (target 0.90):
        #     17R/23NR  realised sensitivity 0.882  (short by 0.018, 100% of runs)
        #     12R/28NR  realised sensitivity 0.833  (short by 0.067, 100% of runs)
        #     25R/25NR  realised sensitivity 0.879  (short by 0.021,  90% of runs)
        # Specificity is fine (+0.001 to +0.027): that axis has more samples.
        #
        # The selection rule is deliberately NOT corrected, because every remedy
        # tested is worse than the bias:
        #   - requiring a lower confidence bound on tpr to clear the target needs
        #     >= 40 responders (at 17, even a perfect 17/17 bounds at 0.838), so
        #     it would reject every threshold and always hit the fallback;
        #   - inflating the nominal target by the measured shortfall fixed
        #     12R/28NR (0.840 -> 0.917) and 25R/25NR (0.880 -> 0.922) but did
        #     nothing for 17R/23NR (0.887 -> 0.882), and cost 6-11 points of
        #     coverage where it did work.
        # The root cause is not fixable in code: realised sensitivity is
        # quantised in steps of 1/n_responders, so with 17 responders the
        # attainable values around the target are 0.8824 and 0.9412 and 0.90
        # simply does not exist. A "90% guarantee" is not well-posed here.
        #
        # So report the realised rate instead of promising the nominal one.
        # ---------------------------------------------------------------
        _resp_mask = (y == 1)
        _nonresp_mask = (y == 0)
        realised_sens_below_t_low = (
            float(1 - np.mean(oof_pred_gray[_resp_mask] == 0)) if _resp_mask.any() else np.nan
        )
        realised_spec_above_t_high = (
            float(1 - np.mean(oof_pred_gray[_nonresp_mask] == 1)) if _nonresp_mask.any() else np.nan
        )

        confident_mask = (oof_pred_gray != -1)
        n_confident = int(np.sum(confident_mask))
        coverage = _safe_div(n_confident, len(y))
        rejection_rate = np.nan if np.isnan(coverage) else 1 - coverage

        # ---------------------------------------------------------------
        # Bootstrap CI for coverage / rejection_rate.
        #
        # These used to be reported as bare point estimates with a NaN CI,
        # although they are among the most uncertain numbers the function
        # produces: coverage depends entirely on t_low/t_high, which are
        # themselves estimated. Against a full re-run-the-pipeline bootstrap
        # (the expensive gold standard) coverage on a 17R/23NR cohort came out
        # at [0.55, 1.00] around a point estimate of 0.85 -- an interval wide
        # enough to change how the gray zone should be read.
        #
        # Cheap middle path: refit ONLY the gray thresholds on each resample,
        # then apply them to the whole cohort. Measured against that gold
        # standard on a 17R/23NR cohort (gold [0.549, 1.000], width 0.451):
        #     refit + evaluate in-sample     [0.450, 1.000]  width 0.550
        #     refit + evaluate full cohort   [0.475, 1.000]  width 0.525  <- used
        #     refit + evaluate out-of-bag    [0.429, 1.000]  width 0.571
        #     bagged refit, full cohort      [0.475, 1.000]  width 0.525  (6.4 s)
        # Out-of-bag is worse, not better: ~15 samples per replicate inject
        # binomial noise that is not part of the uncertainty in a 40-sample
        # statistic. Bagging the refit to mirror the real procedure costs 30x
        # more and lands on the same interval, so it is not worth it.
        #
        # All cheap variants stay WIDER than the gold standard (~16% here). The
        # real procedure bags thresholds over subsamples and then takes a medoid
        # across folds, which stabilises them more than any single refit can;
        # this bootstrap does not reproduce that averaging. The interval is
        # therefore conservative -- the safe direction for a clinical readout,
        # but it should not be quoted as exact.
        # ---------------------------------------------------------------
        # In "pipeline" mode the whole procedure was already re-run per
        # replicate, which subsumes this refit and is strictly better -- reuse it
        # rather than paying for a second, weaker bootstrap.
        _cov_from_pipeline = bool(pipeline_rows) and ("coverage" in standard_ci)

        _idx_all_cov = np.arange(len(y))
        _idx_pos_cov = np.where(y == 1)[0]
        _idx_neg_cov = np.where(y == 0)[0]

        def _one_coverage_boot(child_seed):
            rng_c = np.random.default_rng(child_seed)
            if bootstrap_stratified and len(_idx_pos_cov) > 0 and len(_idx_neg_cov) > 0:
                bi = np.concatenate([
                    rng_c.choice(_idx_pos_cov, size=len(_idx_pos_cov), replace=True),
                    rng_c.choice(_idx_neg_cov, size=len(_idx_neg_cov), replace=True),
                ])
            else:
                bi = rng_c.choice(_idx_all_cov, size=len(_idx_all_cov), replace=True)
            if len(np.unique(y[bi])) < 2:
                return np.nan
            try:
                _tl, _th, _ = _pick_gray_thresholds_from_roc(
                    y[bi], score_oriented[bi],
                    target_sens=gray_target_sensitivity,
                    target_spec=gray_target_specificity,
                )
            except Exception:
                return np.nan
            return float(np.mean(_gray_zone_labels(score_oriented, _tl, _th) != -1))

        if _cov_from_pipeline:
            coverage_ci = standard_ci["coverage"]
            rejection_ci = standard_ci["rejection_rate"]
            n_cov_valid = len(pipeline_rows)
            coverage_ci_source = "pipeline (whole procedure re-run per replicate)"
        else:
            _cov_seeds = np.random.SeedSequence(bootstrap_seed + 2).spawn(int(n_bootstrap))
            _nj = _n_jobs()
            if _nj == 1:
                _cov_vals = [_one_coverage_boot(s) for s in _cov_seeds]
            else:
                try:
                    _cov_vals = Parallel(n_jobs=_nj, prefer="threads")(
                        delayed(_one_coverage_boot)(s) for s in _cov_seeds
                    )
                except Exception:
                    _cov_vals = [_one_coverage_boot(s) for s in _cov_seeds]

            _cov_lo, _cov_hi = _percentile_ci(_cov_vals, ci_level=ci_level)
            coverage_ci = {"ci_low": _cov_lo, "ci_high": _cov_hi}
            # rejection = 1 - coverage, so the bounds swap
            rejection_ci = {
                "ci_low": np.nan if np.isnan(_cov_hi) else 1 - _cov_hi,
                "ci_high": np.nan if np.isnan(_cov_lo) else 1 - _cov_lo,
            }
            n_cov_valid = int(np.sum(np.isfinite(np.asarray(_cov_vals, dtype=float))))
            coverage_ci_source = "gray thresholds refit per replicate"

        if n_confident > 0:
            y_conf = y[confident_mask]
            p_conf = oof_pred_gray[confident_mask].astype(int)  # now only 0/1
            s_conf = score_oriented[confident_mask]

            gray_metrics_conf, cm_gray_conf = _compute_binary_metrics(y_conf, p_conf)
            cm_gray_conf_df = pd.DataFrame(
                cm_gray_conf,
                index=["True_NonResponder(0)", "True_Responder(1)"],
                columns=["Pred_NonResponder(0)", "Pred_Responder(1)"]
            )

            gray_ci = _bootstrap_ci_binary_metrics(
                y_true=y_conf,
                y_pred=p_conf,
                y_score=s_conf if len(np.unique(y_conf)) > 1 else None,
                n_boot=int(n_bootstrap),
                ci_level=ci_level,
                seed=bootstrap_seed + 1,
                stratified=bootstrap_stratified
            )

            # Optional AUC on confident subset
            if "auc" not in gray_ci:
                gray_ci["auc"] = {"ci_low": np.nan, "ci_high": np.nan}
            if len(np.unique(y_conf)) > 1:
                fpr_conf, tpr_conf, _ = roc_curve(y_conf, s_conf)
                auc_conf = float(auc(fpr_conf, tpr_conf))
            else:
                auc_conf = np.nan

            gray_ci["coverage"] = coverage_ci
            gray_ci["rejection_rate"] = rejection_ci

            gray_estimates = {
                **gray_metrics_conf,
                "coverage": coverage,
                "rejection_rate": rejection_rate,
                "realised_sensitivity_below_t_low": realised_sens_below_t_low,
                "realised_specificity_above_t_high": realised_spec_above_t_high,
                "auc": auc_conf
            }

            gray_metrics_table = _metrics_ci_table(
                gray_estimates,
                gray_ci,
                metric_order=[
                    "coverage", "rejection_rate",
                    "realised_sensitivity_below_t_low",
                    "realised_specificity_above_t_high",
                    "auc", "sensitivity", "specificity", "ppv", "npv",
                    "fpr", "fnr", "accuracy", "balanced_accuracy",
                    "risk_false_positive_given_pred_positive",
                    "risk_false_negative_given_pred_negative"
                ],
                population=(
                    f"CONFIDENT SUBSET ONLY (n={n_confident} of {len(y)}, "
                    f"coverage {coverage:.0%}) - self-selected, NOT comparable "
                    f"to the full-cohort table"
                )
            )
        else:
            cm_gray_conf = np.array([[0, 0], [0, 0]])
            cm_gray_conf_df = pd.DataFrame(
                cm_gray_conf,
                index=["True_NonResponder(0)", "True_Responder(1)"],
                columns=["Pred_NonResponder(0)", "Pred_Responder(1)"]
            )
            gray_metrics_conf = {
                "n": 0,
                "coverage": 0.0,
                "rejection_rate": 1.0
            }
            gray_ci = {"coverage": coverage_ci, "rejection_rate": rejection_ci}
            gray_metrics_table = _metrics_ci_table(
                {"coverage": coverage, "rejection_rate": rejection_rate},
                gray_ci,
                metric_order=["coverage", "rejection_rate"],
            )

        # Summary of gray thresholds across folds
        fold_gray_df = pd.DataFrame(fold_gray_thresholds, columns=["t_low", "t_high"]) \
            if len(fold_gray_thresholds) > 0 else pd.DataFrame(columns=["t_low", "t_high"])

        # t_low and t_high are aggregated by INDEPENDENT medians, so the summary
        # pair need not correspond to any single fold's pair and can invert even
        # when no individual fold inverted. Report both the per-fold degeneracies
        # and the aggregate one rather than assuming a well-ordered band.
        median_t_low = float(np.median(fold_gray_df["t_low"])) if len(fold_gray_df) else np.nan
        median_t_high = float(np.median(fold_gray_df["t_high"])) if len(fold_gray_df) else np.nan
        thr_median = float(np.median(fold_thresholds)) if len(fold_thresholds) else np.nan

        # The pair actually reported downstream. Taking the two medians
        # INDEPENDENTLY can yield a (t_low, t_high) that no fold ever produced
        # -- e.g. folds [(1,5), (2,4), (3,9)] give medians (2,5), a band that
        # never existed. Instead take the medoid: the fold whose pair is closest
        # (L1) to the componentwise median, which is by construction a band a
        # real fold produced.
        #
        # NB the original audit also claimed independent medians "can invert
        # even when no fold inverted". That is false: if low_i <= high_i for
        # every fold then the k-th order statistics satisfy low_(k) <= high_(k),
        # so median(low) <= median(high). Verified over 200k random fold sets.
        # Inversion of the final pair therefore only happens when a real fold
        # inverted, which n_folds_gray_inverted already reports.
        if len(fold_gray_df):
            _d = ((fold_gray_df["t_low"] - median_t_low).abs()
                  + (fold_gray_df["t_high"] - median_t_high).abs()).to_numpy()
            _k = int(np.nanargmin(_d))
            t_low_final_cv = float(fold_gray_df["t_low"].iloc[_k])
            t_high_final_cv = float(fold_gray_df["t_high"].iloc[_k])
            gray_pair_source_fold = int(_k + 1)
        else:
            t_low_final_cv = np.nan
            t_high_final_cv = np.nan
            gray_pair_source_fold = None

        gray_zone_results = {
            "targets": {
                "target_sensitivity_for_t_low": gray_target_sensitivity,
                "target_specificity_for_t_high": gray_target_specificity,
            },
            # Requested vs delivered. The targets are in-sample constraints on
            # the training ROC; these are what the out-of-fold labels achieved.
            "realised": {
                "sensitivity_below_t_low": realised_sens_below_t_low,
                "specificity_above_t_high": realised_spec_above_t_high,
                "sensitivity_shortfall": (
                    np.nan if np.isnan(realised_sens_below_t_low)
                    else realised_sens_below_t_low - gray_target_sensitivity
                ),
                "specificity_shortfall": (
                    np.nan if np.isnan(realised_spec_above_t_high)
                    else realised_spec_above_t_high - gray_target_specificity
                ),
                # Smallest achievable step in each rate. If the target is not a
                # multiple of this, it cannot be hit exactly no matter what the
                # selection rule does.
                "sensitivity_granularity": (1.0 / n_pos) if n_pos else np.nan,
                "specificity_granularity": (1.0 / n_neg) if n_neg else np.nan,
                "target_sensitivity_attainable": bool(
                    n_pos and abs(round(gray_target_sensitivity * n_pos)
                                  - gray_target_sensitivity * n_pos) < 1e-9
                ),
                "note": (
                    "Targets constrain the TRAINING ROC; the selection takes the "
                    "boundary threshold, so the realised out-of-fold rate runs "
                    "short by ~0.02-0.07 at these cohort sizes. Quote the "
                    "realised values, not the targets."
                ),
            },
            "fold_gray_thresholds": fold_gray_thresholds,
            "fold_gray_thresholds_df": fold_gray_df,
            # The coherent pair to use (medoid fold). median_* are kept as
            # diagnostics only -- see the note above on independent aggregation.
            "t_low_final": t_low_final_cv,
            "t_high_final": t_high_final_cv,
            "gray_pair_source_fold": gray_pair_source_fold,
            "median_t_low_across_folds": median_t_low,
            "median_t_high_across_folds": median_t_high,

            # Degeneracy flags. inverted => t_low > t_high, so the confident
            # regions overlap and _gray_zone_labels marks the overlap uncertain.
            # threshold_inside => the binary threshold sits within the band; when
            # False, some samples get contradictory binary vs gray labels.
            "gray_zone_inverted_final": bool(t_low_final_cv > t_high_final_cv),
            "gray_zone_inverted_median": bool(median_t_low > median_t_high),
            "n_folds_gray_inverted": int(sum(1 for fd in fold_details if fd.get("gray_zone_inverted_order"))),
            "threshold_inside_gray_zone_final": bool(t_low_final_cv <= thr_median <= t_high_final_cv),
            "threshold_inside_gray_zone_median": bool(median_t_low <= thr_median <= median_t_high),
            "n_folds_threshold_outside_gray_zone": int(
                sum(1 for fd in fold_details if fd.get("threshold_inside_gray_zone") is False)
            ),

            "three_way_table": three_way_table,
            "confident_subset_confusion_matrix": cm_gray_conf,
            "confident_subset_confusion_matrix_df": cm_gray_conf_df,
            "coverage_ci_source": coverage_ci_source,
            "coverage_ci_is_conservative": not _cov_from_pipeline,
            "coverage_ci_method": (
                "gray thresholds refit on each stratified resample and applied "
                "to the whole cohort; ~16% wider than a full re-run-the-pipeline "
                "bootstrap because it does not reproduce the fold-level "
                "averaging that stabilises the real thresholds"
            ),
            "confident_subset_n": n_confident,
            "confident_subset_of_n": int(len(y)),
            "confident_subset_is_self_selected": True,
            "confident_subset_note": (
                "Metrics below describe only the samples the gray zone called "
                "confident. Membership was chosen using these same data, so they "
                "are conditional on that selection and are systematically better "
                "than the full-cohort metrics. Read them with coverage, never as "
                "an improved version of the full-cohort numbers."
            ),
            "confident_subset_metrics": gray_metrics_conf,
            "confident_subset_bootstrap_ci": gray_ci,
            "confident_subset_metrics_table": gray_metrics_table,
            "coverage": coverage,
            "rejection_rate": rejection_rate,
        }

    # -------------------------------------------------------------------------
    # 5) Per-sample out-of-fold table
    # -------------------------------------------------------------------------
    oof_df = None
    if return_oof_predictions:
        oof_df = pd.DataFrame({
            "sample_id": sample_ids,
            "y_true": y,
            "score_raw_original": df.loc[sample_ids, "score_raw"].values,
            "score_used_oriented": score_oriented,
            "pred_binary_oof": oof_pred,
            "threshold_binary_fold": oof_thr,
        }).set_index("sample_id")

        if use_gray_zone:
            oof_df["pred_gray_oof"] = oof_pred_gray  # -1 uncertain, 0, 1
            oof_df["t_low_fold"] = oof_t_low
            oof_df["t_high_fold"] = oof_t_high
            oof_df["is_confident_gray"] = (oof_df["pred_gray_oof"] != -1)

    # -------------------------------------------------------------------------
    # 6) Assemble the result
    # -------------------------------------------------------------------------
    results = {
        "n_samples": int(len(y)),
        "n_positive": int(np.sum(y == 1)),
        "n_negative": int(np.sum(y == 0)),
        "score_orientation": orientation_txt,

        # ROC / AUC -- RESUBSTITUTION over the whole cohort, not out-of-fold.
        "fpr": fpr,
        "tpr": tpr,
        "roc_thresholds": roc_thresholds,
        "auc": roc_auc,
        "auc_is_resubstitution": True,
        "auc_interpretation": (
            "Resubstitution ROC/AUC over the whole cohort. The score is computed "
            "upstream without the labels, so no cross-validation is applied to it "
            "and none is needed IF the signature was fixed in advance. Optimistic "
            "if the signature was selected by looking at these same patients. "
            "Only the threshold-dependent metrics (Se/Sp/PPV/NPV/accuracy) are "
            "cross-validated."
        ),

        # Binary classification from OUT-OF-FOLD predictions
        "confusion_matrix": cm,
        "confusion_matrix_df": cm_df,
        "metrics": standard_metrics,
        "bootstrap_ci": standard_ci,
        # What bootstrap_ci actually brackets. It resamples (y_true, y_pred)
        # with y_pred FROZEN, so it is the sampling variability of the metric
        # given the threshold that was fitted -- it does not include the error
        # in estimating that threshold. Measured against a full
        # re-run-the-pipeline bootstrap the difference is modest (widths within
        # ~0-20%, worst case PPV at 0.80x) because the fold thresholds barely
        # move, but the interval is conditional and should be described as such.
        # coverage / rejection_rate in gray_zone are the exception: those DO
        # refit the thresholds per replicate.
        "bootstrap_ci_interpretation": (
            (
                "Pipeline bootstrap: the whole procedure (fold partition, "
                "subsample bagging, threshold selection, gray zone) was re-run on "
                "each cohort resample, so threshold-estimation and partition "
                "variance are included. Still slightly narrow: duplicated draws "
                "can appear in both a training and a test fold."
            )
            if pipeline_rows else
            (
                "Conditional on the fitted threshold: resamples the cohort with "
                "the out-of-fold predictions held fixed, so threshold-estimation "
                "error is excluded. gray_zone coverage/rejection_rate CIs do "
                "refit. Pass ci_method='pipeline' for the slower honest version."
            )
        ),
        "n_pipeline_bootstrap_valid": (len(pipeline_rows) if pipeline_rows else 0),
        "n_pipeline_bootstrap_requested": (n_pipeline if ci_method == "pipeline" else 0),
        "pipeline_bootstrap_capped": bool(
            ci_method == "pipeline" and int(n_bootstrap) > _PIPELINE_MAX_BOOTSTRAP
        ),
        "metrics_table": standard_metrics_table,

        # Thresholds / folds
        "fold_thresholds": fold_thresholds,
        "fold_details": fold_details,

        # Run configuration. Persisted so downstream consumers
        # (sample_confidence_report, py_plots) re-use the rule this run actually
        # applied instead of falling back to their own default: running the CV
        # with "closest_topleft" and reporting with "youden" used to be silent
        # and unrecoverable from the results dict.
        "config": {
            "ci_method": ci_method,
            "threshold_criterion": threshold_criterion,
            "threshold_aggregation": "median_over_subsamples",
            "n_folds": int(n_folds),
            "n_folds_effective": int(n_folds_eff),
            "n_subsamples": int(n_subsamples),
            "shuffle": bool(shuffle),
            "random_state": random_state,
            "higher_score_is_responder": bool(higher_score_is_responder),
            "use_gray_zone": bool(use_gray_zone),
            "gray_target_sensitivity": gray_target_sensitivity,
            "gray_target_specificity": gray_target_specificity,
            "n_bootstrap": int(n_bootstrap),
            "ci_level": ci_level,
            "bootstrap_seed": bootstrap_seed,
            "bootstrap_stratified": bool(bootstrap_stratified),
        },
    }

    if use_gray_zone:
        results["gray_zone"] = gray_zone_results

    if return_oof_predictions:
        results["oof"] = oof_df

    return results


def sample_confidence_report(
    scores,
    cv_result: dict,
    sample_ids=None,
    threshold_criterion=None,  # None => read from res["config"] (the rule the CV actually used)
    calibration_C: float = 1.0,  # L2 strength for the calibration, on the STANDARDISED score
    use_gray_zone: bool = True,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    bootstrap_seed: int = 123,
    bootstrap_stratified: bool = True,
):
    """
    Build a per-sample confidence report using outputs from evaluate_signature_threshold_cv.

    Parameters
    ----------
    scores : scalar | list | np.ndarray | pd.Series
        Raw scores (same scale as df_scores[score_col] used in v2).
    cv_result : dict
        Output of evaluate_signature_threshold_cv(...), ideally with cv_result["oof"] present.
    sample_ids : list-like, optional
        Optional IDs for the query scores. If None, auto-generated.
    threshold_criterion : str
        "youden" or "closest_topleft" for bootstrap threshold re-estimation.
    use_gray_zone : bool
        If True and gray_zone exists in cv_result, returns gray-zone predictions.
    n_bootstrap : int
        Bootstrap replicates for CI.
    ci_level : float
        e.g. 0.95
    bootstrap_seed : int
        RNG seed.
    bootstrap_stratified : bool
        Stratified bootstrap on the OOF reference set.

    Returns
    -------
    out : dict
        {
          "report": DataFrame (1 row per query sample),
          "bootstrap_summary": dict,
          "reference_info": dict
        }

    Notes
    -----
    - Probabilité calibrée: logistic regression (logit) fitted on OOF scores/labels from cv_result["oof"].
    - Threshold used for binary pred: median(cv_result["fold_thresholds"]).
    - confidence / ambiguity / margin_prob are measured from the DECISION
      BOUNDARY p_at_threshold = P(responder | score = threshold), not from 0.5.
      Youden's threshold only coincides with p = 0.5 on a balanced cohort; on an
      imbalanced one it drifts far from it (~0.20 at 8R/32NR), and measuring
      "confidence" from 0.5 there contradicts the prediction itself.
    - p, risk_FP_if_pred_positive, risk_FN_if_pred_negative and
      conditional_error_risk are left on the true probability scale and are NOT
      recentred: they answer "how likely is this call wrong", which only means
      something as a real probability. conditional_error_risk above 0.5 is
      therefore possible and legitimate; pred_against_majority_prob marks it.
    - Threshold stability: proportion of BOOTSTRAP thresholds below the sample
      score. The CV-based equivalent (proportion of fold thresholds below the
      score) was removed: fold-training sets overlap ~94%, so the fold
      thresholds are near-duplicates and the statistic collapsed to 0/1 for
      ~90% of scores while looking like independent corroboration.
    - Bootstrap CIs are computed from refit-on-bootstrap replicates of the calibration + threshold.
    """

    # -------------------------------------------------------------------------
    # Helpers
    # -------------------------------------------------------------------------
    def _to_1d_array(x):
        if np.isscalar(x):
            return np.array([x], dtype=float)
        if isinstance(x, pd.Series):
            return x.astype(float).to_numpy()
        return np.asarray(x, dtype=float).reshape(-1)

    def _fit_calibration(scores_ref, y_ref):
        """
        Logistic calibration p(y=1|score), standardised and regularised.

        StandardScaler is NOT cosmetic. The L2 penalty acts on the coefficient,
        so without it the amount of shrinkage depends on the arbitrary scale of
        whatever score compute_combined_score() produced: the same data at score
        scale 0.01 / 1 / 100 fits coefficients 0.073 / 0.540 / 0.006 under the
        same C. Standardising first makes C mean the same thing for every
        signature. It changes no probability at a given effective penalty --
        it only makes the penalty comparable.

        C was 1e6 ("quasi non-regularized"), which is fine while the classes
        overlap but collapses once they separate. Measured by stratified-CV log
        loss on the reference set (lower is better):

                              C=1e6    C=10     C=1     C=0.1
          balanced/overlap    0.6693  0.6689  0.6669  0.6717
          PREDIMEL-like       0.6538  0.6530  0.6489  0.6566
          imbalanced 8/32     0.5247  0.5241  0.5198  0.5101
          well separated      0.3768  0.3647  0.3766  0.5131
          near-separable      1.2835  0.1615  0.2317  0.4682   <-- 8x worse

        In that last regime C=1e6 drove the coefficient to 7.13 and saturated
        67% of samples to |p - 0.5| > 0.49, i.e. it reported near-certainty for
        two thirds of the cohort. C=1 lands within 0.01 log loss of the best
        grid value on 4 of 5 cohorts and is never catastrophic, so it is the
        default; calibration_C exposes it. (lbfgs converges in 4-11 iterations
        in every case tested, so max_iter=1000 was never the binding issue.)
        """
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(
                solver="lbfgs",
                C=float(calibration_C),
                max_iter=1000
            )
        )
        model.fit(scores_ref.reshape(-1, 1), y_ref.astype(int))
        return model

    def _gray_label(pred_gray):
        if pred_gray == 0:
            return "non_responder_confident"
        if pred_gray == 1:
            return "responder_confident"
        return "uncertain"

    def _signed_gray_distance(s, t_low, t_high):
        """
        SIGNED distance to the nearest gray-zone boundary:
            > 0  outside the zone (confident) -- how far past the boundary
            < 0  inside the zone (uncertain)  -- how far from the nearest edge
            = 0  exactly ON a boundary: ambiguous by nature, since s == t_high
                 is confident and s == t_low is uncertain, both at distance 0.
                 pred_gray_label is authoritative there, not the sign.

        Unsigned, 0.2 meant both "0.2 past the boundary, confident" and "0.2
        inside the zone, uncertain", so any sort or plot on this column mixed
        the two populations. The sign is taken from _gray_zone_labels rather
        than recomputed, so it cannot disagree with pred_gray_label -- including
        on an inverted zone, where the naive branch order gets it backwards.
        """
        d = float(min(abs(s - t_low), abs(s - t_high)))
        return d if int(_gray_zone_labels(s, t_low, t_high)) != -1 else -d

    # -------------------------------------------------------------------------
    # Validate cv_result
    # -------------------------------------------------------------------------
    if "oof" not in cv_result:
        raise ValueError(
            "cv_result must contain 'oof'. Please run evaluate_signature_threshold_cv(..., return_oof_predictions=True)."
        )

    oof = cv_result["oof"].copy()
    required_oof_cols = {"y_true", "score_raw_original", "score_used_oriented"}
    missing = required_oof_cols - set(oof.columns)
    if missing:
        raise ValueError(f"cv_result['oof'] is missing columns: {sorted(missing)}")

    if "fold_thresholds" not in cv_result or len(cv_result["fold_thresholds"]) == 0:
        raise ValueError("cv_result must contain non-empty 'fold_thresholds'.")

    # Threshold rule: default to whatever the CV run actually used, so the
    # bootstrap thresholds re-estimated below match the fold thresholds this
    # report is comparing them against. An explicit argument still wins.
    cv_config = cv_result.get("config", {}) or {}
    if threshold_criterion is None:
        threshold_criterion = cv_config.get("threshold_criterion", "youden")
    if threshold_criterion not in ("youden", "closest_topleft"):
        raise ValueError("threshold_criterion must be 'youden' or 'closest_topleft'.")

    # -------------------------------------------------------------------------
    # Prepare query scores
    # -------------------------------------------------------------------------
    scores_raw_query = _to_1d_array(scores)
    n_query = len(scores_raw_query)

    if sample_ids is None:
        sample_ids = [f"query_{i+1}" for i in range(n_query)]
    if len(sample_ids) != n_query:
        raise ValueError("sample_ids length must match number of query scores.")

    # Infer orientation (raw -> oriented) using OOF columns
    # If score_used_oriented == score_raw_original => sign=+1, else -1
    raw_ref = oof["score_raw_original"].to_numpy(dtype=float)
    ori_ref = oof["score_used_oriented"].to_numpy(dtype=float)

    # robust sign inference
    diff_plus = np.nanmean(np.abs(ori_ref - raw_ref))
    diff_minus = np.nanmean(np.abs(ori_ref + raw_ref))
    sign = 1.0 if diff_plus <= diff_minus else -1.0

    scores_query_oriented = sign * scores_raw_query

    # Reference set for calibration/bootstraps
    y_ref = oof["y_true"].to_numpy(dtype=int)
    s_ref = oof["score_used_oriented"].to_numpy(dtype=float)

    # -------------------------------------------------------------------------
    # Final thresholds used for direct prediction (from v2 summary)
    # -------------------------------------------------------------------------
    fold_thresholds = np.asarray(cv_result["fold_thresholds"], dtype=float)
    thr_final = float(np.median(fold_thresholds))

    has_gray = use_gray_zone and ("gray_zone" in cv_result) and (cv_result["gray_zone"] is not None)
    if has_gray:
        gray = cv_result["gray_zone"]
        # Try robustly to fetch fold gray thresholds
        if "fold_gray_thresholds" in gray and len(gray["fold_gray_thresholds"]) > 0:
            fold_gray = np.asarray(gray["fold_gray_thresholds"], dtype=float)
            t_low_folds = fold_gray[:, 0]
            t_high_folds = fold_gray[:, 1]
            # Coherent medoid pair from the CV (see gray_zone["t_low_final"]).
            # Independent medians of the two arrays are NOT interchangeable with
            # it: they can describe a band no fold produced. Fall back to them
            # only for result dicts produced before that key existed.
            t_low_final = float(gray.get("t_low_final", np.nan))
            t_high_final = float(gray.get("t_high_final", np.nan))
            if not (np.isfinite(t_low_final) and np.isfinite(t_high_final)):
                t_low_final = float(np.median(t_low_folds))
                t_high_final = float(np.median(t_high_folds))
        else:
            # fallback from medians if present
            t_low_final = float(gray.get("median_t_low_across_folds", np.nan))
            t_high_final = float(gray.get("median_t_high_across_folds", np.nan))
            t_low_folds = np.array([t_low_final], dtype=float)
            t_high_folds = np.array([t_high_final], dtype=float)

        # targets for bootstrap re-estimation if available
        target_sens = gray.get("targets", {}).get("target_sensitivity_for_t_low", 0.90)
        target_spec = gray.get("targets", {}).get("target_specificity_for_t_high", 0.90)
    else:
        t_low_final = np.nan
        t_high_final = np.nan
        t_low_folds = np.array([], dtype=float)
        t_high_folds = np.array([], dtype=float)
        target_sens, target_spec = 0.90, 0.90

    # -------------------------------------------------------------------------
    # Fit calibration on OOF reference set (logit)
    # -------------------------------------------------------------------------
    logit_model = _fit_calibration(s_ref, y_ref)
    p_query = logit_model.predict_proba(scores_query_oriented.reshape(-1, 1))[:, 1]

    # The operating point: the calibrated probability AT the binary threshold.
    # This, not 0.5, is where the decision actually flips.
    p_threshold = float(logit_model.predict_proba(np.array([[thr_final]]))[0, 1])

    # Direct outputs (point estimates)
    pred_binary = (scores_query_oriented >= thr_final).astype(int)
    pred_label = np.where(pred_binary == 1, "responder", "non_responder")

    # Invariant: the calibration is monotone increasing in the oriented score,
    # so "score >= thr_final" and "p >= p_threshold" must be the same statement.
    # If this trips, the calibration has learned a negative slope (orientation
    # inference wrong, or degenerate reference set) and every probability-derived
    # column below would contradict the prediction it is supposed to describe.
    _tol = 1e-9
    _viol = (((p_query > p_threshold + _tol) & (pred_binary == 0)) |
             ((p_query < p_threshold - _tol) & (pred_binary == 1)))
    if np.any(_viol):
        _i = np.where(_viol)[0]
        raise RuntimeError(
            "Calibrated probability contradicts the binary rule -- the logistic "
            "calibration is not monotone in the oriented score.\n"
            f"  threshold (oriented)        : {thr_final!r}\n"
            f"  p at threshold              : {p_threshold!r}\n"
            f"  logit coefficient (scaled)  : {float(logit_model[-1].coef_.ravel()[0])!r}\n"
            f"  reference n / n_pos / n_neg : {len(y_ref)} / {int((y_ref == 1).sum())} / {int((y_ref == 0).sum())}\n"
            f"  offending samples           : {[str(sample_ids[i]) for i in _i[:10]]}\n"
            f"  their oriented scores       : {[float(scores_query_oriented[i]) for i in _i[:10]]}\n"
            f"  their p                     : {[float(p_query[i]) for i in _i[:10]]}\n"
            f"  their pred                  : {[int(pred_binary[i]) for i in _i[:10]]}"
        )

    margin = scores_query_oriented - thr_final  # oriented margin

    # Confidence language is measured from the DECISION BOUNDARY, not from 0.5.
    # Previously confidence = max(p, 1-p) and margin_prob = 2|p - 0.5|, which
    # answer "how far from a coin flip" -- a question nobody asked, and one whose
    # answer disagrees with the prediction whenever p_threshold != 0.5.
    p_centred = _recentre_prob(p_query, p_threshold)
    confidence = np.maximum(p_centred, 1 - p_centred)
    ambiguity = 1 - confidence
    margin_prob = 2 * np.abs(p_centred - 0.5)

    # These stay TRUE calibrated probabilities and are deliberately NOT
    # recentred: "probability this call is wrong" is only meaningful on the real
    # probability scale, and a rescaled version would look like a probability
    # without being one. conditional_error_risk CAN therefore exceed 0.5 -- that
    # is a real property of the Youden operating point on an imbalanced cohort,
    # not an inconsistency, and pred_against_majority_prob below marks it.
    # Only the branch matching the actual prediction is a fact about this
    # sample; the other is "what the risk would have been had we called the
    # other class". Both used to be emitted on every row under names that read
    # as factual, so half of each column was a counterfactual. Report the
    # applicable one and leave the other NaN; conditional_error_risk carries the
    # one that applies, whichever it is.
    risk_fp_if_pred_positive = np.where(pred_binary == 1, 1 - p_query, np.nan)
    risk_fn_if_pred_negative = np.where(pred_binary == 0, p_query, np.nan)
    conditional_error_risk = np.where(pred_binary == 1, 1 - p_query, p_query)

    # Rows where the threshold rule calls a class that the calibrated
    # probability puts in the minority (e.g. "responder" at p = 0.35). Legitimate
    # -- Youden maximises sens+spec, not accuracy -- but the reader must be told.
    pred_against_majority_prob = (((pred_binary == 1) & (p_query < 0.5)) |
                                  ((pred_binary == 0) & (p_query >= 0.5)))

    # Gray zone point estimates
    if has_gray and np.isfinite(t_low_final) and np.isfinite(t_high_final):
        pred_gray = _gray_zone_labels(scores_query_oriented, t_low_final, t_high_final)
        gray_label = np.array([_gray_label(pg) for pg in pred_gray], dtype=object)
        gray_distance = np.array([_signed_gray_distance(s, t_low_final, t_high_final) for s in scores_query_oriented], dtype=float)
    else:
        pred_gray = np.full(n_query, -999, dtype=int)
        gray_label = np.array(["not_available"] * n_query, dtype=object)
        gray_distance = np.full(n_query, np.nan, dtype=float)

    # The binary threshold and the gray band are estimated by independent
    # criteria on independent quantities, so nothing forces them to agree: a row
    # can read pred_label = "responder" and pred_gray_label =
    # "non_responder_confident" at the same time. Flag it per sample rather than
    # letting the reader discover the contradiction by eye.
    binary_gray_conflict = (pred_gray >= 0) & (pred_gray != pred_binary)

    # Threshold stability (CV-based from fold thresholds)
    # q = P(score >= T)
    # threshold_stability_cv / pred_stability_cv / gray_confident_cv USED TO BE
    # COMPUTED HERE and are deliberately gone.
    #
    # They measured agreement across fold_thresholds, but any two fold-training
    # sets share ~94% of their samples, so the fold thresholds are near-duplicates
    # by construction. Measured over a 200-point sweep of the score range:
    #
    #                            distinct values   pinned at exactly 0 or 1
    #   threshold_stability_cv        2 - 4              89% - 97%
    #   threshold_stability_boot     11 - 13             65% - 73%
    #   gray_confident_cv             4 - 5              82% - 96%
    #   gray_confident_boot          14 - 16             31% - 66%
    #
    # 17 folds produced only TWO distinct thresholds. On a 20R/20NR cohort the
    # CV statistic was strictly between 0.05 and 0.95 over 0% of the observed
    # score range -- a binary indicator wearing the costume of a proportion. And
    # where the bootstrap said genuinely uncertain (0.07-0.83), the CV version
    # read 0.00-1.00, i.e. it claimed certainty precisely where there was none.
    #
    # Sat next to a calibrated probability it read as independent corroboration.
    # It was not: it is a degenerate approximation of the *_boot columns, which
    # refit the threshold on real resamples and are kept. Do not reintroduce
    # these without re-running that comparison.

    # -------------------------------------------------------------------------
    # Bootstrap: calibration + threshold(s) + CIs
    # -------------------------------------------------------------------------
    idx_all = np.arange(len(y_ref))
    idx_pos = np.where(y_ref == 1)[0]
    idx_neg = np.where(y_ref == 0)[0]

    n_bootstrap = int(n_bootstrap)
    # Store bootstrap replicate outputs (per query sample)
    p_boot = np.full((n_bootstrap, n_query), np.nan, dtype=float)
    thr_boot = np.full(n_bootstrap, np.nan, dtype=float)
    # Operating point per replicate: the threshold and the calibration are both
    # refit on each resample, so the point the confidence is measured from moves
    # too. Using the point-estimate p_threshold here would mix a bootstrapped
    # numerator with a fixed centre and understate the spread.
    p_thr_boot = np.full(n_bootstrap, np.nan, dtype=float)

    if has_gray:
        t_low_boot = np.full(n_bootstrap, np.nan, dtype=float)
        t_high_boot = np.full(n_bootstrap, np.nan, dtype=float)
    else:
        t_low_boot = None
        t_high_boot = None

    # One bootstrap replicate (refit logistic calibration + threshold(s)).
    # Each replicate gets its own seed so the result is independent of the
    # execution order (reproducible across runs and worker counts).
    def _one_conf_boot(child_seed):
        rng_b = np.random.default_rng(child_seed)
        if bootstrap_stratified and (len(idx_pos) > 0) and (len(idx_neg) > 0):
            boot_pos = rng_b.choice(idx_pos, size=len(idx_pos), replace=True)
            boot_neg = rng_b.choice(idx_neg, size=len(idx_neg), replace=True)
            boot_idx = np.concatenate([boot_pos, boot_neg])
            rng_b.shuffle(boot_idx)
        else:
            boot_idx = rng_b.choice(idx_all, size=len(idx_all), replace=True)

        yb = y_ref[boot_idx]
        sb = s_ref[boot_idx]

        p_vec = np.full(n_query, np.nan, dtype=float)
        thr_v = np.nan
        tl_v = np.nan
        th_v = np.nan
        p_thr_v = np.nan

        # Need both classes for ROC/logit
        model_b = None
        if len(np.unique(yb)) >= 2:
            try:
                model_b = _fit_calibration(sb, yb)
                p_vec = model_b.predict_proba(scores_query_oriented.reshape(-1, 1))[:, 1]
            except Exception:
                model_b = None
            try:
                thr_v, _ = _pick_single_threshold_from_roc(yb, sb, criterion=threshold_criterion)
            except Exception:
                thr_v = np.nan
            if (model_b is not None) and np.isfinite(thr_v):
                try:
                    p_thr_v = float(model_b.predict_proba(np.array([[thr_v]]))[0, 1])
                except Exception:
                    p_thr_v = np.nan
            if has_gray:
                try:
                    tl_v, th_v, _ = _pick_gray_thresholds_from_roc(
                        yb, sb, target_sens=target_sens, target_spec=target_spec)
                except Exception:
                    tl_v, th_v = np.nan, np.nan
        return p_vec, thr_v, tl_v, th_v, p_thr_v

    _child_seeds = np.random.SeedSequence(bootstrap_seed).spawn(n_bootstrap)
    _n = _n_jobs()
    if _n == 1:
        _boot_results = [_one_conf_boot(s) for s in _child_seeds]
    else:
        try:
            _boot_results = Parallel(n_jobs=_n, prefer="threads")(
                delayed(_one_conf_boot)(s) for s in _child_seeds
            )
        except Exception:
            _boot_results = [_one_conf_boot(s) for s in _child_seeds]

    for b, (_p_vec, _thr_v, _tl_v, _th_v, _p_thr_v) in enumerate(_boot_results):
        p_boot[b, :] = _p_vec
        thr_boot[b] = _thr_v
        p_thr_boot[b] = _p_thr_v
        if has_gray:
            t_low_boot[b] = _tl_v
            t_high_boot[b] = _th_v

    # Derived bootstrap arrays. Each replicate is recentred on ITS OWN
    # operating point, matching how the point estimates above are built.
    p_centred_boot = _recentre_prob(p_boot, p_thr_boot.reshape(-1, 1))
    confidence_boot = np.maximum(p_centred_boot, 1 - p_centred_boot)
    ambiguity_boot = 1 - confidence_boot
    margin_prob_boot = 2 * np.abs(p_centred_boot - 0.5)

    # Masked like the point estimates above, so the CI of a counterfactual is
    # NaN rather than a plausible-looking interval.
    _is_pos = (pred_binary.reshape(1, -1) == 1)
    risk_fp_boot = np.where(_is_pos, 1 - p_boot, np.nan)
    risk_fn_boot = np.where(_is_pos, np.nan, p_boot)
    # conditional risk based on FIXED point-estimate prediction
    cond_err_boot = np.where(_is_pos, 1 - p_boot, p_boot)

    # threshold stability (bootstrap-based)
    # q_boot = P(score >= T_boot) estimated by mean indicator across bootstrap thresholds
    valid_thr_mask = np.isfinite(thr_boot)
    thr_boot_valid = thr_boot[valid_thr_mask]
    if thr_boot_valid.size > 0:
        threshold_stability_boot = np.array(
            [np.mean(s >= thr_boot_valid) for s in scores_query_oriented],
            dtype=float
        )
        pred_stability_boot = np.where(pred_binary == 1, threshold_stability_boot, 1 - threshold_stability_boot)
        margin_to_thr_boot = scores_query_oriented.reshape(1, -1) - thr_boot_valid.reshape(-1, 1)
    else:
        threshold_stability_boot = np.full(n_query, np.nan)
        pred_stability_boot = np.full(n_query, np.nan)
        margin_to_thr_boot = np.full((1, n_query), np.nan)

    # gray-zone bootstrap stability
    if has_gray and (t_low_boot is not None) and (t_high_boot is not None):
        valid_gray_mask = np.isfinite(t_low_boot) & np.isfinite(t_high_boot)
        t_low_valid = t_low_boot[valid_gray_mask]
        t_high_valid = t_high_boot[valid_gray_mask]
        if t_low_valid.size > 0:
            gray_pred_boot = np.full((t_low_valid.size, n_query), -1, dtype=int)
            for i, s in enumerate(scores_query_oriented):
                # Same rule as the point estimate: scalar score against the
                # arrays of bootstrap thresholds (broadcasts).
                gray_pred_boot[:, i] = _gray_zone_labels(s, t_low_valid, t_high_valid)

            gray_confident_boot = np.mean(gray_pred_boot != -1, axis=0)
            gray_same_label_boot = np.mean(gray_pred_boot == pred_gray.reshape(1, -1), axis=0)
            # pred_gray == -999 means the point-estimate zone was unavailable, so
            # there is no label for the replicates to agree WITH. The comparison
            # above then never matches and reports 0.0, which is indistinguishable
            # from "the replicates genuinely disagree every time". Must be NaN.
            gray_same_label_boot = np.where(pred_gray == -999, np.nan, gray_same_label_boot)
        else:
            gray_confident_boot = np.full(n_query, np.nan)
            gray_same_label_boot = np.full(n_query, np.nan)
    else:
        gray_confident_boot = np.full(n_query, np.nan)
        gray_same_label_boot = np.full(n_query, np.nan)

    # -------------------------------------------------------------------------
    # Build CI table per sample
    # -------------------------------------------------------------------------
    alpha = 1 - ci_level

    rows = []
    for i, sid in enumerate(sample_ids):
        # point metrics
        row = {
            "sample_id": sid,
            "score_raw": float(scores_raw_query[i]),
            "score_oriented": float(scores_query_oriented[i]),

            # binary pred
            "threshold_used": float(thr_final),
            "pred": int(pred_binary[i]),
            "pred_label": pred_label[i],
            "margin": float(margin[i]),

            # calibrated probability + derived
            "p": float(p_query[i]),  # alias
            "p_calibrated_logit": float(p_query[i]),
            # The operating point: p at the decision threshold. Same for every
            # row, but it is the reference the three columns below are measured
            # from, so it belongs next to them rather than buried in a summary.
            "p_at_threshold": float(p_threshold),
            "confidence": float(confidence[i]),
            "ambiguity": float(ambiguity[i]),
            "margin_prob": float(margin_prob[i]),
            "pred_against_majority_prob": bool(pred_against_majority_prob[i]),

            # risks
            "risk_FP_if_pred_positive": float(risk_fp_if_pred_positive[i]),
            "risk_FN_if_pred_negative": float(risk_fn_if_pred_negative[i]),
            "conditional_error_risk": float(conditional_error_risk[i]),

            # threshold stability (bootstrap only -- the CV-based versions were
            # degenerate by construction, see the note above)
            "threshold_stability_boot": float(threshold_stability_boot[i]),
            "pred_stability_boot": float(pred_stability_boot[i]),

            # gray zone
            "gray_t_low": float(t_low_final) if np.isfinite(t_low_final) else np.nan,
            "gray_t_high": float(t_high_final) if np.isfinite(t_high_final) else np.nan,
            "pred_gray": int(pred_gray[i]) if pred_gray[i] != -999 else np.nan,
            "pred_gray_label": gray_label[i],
            "binary_gray_conflict": bool(binary_gray_conflict[i]),
            "gray_distance_to_boundary": float(gray_distance[i]),
            "gray_confident_boot": float(gray_confident_boot[i]),
            "gray_label_stability_boot": float(gray_same_label_boot[i]),
        }

        # Bootstrap CIs (probability + derived)
        p_ci = _percentile_ci(p_boot[:, i], ci_level=ci_level)
        conf_ci = _percentile_ci(confidence_boot[:, i], ci_level=ci_level)
        amb_ci = _percentile_ci(ambiguity_boot[:, i], ci_level=ci_level)
        mp_ci = _percentile_ci(margin_prob_boot[:, i], ci_level=ci_level)
        rfp_ci = _percentile_ci(risk_fp_boot[:, i], ci_level=ci_level)
        rfn_ci = _percentile_ci(risk_fn_boot[:, i], ci_level=ci_level)
        cer_ci = _percentile_ci(cond_err_boot[:, i], ci_level=ci_level)

        # Threshold/margin CI from bootstrap thresholds
        if margin_to_thr_boot.size > 0:
            margin_thr_ci = _percentile_ci(margin_to_thr_boot[:, i], ci_level=ci_level)
        else:
            margin_thr_ci = (np.nan, np.nan)

        row.update({
            f"p_ci_low_{int(ci_level*100)}": p_ci[0],
            f"p_ci_high_{int(ci_level*100)}": p_ci[1],

            f"confidence_ci_low_{int(ci_level*100)}": conf_ci[0],
            f"confidence_ci_high_{int(ci_level*100)}": conf_ci[1],

            f"ambiguity_ci_low_{int(ci_level*100)}": amb_ci[0],
            f"ambiguity_ci_high_{int(ci_level*100)}": amb_ci[1],

            f"margin_prob_ci_low_{int(ci_level*100)}": mp_ci[0],
            f"margin_prob_ci_high_{int(ci_level*100)}": mp_ci[1],

            f"risk_FP_if_pred_positive_ci_low_{int(ci_level*100)}": rfp_ci[0],
            f"risk_FP_if_pred_positive_ci_high_{int(ci_level*100)}": rfp_ci[1],

            f"risk_FN_if_pred_negative_ci_low_{int(ci_level*100)}": rfn_ci[0],
            f"risk_FN_if_pred_negative_ci_high_{int(ci_level*100)}": rfn_ci[1],

            f"conditional_error_risk_ci_low_{int(ci_level*100)}": cer_ci[0],
            f"conditional_error_risk_ci_high_{int(ci_level*100)}": cer_ci[1],

            f"margin_vs_boot_threshold_ci_low_{int(ci_level*100)}": margin_thr_ci[0],
            f"margin_vs_boot_threshold_ci_high_{int(ci_level*100)}": margin_thr_ci[1],
        })

        rows.append(row)

    report = pd.DataFrame(rows).set_index("sample_id")

    # -------------------------------------------------------------------------
    # Global bootstrap summary (useful debug / reproducibility)
    # -------------------------------------------------------------------------
    bootstrap_summary = {
        "n_bootstrap_requested": int(n_bootstrap),
        "n_bootstrap_valid_threshold": int(np.sum(np.isfinite(thr_boot))),
        "n_bootstrap_valid_gray": int(np.sum(np.isfinite(t_low_boot) & np.isfinite(t_high_boot))) if has_gray else 0,
        "threshold_boot_ci": {
            "ci_low": _percentile_ci(thr_boot, ci_level=ci_level)[0],
            "ci_high": _percentile_ci(thr_boot, ci_level=ci_level)[1],
        },
        "p_at_threshold_boot_ci": {
            "ci_low": _percentile_ci(p_thr_boot, ci_level=ci_level)[0],
            "ci_high": _percentile_ci(p_thr_boot, ci_level=ci_level)[1],
        },
    }

    if has_gray:
        bootstrap_summary["gray_t_low_boot_ci"] = {
            "ci_low": _percentile_ci(t_low_boot, ci_level=ci_level)[0],
            "ci_high": _percentile_ci(t_low_boot, ci_level=ci_level)[1],
        }
        bootstrap_summary["gray_t_high_boot_ci"] = {
            "ci_low": _percentile_ci(t_high_boot, ci_level=ci_level)[0],
            "ci_high": _percentile_ci(t_high_boot, ci_level=ci_level)[1],
        }

    reference_info = {
        "orientation_sign_raw_to_oriented": float(sign),  # +1 or -1
        "score_orientation_text": cv_result.get("score_orientation", "unknown"),
        "reference_n": int(len(y_ref)),
        "reference_n_pos": int(np.sum(y_ref == 1)),
        "reference_n_neg": int(np.sum(y_ref == 0)),
        "threshold_final_median_fold": float(thr_final),
        "p_at_threshold": float(p_threshold),
        "n_pred_against_majority_prob": int(np.sum(pred_against_majority_prob)),
        "calibration_C": float(calibration_C),
        "calibration_coef_standardised": float(logit_model[-1].coef_.ravel()[0]),
        "calibration_n_iter": int(np.ravel(logit_model[-1].n_iter_)[0]),
        "gray_zone_available": bool(has_gray),
        "gray_t_low_final_median_fold": float(t_low_final) if np.isfinite(t_low_final) else np.nan,
        "gray_t_high_final_median_fold": float(t_high_final) if np.isfinite(t_high_final) else np.nan,
        "gray_zone_inverted": bool(np.isfinite(t_low_final) and np.isfinite(t_high_final)
                                   and t_low_final > t_high_final),
        "threshold_inside_gray_zone": bool(np.isfinite(t_low_final) and np.isfinite(t_high_final)
                                           and t_low_final <= thr_final <= t_high_final),
        "n_binary_gray_conflicts": int(np.sum(binary_gray_conflict)),
    }

    return {
        "report": report,
        "bootstrap_summary": bootstrap_summary,
        "reference_info": reference_info,
    }

