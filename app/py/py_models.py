# contains : 
# nested_cv_signature() : Evaluates a signature for a specific model
# sample_confidence_report() : Evaluates the confidence we can have in a specific prediction


import os
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc, confusion_matrix
from joblib import Parallel, delayed


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



def nested_cv_signature(
    df_scores: pd.DataFrame,
    sample_ID_responders,
    sample_ID_non_responders,
    score_col: str = None,
    # CV
    n_outer: int = 5,
    n_inner: int = 5,
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
    # Outputs
    return_oof_predictions: bool = True,
):
    """
    Nested CV for a fixed per-sample score (1D) with:
      - single-threshold classification (ROC-derived threshold selected in inner CV)
      - optional gray zone (double threshold)
      - OOF ROC/AUC, confusion matrices, Se/Sp/PPV/NPV
      - bootstrap confidence intervals

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
    def _safe_div(num, den):
        return np.nan if den == 0 else num / den

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

    def _pick_single_threshold_from_roc(y_true, scores, criterion="youden"):
        """
        scores oriented so that higher = more likely responder.
        Returns chosen threshold + ROC arrays.
        """
        fpr, tpr, thresholds = roc_curve(y_true, scores)

        if criterion == "youden":
            # J = sensitivity + specificity -1 = tpr - fpr
            j = tpr - fpr
            idx = int(np.nanargmax(j))
        elif criterion == "closest_topleft":
            d = np.sqrt((fpr - 0.0) ** 2 + (tpr - 1.0) ** 2)
            idx = int(np.nanargmin(d))
        else:
            raise ValueError("threshold_criterion must be 'youden' or 'closest_topleft'.")

        return float(thresholds[idx]), {
            "fpr": fpr,
            "tpr": tpr,
            "thresholds": thresholds,
            "best_idx": idx,
        }

    def _pick_gray_thresholds_from_roc(
        y_true, scores,
        target_sens=0.90, target_spec=0.90
    ):
        """
        Returns (t_low, t_high) based on ROC on oriented scores (higher => responder):
          - t_low: highest threshold achieving target sensitivity
                   -> scores <= t_low are "confident negatives" (few FN expected)
          - t_high: lowest threshold achieving target specificity
                   -> scores >= t_high are "confident positives" (few FP expected)

        In practice:
          score <= t_low  => class 0 (non-responder)
          score >= t_high => class 1 (responder)
          else            => uncertain (-1)
        """
        fpr, tpr, thresholds = roc_curve(y_true, scores)
        spec = 1 - fpr

        # t_low from target sensitivity
        idx_sens_valid = np.where(tpr >= target_sens)[0]
        if len(idx_sens_valid) == 0:
            # Fallback: pick threshold with max sensitivity (usually very low threshold)
            idx_low = int(np.nanargmax(tpr))
        else:
            # Among thresholds with target sensitivity, maximize specificity
            idx_low = idx_sens_valid[int(np.nanargmax(spec[idx_sens_valid]))]

        # t_high from target specificity
        idx_spec_valid = np.where(spec >= target_spec)[0]
        if len(idx_spec_valid) == 0:
            # Fallback: pick threshold with max specificity (usually very high threshold)
            idx_high = int(np.nanargmax(spec))
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

    def _apply_gray_zone(scores, t_low, t_high):
        """
        Returns pred in {-1,0,1}:
          0 if score <= t_low
          1 if score >= t_high
         -1 otherwise (uncertain)
        """
        pred = np.full(len(scores), -1, dtype=int)
        pred[scores <= t_low] = 0
        pred[scores >= t_high] = 1

        # If thresholds overlap inversely (t_low > t_high), all points could be assigned
        # by both conditions; positives overwrite negatives above.
        # We handle this with a flag elsewhere.
        return pred

    def _percentile_ci(values, ci_level=0.95):
        arr = np.asarray(values, dtype=float)
        arr = arr[~np.isnan(arr)]
        if arr.size == 0:
            return (np.nan, np.nan)
        alpha = 1 - ci_level
        lo = np.quantile(arr, alpha / 2)
        hi = np.quantile(arr, 1 - alpha / 2)
        return float(lo), float(hi)

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

    def _dict_to_ci_table(point_estimates: dict, ci_dict: dict, metric_order=None):
        rows = []
        if metric_order is None:
            metric_order = list(point_estimates.keys())
        for m in metric_order:
            if m not in point_estimates:
                continue
            pe = point_estimates[m]
            ci = ci_dict.get(m, {"ci_low": np.nan, "ci_high": np.nan})
            rows.append({
                "metric": m,
                "estimate": pe,
                "ci_low": ci["ci_low"],
                "ci_high": ci["ci_high"],
            })
        return pd.DataFrame(rows)

    # -------------------------------------------------------------------------
    # 1) Validation + data prep
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
        ex = list(overlap)[:10]
        raise ValueError(f"Some IDs are in both classes (e.g., {ex}).")

    all_ids = list(responders | non_responders)
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

    X = df["score"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=int)
    sample_ids = df.index.to_numpy()

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos < 2 or n_neg < 2:
        raise ValueError("Need at least 2 samples per class.")

    n_outer_eff = min(n_outer, n_pos, n_neg)
    if n_outer_eff < 2:
        raise ValueError("Effective n_outer < 2 (not enough samples per class).")

    # -------------------------------------------------------------------------
    # 2) Nested CV
    # -------------------------------------------------------------------------
    outer_cv = StratifiedKFold(
        n_splits=int(n_outer_eff), shuffle=shuffle, random_state=random_state
    )

    # OOF containers
    oof_score = np.full(len(y), np.nan, dtype=float)  # oriented score
    oof_pred = np.full(len(y), -1, dtype=int)         # standard binary pred
    oof_thr = np.full(len(y), np.nan, dtype=float)    # threshold used in that outer fold

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

    for fold_id, (outer_train_idx, outer_test_idx) in enumerate(outer_cv.split(X, y), start=1):
        X_train, y_train = X[outer_train_idx], y[outer_train_idx]
        X_test, y_test = X[outer_test_idx], y[outer_test_idx]

        # Inner CV for threshold stability (collect thresholds across inner train folds, aggregate by median)
        n_pos_train = int((y_train == 1).sum())
        n_neg_train = int((y_train == 0).sum())
        n_inner_eff = min(n_inner, n_pos_train, n_neg_train)

        inner_thresholds_single = []
        inner_thresholds_low = []
        inner_thresholds_high = []

        inner_mode = "nested_inner_cv"
        if n_inner_eff < 2:
            # fallback: choose directly on outer train
            inner_mode = "fallback_outer_train_only"
            thr_single, _ = _pick_single_threshold_from_roc(y_train, X_train, criterion=threshold_criterion)
            inner_thresholds_single = [thr_single]

            if use_gray_zone:
                t_low, t_high, _ = _pick_gray_thresholds_from_roc(
                    y_train, X_train,
                    target_sens=gray_target_sensitivity,
                    target_spec=gray_target_specificity
                )
                inner_thresholds_low = [t_low]
                inner_thresholds_high = [t_high]
        else:
            inner_cv = StratifiedKFold(
                n_splits=int(n_inner_eff),
                shuffle=shuffle,
                random_state=random_state + fold_id
            )
            for inner_train_idx, inner_val_idx in inner_cv.split(X_train, y_train):
                # Threshold is estimated from inner TRAIN only (clean logic)
                Xi_train = X_train[inner_train_idx]
                yi_train = y_train[inner_train_idx]

                thr_i, _ = _pick_single_threshold_from_roc(
                    yi_train, Xi_train, criterion=threshold_criterion
                )
                inner_thresholds_single.append(thr_i)

                if use_gray_zone:
                    t_low_i, t_high_i, _ = _pick_gray_thresholds_from_roc(
                        yi_train, Xi_train,
                        target_sens=gray_target_sensitivity,
                        target_spec=gray_target_specificity
                    )
                    inner_thresholds_low.append(t_low_i)
                    inner_thresholds_high.append(t_high_i)

        # Aggregate inner thresholds -> outer-fold threshold(s)
        thr_outer = float(np.median(inner_thresholds_single))
        fold_thresholds.append(thr_outer)

        # Apply standard binary threshold on outer test
        y_pred_test = (X_test >= thr_outer).astype(int)

        # Store OOF
        oof_score[outer_test_idx] = X_test
        oof_pred[outer_test_idx] = y_pred_test
        oof_thr[outer_test_idx] = thr_outer

        # Gray zone
        gray_info = {}
        if use_gray_zone:
            t_low_outer = float(np.median(inner_thresholds_low))
            t_high_outer = float(np.median(inner_thresholds_high))

            pred_gray_test = _apply_gray_zone(X_test, t_low_outer, t_high_outer)

            oof_pred_gray[outer_test_idx] = pred_gray_test
            oof_t_low[outer_test_idx] = t_low_outer
            oof_t_high[outer_test_idx] = t_high_outer

            fold_gray_thresholds.append((t_low_outer, t_high_outer))
            gray_info = {
                "t_low_outer": t_low_outer,
                "t_high_outer": t_high_outer,
                "gray_zone_inverted_order": bool(t_low_outer > t_high_outer),
                "inner_thresholds_low": [float(t) for t in inner_thresholds_low],
                "inner_thresholds_high": [float(t) for t in inner_thresholds_high],
            }

        fold_details.append({
            "fold": fold_id,
            "n_train": int(len(outer_train_idx)),
            "n_test": int(len(outer_test_idx)),
            "n_pos_test": int(np.sum(y_test == 1)),
            "n_neg_test": int(np.sum(y_test == 0)),
            "threshold_outer": thr_outer,
            "inner_thresholds_single": [float(t) for t in inner_thresholds_single],
            "inner_mode": inner_mode,
            **gray_info
        })

    # Basic check
    if np.isnan(oof_score).any() or np.any(oof_pred < 0):
        raise RuntimeError("OOF predictions are incomplete.")

    if use_gray_zone and np.any(oof_pred_gray == -999):
        raise RuntimeError("OOF gray-zone predictions are incomplete.")

    # -------------------------------------------------------------------------
    # 3) OOF ROC + standard confusion matrix + standard metrics
    # -------------------------------------------------------------------------
    fpr, tpr, roc_thresholds = roc_curve(y, oof_score)
    roc_auc = float(auc(fpr, tpr))

    standard_metrics, cm = _compute_binary_metrics(y, oof_pred)
    cm_df = pd.DataFrame(
        cm,
        index=["True_NonResponder(0)", "True_Responder(1)"],
        columns=["Pred_NonResponder(0)", "Pred_Responder(1)"]
    )

    # Bootstrap CI for standard metrics (+ AUC)
    standard_ci = _bootstrap_ci_binary_metrics(
        y_true=y,
        y_pred=oof_pred,
        y_score=oof_score,
        n_boot=int(n_bootstrap),
        ci_level=ci_level,
        seed=bootstrap_seed,
        stratified=bootstrap_stratified
    )
    standard_ci["auc"] = standard_ci.get("auc", {"ci_low": np.nan, "ci_high": np.nan})

    standard_metric_estimates = {**standard_metrics, "auc": roc_auc}
    standard_metrics_table = _dict_to_ci_table(
        standard_metric_estimates,
        standard_ci,
        metric_order=[
            "auc", "sensitivity", "specificity", "ppv", "npv",
            "fpr", "fnr", "accuracy", "balanced_accuracy",
            "risk_false_positive_given_pred_positive",
            "risk_false_negative_given_pred_negative",
            "prevalence"
        ]
    )

    # prevalence is deterministic here (bootstrap CI omitted unless you want to add)
    # fill prevalence CI with NaN (already)
    # -------------------------------------------------------------------------
    # 4) Gray zone outputs
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

        confident_mask = (oof_pred_gray != -1)
        n_confident = int(np.sum(confident_mask))
        coverage = _safe_div(n_confident, len(y))
        rejection_rate = np.nan if np.isnan(coverage) else 1 - coverage

        if n_confident > 0:
            y_conf = y[confident_mask]
            p_conf = oof_pred_gray[confident_mask].astype(int)  # now only 0/1
            s_conf = oof_score[confident_mask]

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

            gray_estimates = {
                **gray_metrics_conf,
                "coverage": coverage,
                "rejection_rate": rejection_rate,
                "auc": auc_conf
            }

            gray_metrics_table = _dict_to_ci_table(
                gray_estimates,
                gray_ci,
                metric_order=[
                    "coverage", "rejection_rate",
                    "auc", "sensitivity", "specificity", "ppv", "npv",
                    "fpr", "fnr", "accuracy", "balanced_accuracy",
                    "risk_false_positive_given_pred_positive",
                    "risk_false_negative_given_pred_negative"
                ]
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
            gray_ci = {}
            gray_metrics_table = pd.DataFrame(
                columns=["metric", "estimate", "ci_low", "ci_high"]
            )

        # Summary of gray thresholds across folds
        fold_gray_df = pd.DataFrame(fold_gray_thresholds, columns=["t_low", "t_high"]) \
            if len(fold_gray_thresholds) > 0 else pd.DataFrame(columns=["t_low", "t_high"])

        gray_zone_results = {
            "targets": {
                "target_sensitivity_for_t_low": gray_target_sensitivity,
                "target_specificity_for_t_high": gray_target_specificity,
            },
            "fold_gray_thresholds": fold_gray_thresholds,
            "fold_gray_thresholds_df": fold_gray_df,
            "median_t_low_across_folds": float(np.median(fold_gray_df["t_low"])) if len(fold_gray_df) else np.nan,
            "median_t_high_across_folds": float(np.median(fold_gray_df["t_high"])) if len(fold_gray_df) else np.nan,
            "three_way_table": three_way_table,
            "confident_subset_confusion_matrix": cm_gray_conf,
            "confident_subset_confusion_matrix_df": cm_gray_conf_df,
            "confident_subset_metrics": gray_metrics_conf,
            "confident_subset_bootstrap_ci": gray_ci,
            "confident_subset_metrics_table": gray_metrics_table,
            "coverage": coverage,
            "rejection_rate": rejection_rate,
        }

    # -------------------------------------------------------------------------
    # 5) OOF dataframe
    # -------------------------------------------------------------------------
    oof_df = None
    if return_oof_predictions:
        oof_df = pd.DataFrame({
            "sample_id": sample_ids,
            "y_true": y,
            "score_raw_original": df.loc[sample_ids, "score_raw"].values,
            "score_used_oriented": oof_score,
            "pred_binary_oof": oof_pred,
            "threshold_binary_outer_fold": oof_thr,
        }).set_index("sample_id")

        if use_gray_zone:
            oof_df["pred_gray_oof"] = oof_pred_gray  # -1 uncertain, 0, 1
            oof_df["t_low_outer_fold"] = oof_t_low
            oof_df["t_high_outer_fold"] = oof_t_high
            oof_df["is_confident_gray"] = (oof_df["pred_gray_oof"] != -1)

    # -------------------------------------------------------------------------
    # 6) Return
    # -------------------------------------------------------------------------
    results = {
        "n_samples": int(len(y)),
        "n_positive": int(np.sum(y == 1)),
        "n_negative": int(np.sum(y == 0)),
        "score_orientation": orientation_txt,

        # ROC / AUC (OOF)
        "fpr": fpr,
        "tpr": tpr,
        "roc_thresholds": roc_thresholds,
        "auc": roc_auc,

        # Standard binary classification (OOF)
        "confusion_matrix": cm,
        "confusion_matrix_df": cm_df,
        "metrics": standard_metrics,
        "bootstrap_ci": standard_ci,
        "metrics_table": standard_metrics_table,

        # Thresholds / folds
        "fold_thresholds": fold_thresholds,
        "fold_details": fold_details,
    }

    if use_gray_zone:
        results["gray_zone"] = gray_zone_results

    if return_oof_predictions:
        results["oof"] = oof_df

    return results


def sample_confidence_report(
    scores,
    res_nested_cv_signature: dict,
    sample_ids=None,
    threshold_criterion: str = "youden",  # same logic as v2
    use_gray_zone: bool = True,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    bootstrap_seed: int = 123,
    bootstrap_stratified: bool = True,
):
    """
    Build a per-sample confidence report using outputs from nested_cv_signature.

    Parameters
    ----------
    scores : scalar | list | np.ndarray | pd.Series
        Raw scores (same scale as df_scores[score_col] used in v2).
    res_nested_cv_signature : dict
        Output of nested_cv_signature(...), ideally with res_nested_cv_signature["oof"] present.
    sample_ids : list-like, optional
        Optional IDs for the query scores. If None, auto-generated.
    threshold_criterion : str
        "youden" or "closest_topleft" for bootstrap threshold re-estimation.
    use_gray_zone : bool
        If True and gray_zone exists in res_nested_cv_signature, returns gray-zone predictions.
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
    - Probabilité calibrée: logistic regression (logit) fitted on OOF scores/labels from res_nested_cv_signature["oof"].
    - Threshold used for binary pred: median(res_nested_cv_signature["fold_thresholds"]).
    - Threshold stability:
        * CV-based: proportion of outer-fold thresholds below the sample score
        * Bootstrap-based: proportion of bootstrap thresholds below the sample score
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

    def _safe_div(a, b):
        return np.nan if b == 0 else a / b

    def _percentile_ci(arr, ci=0.95):
        arr = np.asarray(arr, dtype=float)
        arr = arr[~np.isnan(arr)]
        if arr.size == 0:
            return (np.nan, np.nan)
        alpha = 1 - ci
        return (float(np.quantile(arr, alpha / 2)),
                float(np.quantile(arr, 1 - alpha / 2)))

    def _fit_logit(scores_ref, y_ref):
        """
        Simple logistic calibration p(y=1|score).
        """
        model = LogisticRegression(
            solver="lbfgs",
            C=1e6,           # quasi non-regularized
            max_iter=1000
        )
        model.fit(scores_ref.reshape(-1, 1), y_ref.astype(int))
        return model

    def _pick_threshold_from_scores(y_true, scores_ref, criterion="youden"):
        fpr, tpr, thresholds = roc_curve(y_true, scores_ref)
        if criterion == "youden":
            idx = int(np.nanargmax(tpr - fpr))
        elif criterion == "closest_topleft":
            d = np.sqrt((fpr - 0.0) ** 2 + (tpr - 1.0) ** 2)
            idx = int(np.nanargmin(d))
        else:
            raise ValueError("threshold_criterion must be 'youden' or 'closest_topleft'.")
        return float(thresholds[idx])

    def _pick_gray_thresholds_from_scores(y_true, scores_ref, target_sens, target_spec):
        fpr, tpr, thresholds = roc_curve(y_true, scores_ref)
        spec = 1 - fpr

        # t_low: target sensitivity
        valid_sens = np.where(tpr >= target_sens)[0]
        if len(valid_sens) == 0:
            idx_low = int(np.nanargmax(tpr))
        else:
            idx_low = valid_sens[int(np.nanargmax(spec[valid_sens]))]

        # t_high: target specificity
        valid_spec = np.where(spec >= target_spec)[0]
        if len(valid_spec) == 0:
            idx_high = int(np.nanargmax(spec))
        else:
            idx_high = valid_spec[int(np.nanargmax(tpr[valid_spec]))]

        return float(thresholds[idx_low]), float(thresholds[idx_high])

    def _gray_pred(s, t_low, t_high):
        """
        Returns:
          0 if s <= t_low
          1 if s >= t_high
         -1 otherwise
        """
        if s <= t_low:
            return 0
        if s >= t_high:
            return 1
        return -1

    def _gray_label(pred_gray):
        if pred_gray == 0:
            return "non_responder_confident"
        if pred_gray == 1:
            return "responder_confident"
        return "uncertain"

    def _gray_distance(s, t_low, t_high):
        if s <= t_low:
            return t_low - s
        if s >= t_high:
            return s - t_high
        return min(s - t_low, t_high - s)

    # -------------------------------------------------------------------------
    # Validate res_nested_cv_signature
    # -------------------------------------------------------------------------
    if "oof" not in res_nested_cv_signature:
        raise ValueError(
            "res_nested_cv_signature must contain 'oof'. Please run nested_cv_signature_v2(..., return_oof_predictions=True)."
        )

    oof = res_nested_cv_signature["oof"].copy()
    required_oof_cols = {"y_true", "score_raw_original", "score_used_oriented"}
    missing = required_oof_cols - set(oof.columns)
    if missing:
        raise ValueError(f"res_nested_cv_signature['oof'] is missing columns: {missing}")

    if "fold_thresholds" not in res_nested_cv_signature or len(res_nested_cv_signature["fold_thresholds"]) == 0:
        raise ValueError("res_nested_cv_signature must contain non-empty 'fold_thresholds'.")

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
    fold_thresholds = np.asarray(res_nested_cv_signature["fold_thresholds"], dtype=float)
    thr_final = float(np.median(fold_thresholds))

    has_gray = use_gray_zone and ("gray_zone" in res_nested_cv_signature) and (res_nested_cv_signature["gray_zone"] is not None)
    if has_gray:
        gz = res_nested_cv_signature["gray_zone"]
        # Try robustly to fetch fold gray thresholds
        if "fold_gray_thresholds" in gz and len(gz["fold_gray_thresholds"]) > 0:
            fold_gray = np.asarray(gz["fold_gray_thresholds"], dtype=float)
            t_low_folds = fold_gray[:, 0]
            t_high_folds = fold_gray[:, 1]
            t_low_final = float(np.median(t_low_folds))
            t_high_final = float(np.median(t_high_folds))
        else:
            # fallback from medians if present
            t_low_final = float(gz.get("median_t_low_across_folds", np.nan))
            t_high_final = float(gz.get("median_t_high_across_folds", np.nan))
            t_low_folds = np.array([t_low_final], dtype=float)
            t_high_folds = np.array([t_high_final], dtype=float)

        # targets for bootstrap re-estimation if available
        target_sens = gz.get("targets", {}).get("target_sensitivity_for_t_low", 0.90)
        target_spec = gz.get("targets", {}).get("target_specificity_for_t_high", 0.90)
    else:
        t_low_final = np.nan
        t_high_final = np.nan
        t_low_folds = np.array([], dtype=float)
        t_high_folds = np.array([], dtype=float)
        target_sens, target_spec = 0.90, 0.90

    # -------------------------------------------------------------------------
    # Fit calibration on OOF reference set (logit)
    # -------------------------------------------------------------------------
    logit_model = _fit_logit(s_ref, y_ref)
    p_query = logit_model.predict_proba(scores_query_oriented.reshape(-1, 1))[:, 1]

    # Direct outputs (point estimates)
    pred_binary = (scores_query_oriented >= thr_final).astype(int)
    pred_label = np.where(pred_binary == 1, "responder", "non_responder")

    margin = scores_query_oriented - thr_final  # oriented margin
    confidence = np.maximum(p_query, 1 - p_query)
    ambiguity = 1 - confidence
    margin_prob = 2 * np.abs(p_query - 0.5)

    risk_fp_if_pred_positive = 1 - p_query
    risk_fn_if_pred_negative = p_query
    conditional_error_risk = np.where(pred_binary == 1, risk_fp_if_pred_positive, risk_fn_if_pred_negative)

    # Gray zone point estimates
    if has_gray and np.isfinite(t_low_final) and np.isfinite(t_high_final):
        pred_gray = np.array([_gray_pred(s, t_low_final, t_high_final) for s in scores_query_oriented], dtype=int)
        gray_label = np.array([_gray_label(pg) for pg in pred_gray], dtype=object)
        gray_distance = np.array([_gray_distance(s, t_low_final, t_high_final) for s in scores_query_oriented], dtype=float)
    else:
        pred_gray = np.full(n_query, -999, dtype=int)
        gray_label = np.array(["not_available"] * n_query, dtype=object)
        gray_distance = np.full(n_query, np.nan, dtype=float)

    # Threshold stability (CV-based from fold thresholds)
    # q = P(score >= T)
    threshold_stability_cv = np.array([np.mean(s >= fold_thresholds) for s in scores_query_oriented], dtype=float)
    pred_stability_cv = np.where(pred_binary == 1, threshold_stability_cv, 1 - threshold_stability_cv)

    if has_gray and len(t_low_folds) > 0 and len(t_high_folds) > 0:
        gray_confident_cv = np.array([
            np.mean((s <= t_low_folds) | (s >= t_high_folds)) for s in scores_query_oriented
        ], dtype=float)
    else:
        gray_confident_cv = np.full(n_query, np.nan)

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

        # Need both classes for ROC/logit
        if len(np.unique(yb)) >= 2:
            try:
                model_b = _fit_logit(sb, yb)
                p_vec = model_b.predict_proba(scores_query_oriented.reshape(-1, 1))[:, 1]
            except Exception:
                pass
            try:
                thr_v = _pick_threshold_from_scores(yb, sb, criterion=threshold_criterion)
            except Exception:
                thr_v = np.nan
            if has_gray:
                try:
                    tl_v, th_v = _pick_gray_thresholds_from_scores(
                        yb, sb, target_sens=target_sens, target_spec=target_spec)
                except Exception:
                    tl_v, th_v = np.nan, np.nan
        return p_vec, thr_v, tl_v, th_v

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

    for b, (_p_vec, _thr_v, _tl_v, _th_v) in enumerate(_boot_results):
        p_boot[b, :] = _p_vec
        thr_boot[b] = _thr_v
        if has_gray:
            t_low_boot[b] = _tl_v
            t_high_boot[b] = _th_v

    # Derived bootstrap arrays
    confidence_boot = np.maximum(p_boot, 1 - p_boot)
    ambiguity_boot = 1 - confidence_boot
    margin_prob_boot = 2 * np.abs(p_boot - 0.5)

    risk_fp_boot = 1 - p_boot
    risk_fn_boot = p_boot
    # conditional risk based on FIXED point-estimate prediction
    cond_err_boot = np.where(pred_binary.reshape(1, -1) == 1, risk_fp_boot, risk_fn_boot)

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
        tlv = t_low_boot[valid_gray_mask]
        thv = t_high_boot[valid_gray_mask]
        if tlv.size > 0:
            gray_pred_boot = np.full((tlv.size, n_query), -1, dtype=int)
            for i, s in enumerate(scores_query_oriented):
                gray_pred_boot[:, i] = np.where(s <= tlv, 0, np.where(s >= thv, 1, -1))

            gray_confident_boot = np.mean(gray_pred_boot != -1, axis=0)
            gray_same_label_boot = np.mean(gray_pred_boot == pred_gray.reshape(1, -1), axis=0)
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
            "confidence": float(confidence[i]),
            "ambiguity": float(ambiguity[i]),
            "margin_prob": float(margin_prob[i]),

            # risks
            "risk_FP_if_pred_positive": float(risk_fp_if_pred_positive[i]),
            "risk_FN_if_pred_negative": float(risk_fn_if_pred_negative[i]),
            "conditional_error_risk": float(conditional_error_risk[i]),

            # threshold stability
            "threshold_stability_cv": float(threshold_stability_cv[i]),
            "pred_stability_cv": float(pred_stability_cv[i]),
            "threshold_stability_boot": float(threshold_stability_boot[i]),
            "pred_stability_boot": float(pred_stability_boot[i]),

            # gray zone
            "gray_t_low": float(t_low_final) if np.isfinite(t_low_final) else np.nan,
            "gray_t_high": float(t_high_final) if np.isfinite(t_high_final) else np.nan,
            "pred_gray": int(pred_gray[i]) if pred_gray[i] != -999 else np.nan,
            "pred_gray_label": gray_label[i],
            "gray_distance_to_boundary": float(gray_distance[i]),
            "gray_confident_cv": float(gray_confident_cv[i]),
            "gray_confident_boot": float(gray_confident_boot[i]),
            "gray_label_stability_boot": float(gray_same_label_boot[i]),
        }

        # Bootstrap CIs (probability + derived)
        p_ci = _percentile_ci(p_boot[:, i], ci=ci_level)
        conf_ci = _percentile_ci(confidence_boot[:, i], ci=ci_level)
        amb_ci = _percentile_ci(ambiguity_boot[:, i], ci=ci_level)
        mp_ci = _percentile_ci(margin_prob_boot[:, i], ci=ci_level)
        rfp_ci = _percentile_ci(risk_fp_boot[:, i], ci=ci_level)
        rfn_ci = _percentile_ci(risk_fn_boot[:, i], ci=ci_level)
        cer_ci = _percentile_ci(cond_err_boot[:, i], ci=ci_level)

        # Threshold/margin CI from bootstrap thresholds
        if margin_to_thr_boot.size > 0:
            margin_thr_ci = _percentile_ci(margin_to_thr_boot[:, i], ci=ci_level)
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
            "ci_low": _percentile_ci(thr_boot, ci=ci_level)[0],
            "ci_high": _percentile_ci(thr_boot, ci=ci_level)[1],
        },
    }

    if has_gray:
        bootstrap_summary["gray_t_low_boot_ci"] = {
            "ci_low": _percentile_ci(t_low_boot, ci=ci_level)[0],
            "ci_high": _percentile_ci(t_low_boot, ci=ci_level)[1],
        }
        bootstrap_summary["gray_t_high_boot_ci"] = {
            "ci_low": _percentile_ci(t_high_boot, ci=ci_level)[0],
            "ci_high": _percentile_ci(t_high_boot, ci=ci_level)[1],
        }

    reference_info = {
        "orientation_sign_raw_to_oriented": float(sign),  # +1 or -1
        "score_orientation_text": res_nested_cv_signature.get("score_orientation", "unknown"),
        "reference_n": int(len(y_ref)),
        "reference_n_pos": int(np.sum(y_ref == 1)),
        "reference_n_neg": int(np.sum(y_ref == 0)),
        "threshold_final_median_fold": float(thr_final),
        "gray_zone_available": bool(has_gray),
        "gray_t_low_final_median_fold": float(t_low_final) if np.isfinite(t_low_final) else np.nan,
        "gray_t_high_final_median_fold": float(t_high_final) if np.isfinite(t_high_final) else np.nan,
    }

    return {
        "report": report,
        "bootstrap_summary": bootstrap_summary,
        "reference_info": reference_info,
    }

