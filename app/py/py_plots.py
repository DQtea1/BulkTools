### Contains visualisation functions ### 
# my_Box_Wilcox() : boxplot + wilcox between two conditions
# myROC_AUC() : Draws ROC curve

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from matplotlib.patches import Rectangle
from pathlib import Path

import py.py_models as pyM


def plot_sample_signature_confidence(
    query_scores,
    res_v2: dict,
    query_ids=None,
    confidence_res: dict = None,
    confidence_kwargs: dict = None,
    # Threshold bootstrap density (visualization only)
    n_bootstrap_threshold_plot: int = 1000,
    threshold_criterion: str = "youden",
    bootstrap_seed: int = 123,
    bootstrap_stratified: bool = True,
    # Text / style
    language: str = "fr",   # "fr" or "en"
    figsize=(13, 8),
    dpi=250,
    jitter_strength: float = 0.16,
    density_grid_points: int = 500,
    title: str = None,
    subtitle: str = None,
    point_size_ref: int = 24,
    point_size_query: int = 85,
    alpha_ref: float = 0.60,
    alpha_density_fill: float = 0.10,
    alpha_threshold_density: float = 0.10,
    # Save/show
    save_path: str = None,  # e.g. "figure.png" or "figure.pdf"
    show: bool = False,
    random_state: int = 123,
):
    """
    Publication-ready plot for gene signature score confidence visualization.

    Parameters
    ----------
    query_scores : scalar or list-like
        Raw score(s) to highlight in green.
    res_v2 : dict
        Output of nested_cv_signature_v2(..., return_oof_predictions=True).
    query_ids : list-like, optional
        IDs for query samples. If None, auto-generated.
    confidence_res : dict, optional
        Output of pyM.sample_confidence_report(...). If None, computed internally.
    confidence_kwargs : dict, optional
        kwargs passed to pyM.sample_confidence_report when confidence_res is None.
    save_path : str, optional
        Path to save figure (.png/.pdf/.svg supported by matplotlib).
    """

    # ---------------------------------------------------------------------
    # Helpers
    # ---------------------------------------------------------------------
    def _to_1d(x):
        if np.isscalar(x):
            return np.array([float(x)], dtype=float)
        return np.asarray(x, dtype=float).reshape(-1)

    def _infer_orientation_sign_from_oof(oof_df):
        raw_ref = oof_df["score_raw_original"].to_numpy(dtype=float)
        ori_ref = oof_df["score_used_oriented"].to_numpy(dtype=float)
        diff_plus = np.nanmean(np.abs(ori_ref - raw_ref))
        diff_minus = np.nanmean(np.abs(ori_ref + raw_ref))
        return 1.0 if diff_plus <= diff_minus else -1.0  # oriented = sign * raw

    def _kde_gaussian_manual(x, grid, bandwidth=None):
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        grid = np.asarray(grid, dtype=float)

        if x.size == 0:
            return np.zeros_like(grid)

        if x.size == 1:
            bw = 0.1 if bandwidth is None else float(bandwidth)
            z = (grid - x[0]) / bw
            return np.exp(-0.5 * z*z) / (np.sqrt(2*np.pi)*bw)

        if bandwidth is None:
            std = np.std(x, ddof=1)
            iqr = np.subtract(*np.percentile(x, [75, 25]))
            sigma = std if np.isfinite(std) and std > 0 else np.nan
            if np.isfinite(iqr) and iqr > 0:
                alt = iqr / 1.34
                sigma = np.nanmin([sigma, alt]) if np.isfinite(sigma) else alt
            if not np.isfinite(sigma) or sigma <= 0:
                span = np.max(x) - np.min(x)
                sigma = span / 6 if span > 0 else 1.0
            bandwidth = 0.9 * sigma * (x.size ** (-1/5))
            if not np.isfinite(bandwidth) or bandwidth <= 0:
                bandwidth = max(np.std(x), 1e-3)
                if bandwidth <= 0:
                    bandwidth = 1e-3

        bw = float(bandwidth)
        z = (grid[:, None] - x[None, :]) / bw
        dens = np.mean(np.exp(-0.5*z*z) / (np.sqrt(2*np.pi)*bw), axis=1)
        return dens

    def _scale_density(d, height):
        d = np.asarray(d, dtype=float)
        if d.size == 0 or not np.any(np.isfinite(d)):
            return np.zeros_like(d)
        m = np.nanmax(d)
        return np.zeros_like(d) if (not np.isfinite(m) or m <= 0) else (d / m) * height

    def _pick_threshold(y_true, scores_ref, criterion="youden"):
        fpr, tpr, thresholds = roc_curve(y_true, scores_ref)
        if criterion == "youden":
            idx = int(np.nanargmax(tpr - fpr))
        elif criterion == "closest_topleft":
            d = np.sqrt((fpr - 0)**2 + (tpr - 1)**2)
            idx = int(np.nanargmin(d))
        else:
            raise ValueError("threshold_criterion must be 'youden' or 'closest_topleft'.")
        return float(thresholds[idx])

    def _bootstrap_thresholds(y_ref, s_ref, n_boot=1000, seed=123, stratified=True, criterion="youden"):
        rng = np.random.default_rng(seed)
        y_ref = np.asarray(y_ref, dtype=int)
        s_ref = np.asarray(s_ref, dtype=float)

        idx_all = np.arange(len(y_ref))
        idx_pos = np.where(y_ref == 1)[0]
        idx_neg = np.where(y_ref == 0)[0]

        out = np.full(n_boot, np.nan, dtype=float)

        for b in range(n_boot):
            if stratified and len(idx_pos) > 0 and len(idx_neg) > 0:
                boot_pos = rng.choice(idx_pos, size=len(idx_pos), replace=True)
                boot_neg = rng.choice(idx_neg, size=len(idx_neg), replace=True)
                idx = np.concatenate([boot_pos, boot_neg])
                rng.shuffle(idx)
            else:
                idx = rng.choice(idx_all, size=len(idx_all), replace=True)

            yb = y_ref[idx]
            sb = s_ref[idx]
            if len(np.unique(yb)) < 2:
                continue
            try:
                out[b] = _pick_threshold(yb, sb, criterion=criterion)
            except Exception:
                pass
        return out[np.isfinite(out)]

    def _find_ci_cols(columns, prefix):
        low = [c for c in columns if c.startswith(f"{prefix}_ci_low_")]
        high = [c for c in columns if c.startswith(f"{prefix}_ci_high_")]
        return (low[0] if low else None, high[0] if high else None)

    # ---------------------------------------------------------------------
    # Labels
    # ---------------------------------------------------------------------
    L = {
        "fr": {
            "nr": "Non-répondeurs",
            "r": "Répondeurs",
            "query": "Samples à prédire",
            "score_axis": "Score de signature",
            "prob_axis": "Probabilité calibrée P(réponse=1 | score)",
            "gray_zone": "Zone grise",
            "threshold": "Seuil (binaire)",
            "thr_boot": "Densité thresholds bootstrap",
            "dens_nr": "Densité scores non-répondeurs",
            "dens_r": "Densité scores répondeurs",
            "pred_nr": "Non-répondeur",
            "pred_r": "Répondeur",
            "title": "Confiance de prédiction sur score de signature",
            "subtitle_default": "Scores projetés, densités par classe, seuil, zone grise et proba calibrée",
            "pred": "Prédiction",
            "p": "Proba de Réponse",
            "conf": "Confiance",
            "amb": "Ambiguïté",
            "margin": "Marge (score - seuil)",
            "cer": "Risque d'erreur conditionnel",
            "gray": "Zone grise"
        },
        "en": {
            "nr": "Non-responders",
            "r": "Responders",
            "query": "Query samples",
            "score_axis": "Signature score",
            "prob_axis": "Calibrated probability P(response=1 | score)",
            "gray_zone": "Gray zone",
            "threshold": "Threshold (binary)",
            "thr_boot": "Bootstrap threshold density",
            "dens_nr": "Non-responder score density",
            "dens_r": "Responder score density",
            "pred_nr": "Non-responder",
            "pred_r": "Responder",
            "title": "Prediction confidence on signature score",
            "subtitle_default": "Projected scores, class densities, threshold, gray zone and calibrated probability",
            "pred": "Prediction",
            "p": "Response Proba",
            "conf": "Confidence",
            "amb": "Ambiguity",
            "margin": "Margin (score - threshold)",
            "cer": "Conditional error risk",
            "gray": "Gray zone"
        }
    }["fr" if language.lower().startswith("fr") else "en"]

    # ---------------------------------------------------------------------
    # Validate res_v2
    # ---------------------------------------------------------------------
    if "oof" not in res_v2:
        raise ValueError("res_v2 must contain 'oof'. Run nested_cv_signature_v2(..., return_oof_predictions=True).")

    oof = res_v2["oof"].copy()
    needed = {"y_true", "score_raw_original", "score_used_oriented"}
    if not needed.issubset(oof.columns):
        raise ValueError(f"res_v2['oof'] must contain columns {needed}")

    if "fold_thresholds" not in res_v2 or len(res_v2["fold_thresholds"]) == 0:
        raise ValueError("res_v2 must contain non-empty 'fold_thresholds'.")

    # ---------------------------------------------------------------------
    # Query scores & confidence report
    # ---------------------------------------------------------------------
    query_scores = _to_1d(query_scores)
    n_query = len(query_scores)
    if query_ids is None:
        query_ids = [f"query_{i+1}" for i in range(n_query)]
    if len(query_ids) != n_query:
        raise ValueError("query_ids length must match query_scores length.")

    if confidence_res is None:
        if "pyM.sample_confidence_report" not in globals():
            raise RuntimeError("pyM.sample_confidence_report not found. Pass confidence_res or define the function first.")
        ck = {} if confidence_kwargs is None else dict(confidence_kwargs)
        ck.setdefault("n_bootstrap", 1000)
        ck.setdefault("threshold_criterion", threshold_criterion)
        ck.setdefault("use_gray_zone", True)
        confidence_res = pyM.sample_confidence_report(
            scores=query_scores,
            res_v2=res_v2,
            sample_ids=query_ids,
            **ck
        )

    rep = confidence_res["report"].copy()

    # ---------------------------------------------------------------------
    # Axis orientation / thresholds on raw scale
    # ---------------------------------------------------------------------
    sign = _infer_orientation_sign_from_oof(oof)  # oriented = sign * raw
    thr_oriented = float(np.median(np.asarray(res_v2["fold_thresholds"], dtype=float)))
    thr_raw = sign * thr_oriented

    gray_available = ("gray_zone" in res_v2) and (res_v2["gray_zone"] is not None)
    if gray_available and "fold_gray_thresholds" in res_v2["gray_zone"] and len(res_v2["gray_zone"]["fold_gray_thresholds"]) > 0:
        fg = np.asarray(res_v2["gray_zone"]["fold_gray_thresholds"], dtype=float)
        t_low_raw = sign * float(np.median(fg[:, 0]))
        t_high_raw = sign * float(np.median(fg[:, 1]))
    else:
        t_low_raw = np.nan
        t_high_raw = np.nan
        gray_available = False

    # ---------------------------------------------------------------------
    # Threshold bootstrap density (visual)
    # ---------------------------------------------------------------------
    y_ref = oof["y_true"].to_numpy(dtype=int)
    s_ref_oriented = oof["score_used_oriented"].to_numpy(dtype=float)
    thr_boot_oriented = _bootstrap_thresholds(
        y_ref=y_ref,
        s_ref=s_ref_oriented,
        n_boot=int(n_bootstrap_threshold_plot),
        seed=bootstrap_seed,
        stratified=bootstrap_stratified,
        criterion=threshold_criterion
    )
    thr_boot_raw = sign * thr_boot_oriented if len(thr_boot_oriented) else np.array([])

    # ---------------------------------------------------------------------
    # Main data (raw scores)
    # ---------------------------------------------------------------------
    scores_raw = oof["score_raw_original"].to_numpy(dtype=float)
    y = oof["y_true"].to_numpy(dtype=int)

    s_non = scores_raw[y == 0]
    s_resp = scores_raw[y == 1]

    x_candidates = np.concatenate([
        scores_raw,
        query_scores,
        np.array([thr_raw]),
        np.array([t_low_raw, t_high_raw]) if gray_available else np.array([])
    ])
    if thr_boot_raw.size:
        x_candidates = np.concatenate([x_candidates, thr_boot_raw])

    x_candidates = x_candidates[np.isfinite(x_candidates)]
    x_min, x_max = np.min(x_candidates), np.max(x_candidates)
    span = x_max - x_min if x_max > x_min else 1.0
    pad = 0.08 * span
    xlim = (x_min - pad, x_max + pad)

    grid = np.linspace(*xlim, density_grid_points)
    d_non = _kde_gaussian_manual(s_non, grid)
    d_resp = _kde_gaussian_manual(s_resp, grid)
    d_thr = _kde_gaussian_manual(thr_boot_raw, grid) if thr_boot_raw.size > 1 else np.zeros_like(grid)

    # Vertical layout in top axis
    rect_y0, rect_h = 0.00, 1.00
    thr_base, thr_h = 1.10, 0.55
    cls_base, cls_h = 1.85, 1.10

    y_thr_curve = thr_base + _scale_density(d_thr, thr_h)
    y_non_curve = cls_base + _scale_density(d_non, cls_h)
    y_resp_curve = cls_base + _scale_density(d_resp, cls_h)

    # Jitter points inside rectangle
    rng = np.random.default_rng(random_state)
    y_pts = rect_y0 + rect_h * 0.5 + rng.uniform(
        -jitter_strength * rect_h, jitter_strength * rect_h, size=len(scores_raw)
    )

    # ---------------------------------------------------------------------
    # Figure
    # ---------------------------------------------------------------------
    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = fig.add_gridspec(2, 1, height_ratios=[4.7, 1.35], hspace=0.28)
    ax = fig.add_subplot(gs[0, 0])
    axp = fig.add_subplot(gs[1, 0])

    # Top panel background rectangle
    ax.add_patch(Rectangle(
        (xlim[0], rect_y0), xlim[1] - xlim[0], rect_h,
        facecolor="0.92", edgecolor="0.30", linewidth=0.9, alpha=0.45
    ))

    # Gray zone (shaded)
    if gray_available and np.isfinite(t_low_raw) and np.isfinite(t_high_raw):
        x0, x1 = min(t_low_raw, t_high_raw), max(t_low_raw, t_high_raw)
        ax.add_patch(Rectangle(
            (x0, rect_y0), x1 - x0, rect_h,
            facecolor="gold", edgecolor=None, alpha=0.20
        ))
        ax.text((x0 + x1) / 2, rect_y0 + rect_h + 0.03, L["gray_zone"],
                ha="center", va="bottom", fontsize=10, color="darkgoldenrod")

    # Reference jitter points
    mask_non = (y == 0)
    mask_resp = (y == 1)
    ax.scatter(scores_raw[mask_non], y_pts[mask_non], s=point_size_ref,
               color="red", alpha=alpha_ref, edgecolor="none", label=L["nr"])
    ax.scatter(scores_raw[mask_resp], y_pts[mask_resp], s=point_size_ref,
               color="blue", alpha=alpha_ref, edgecolor="none", label=L["r"])

    # Query points (green)
    # slight vertical offsets if multiple points
    if n_query == 1:
        yq = np.array([rect_y0 + rect_h * 0.5])
    else:
        offsets = np.linspace(-0.10, 0.10, n_query)
        yq = rect_y0 + rect_h * 0.5 + offsets
    ax.scatter(query_scores, yq, s=point_size_query, color="green",
               edgecolor="black", linewidth=0.8, zorder=6, label=L["query"])

    # Label query IDs (small)
    for xi, yi, qid in zip(query_scores, yq, query_ids):
        ax.text(xi, yi + 0.11, str(qid), color="green", fontsize=8.5,
                ha="center", va="bottom")

    # Class densities
    ax.plot(grid, y_non_curve, color="red", lw=2.0, label=L["dens_nr"])
    ax.fill_between(grid, cls_base, y_non_curve, color="red", alpha=alpha_density_fill)

    ax.plot(grid, y_resp_curve, color="blue", lw=2.0, label=L["dens_r"])
    ax.fill_between(grid, cls_base, y_resp_curve, color="blue", alpha=alpha_density_fill)

    # Threshold bootstrap density
    if thr_boot_raw.size > 1:
        ax.plot(grid, y_thr_curve, color="0.35", lw=2.0, ls="--", label=L["thr_boot"])
        ax.fill_between(grid, thr_base, y_thr_curve, color="0.35", alpha=alpha_threshold_density)

    # Vertical lines: threshold + gray boundaries + query
    ax.axvline(thr_raw, color="black", lw=2.0, ls="-", label=L["threshold"])
    if gray_available and np.isfinite(t_low_raw) and np.isfinite(t_high_raw):
        ax.axvline(t_low_raw, color="darkgoldenrod", lw=1.6, ls=":")
        ax.axvline(t_high_raw, color="darkgoldenrod", lw=1.6, ls=":")

    for qs in query_scores:
        ax.axvline(qs, color="green", lw=1.2, ls="-.", alpha=0.65)

    # Cosmetics top panel
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.12, cls_base + cls_h + 0.30)
    ax.set_yticks([])
    ax.set_xlabel(L["score_axis"])
    ax.grid(axis="x", alpha=0.15)

    if title is None:
        title = L["title"]
    ax.set_title(title, fontsize=13, weight="bold", pad=12)

    if subtitle is None:
        subtitle = L["subtitle_default"]
    ax.text(0.0, 1.02, subtitle, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=9.5, color="0.35")

    # Compact legend
    ax.legend(loc="upper left", fontsize=8.8, frameon=True, ncol=2)

    # ---------------------------------------------------------------------
    # Bottom panel: calibrated probability + CI for query samples
    # ---------------------------------------------------------------------
    axp.set_xlim(0, 1)
    axp.set_ylim(0, 1)
    axp.set_yticks([])
    axp.set_xlabel(L["prob_axis"])

    # Background bar + midline
    axp.add_patch(Rectangle((0, 0.30), 1.0, 0.40, facecolor="0.92", edgecolor="0.30", alpha=0.45))
    axp.axvline(0.5, color="black", lw=1.0, ls="--", alpha=0.7)

    # plot each query sample probability + CI
    y_levels = np.linspace(0.62, 0.38, n_query) if n_query > 1 else np.array([0.50])

    p_rows = []
    for i, qid in enumerate(query_ids):
        r = rep.loc[qid]
        p_val = float(r["p_calibrated_logit"] if "p_calibrated_logit" in r.index else r["p"])
        low_col, high_col = _find_ci_cols(r.index, "p")
        p_low = float(r[low_col]) if low_col else np.nan
        p_high = float(r[high_col]) if high_col else np.nan

        y0 = y_levels[i]

        if np.isfinite(p_low) and np.isfinite(p_high):
            axp.plot([p_low, p_high], [y0, y0], color="green", lw=5, alpha=0.25, solid_capstyle="round")
            axp.plot([p_low, p_high], [y0, y0], color="green", lw=2.0)

        axp.scatter([p_val], [y0], s=65, color="green", edgecolor="black", linewidth=0.7, zorder=5)

        label_txt = f"{qid}: p={p_val:.3f}"
        if np.isfinite(p_low) and np.isfinite(p_high):
            label_txt += f" [{p_low:.3f}, {p_high:.3f}]"
        axp.text(0.01, y0 + 0.06, label_txt, ha="left", va="bottom", fontsize=8.5, transform=axp.transData)

        p_rows.append((qid, p_val, p_low, p_high))

    # Left/right class labels
    axp.text(0.00, 0.05, L["pred_nr"], color="red", ha="left", va="bottom", fontsize=9)
    axp.text(1.00, 0.05, L["pred_r"], color="blue", ha="right", va="bottom", fontsize=9)

    # Remove spines
    for sp in ["top", "right", "left"]:
        axp.spines[sp].set_visible(False)

    # ---------------------------------------------------------------------
    # Side annotation table (first query sample summary only, compact)
    # ---------------------------------------------------------------------
    r0 = rep.loc[query_ids[0]]
    low_col, high_col = _find_ci_cols(r0.index, "p")
    p_low0 = float(r0[low_col]) if low_col else np.nan
    p_high0 = float(r0[high_col]) if high_col else np.nan

    pred_lbl = str(r0.get("pred_label", "NA"))
    conf0 = float(r0.get("confidence", np.nan))
    amb0 = float(r0.get("ambiguity", np.nan))
    margin0 = float(r0.get("margin", np.nan))
    cer0 = float(r0.get("conditional_error_risk", np.nan))
    gz_lbl0 = str(r0.get("pred_gray_label", "not_available"))

    info = (f"{L['pred']}: {pred_lbl}\n")
    if gray_available:
        info += f"{L['gray']}: {gz_lbl0}\n"
    info += f"{L['p']}: {float(r0.get('p', r0.get('p_calibrated_logit', np.nan))):.3f}"
    
    if np.isfinite(p_low0) and np.isfinite(p_high0):
        info += f" [{p_low0:.3f}, {p_high0:.3f}]"
    info += (
        # f"\n{L['conf']}: {conf0:.3f}"
        # f"\n{L['amb']}: {amb0:.3f}"
        f"\n{L['margin']}: {margin0:.3f}"
        # f"\n{L['cer']}: {cer0:.3f}"
    )

    ax.text(
        0.995, 0.98, info,
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9.2,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.6", alpha=0.92)
    )

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", dpi=dpi)

    if show:
        plt.show()

    return {
        "fig": fig,
        "axes": (ax, axp),
        "confidence_report": rep,
        "threshold_boot_raw": thr_boot_raw,
        "threshold_raw": thr_raw,
        "gray_t_low_raw": t_low_raw,
        "gray_t_high_raw": t_high_raw,
        "save_path" : save_path
    }


def my_Box_Wilcox(df, responders, non_responders, score_col="score",
                  out_dir=None, return_mode="path", random_seed=42,
                  show_points=True):
    """
    Compare scores between responders and non-responders (Mann-Whitney U),
    and return either a saved plot path (recommended for R/Shiny) or a matplotlib figure.

    Parameters
    ----------
    df : pandas.DataFrame or pandas.Series
        If DataFrame, `score_col` is used. If Series, values are used directly.
        Index must contain sample IDs.
    responders : list-like
        IDs of responder samples.
    non_responders : list-like
        IDs of non-responder samples.
    score_col : str
        Column name if df is a DataFrame.
    out_dir : str or None
        Output PNG path if return_mode='path'. If None, a temporary file is created.
    return_mode : {'path', 'figure'}
        'path' -> saves PNG and returns path + stats
        'figure' -> returns fig object + stats
    random_seed : int
        Seed for jittered points.
    show_points : bool
        Whether to overlay jittered individual points.

    Returns
    -------
    dict
        Keys:
          - p_value
          - u_stat
          - n_responders
          - n_non_responders
          - plot_path (if return_mode='path')
          - figure (if return_mode='figure')
    """
    import os
    import tempfile
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import mannwhitneyu

    # --- get a 1D scores Series no matter what df is ---
    if isinstance(df, pd.DataFrame):
        if score_col not in df.columns:
            raise ValueError(f"Column '{score_col}' not found in DataFrame columns: {list(df.columns)}")
        scores = df[score_col].copy()
    else:
        scores = pd.Series(df).copy()  # supports Series-like input too

    if scores.index is None:
        raise ValueError("Scores must have an index containing sample IDs.")

    # Make IDs comparable (str vs int mismatch)
    scores.index = scores.index.astype(str)
    responders_ids = [str(x) for x in responders]
    non_responders_ids = [str(x) for x in non_responders]

    # Check overlap
    overlap = set(responders_ids) & set(non_responders_ids)
    if overlap:
        raise ValueError(f"{len(overlap)} IDs are in BOTH lists (example: {list(overlap)[:5]})")

    # Get the scores
    r  = scores.reindex(responders_ids).astype(float).dropna().to_numpy()
    nr = scores.reindex(non_responders_ids).astype(float).dropna().to_numpy()

    # Basic sanity checks
    if len(r) == 0 or len(nr) == 0:
        raise ValueError(
            f"Empty group after reindex/dropna: responders={len(r)}, non_responders={len(nr)}. "
            "Check IDs and score index."
        )

    # Wilcoxon-type test for independent groups (Mann–Whitney U)
    u_stat, pval = mannwhitneyu(r, nr, alternative="two-sided")

    # Create plot
    rng = np.random.default_rng(random_seed)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.boxplot([r, nr], tick_labels=["Responders", "Non-responders"], showfliers=False)
    ax.set_ylabel("Score")
    ax.set_title("Score by response status")

    if show_points:
        ax.plot(rng.normal(1, 0.04, len(r)),  r,  "o", ms=3, alpha=0.5)
        ax.plot(rng.normal(2, 0.04, len(nr)), nr, "o", ms=3, alpha=0.5)

    ax.text(
        0.5, 0.9,
        f"Mann–Whitney p = {pval:.3g}",
        transform=ax.transAxes,
        ha="center"
    )

    fig.tight_layout()

    result = {
        "p_value": float(pval),
        "u_stat": float(u_stat),
        "n_responders": int(len(r)),
        "n_non_responders": int(len(nr)),
    }

    if return_mode == "figure":
        result["figure"] = fig
        return result

    # Default: save PNG and return path (best for R/Shiny)
    if out_dir is None:
        out_dir = tempfile.mkdtemp()

    out_dir = Path(out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    save_path = out_dir / "wilcox_boxplot.png"
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    result["save_path"] = str(save_path)
    return {"plot": fig, "save_path": str(save_path), **result}


def myROC_AUC(df, responders, non_responders,
              threshold_conf_mat=0.5,
              score_col="score",
              out_dir=None):

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from pathlib import Path
    from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, ConfusionMatrixDisplay

    # --- get scores as a 1D Series ---
    scores = df[score_col] if isinstance(df, pd.DataFrame) else pd.Series(df)
    scores = scores.copy()
    scores.index = scores.index.astype(str)

    responders_ids = pd.Index([str(x) for x in responders])
    non_responders_ids = pd.Index([str(x) for x in non_responders])

    overlap = responders_ids.intersection(non_responders_ids)
    if len(overlap) > 0:
        raise ValueError(f"{len(overlap)} IDs are in BOTH lists (example: {list(overlap[:5])})")

    # --- build labels y aligned to scores.index ---
    y = pd.Series(np.nan, index=scores.index, dtype="float64")
    y.loc[scores.index.intersection(responders_ids)] = 1
    y.loc[scores.index.intersection(non_responders_ids)] = 0

    # Keep only patients with both a score and a label
    mask = y.notna() & scores.notna()
    y_true = y.loc[mask].astype(int).to_numpy()
    y_score = scores.loc[mask].astype(float).to_numpy()

    print(f"Using N={len(y_true)} (responders={y_true.sum()}, non-responders={(y_true==0).sum()})")

    # --- Confusion matrix (same length arrays!) ---
    # If "higher score => responder", use >= threshold for class 1
    y_pred = (y_score >= threshold_conf_mat).astype(int)

    conf_mat = confusion_matrix(y_true, y_pred, labels=[0, 1])
    fig_cm, ax_cm = plt.subplots(figsize=(4.5, 4))
    disp = ConfusionMatrixDisplay(confusion_matrix=conf_mat,
                                  display_labels=["Non-Responder", "Responder"])
    disp.plot(ax=ax_cm, values_format="d")
    fig_cm.tight_layout()

    # --- ROC + AUC ---
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc = roc_auc_score(y_true, y_score)

    fig_roc, ax = plt.subplots(figsize=(5, 4))
    ax.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC curve")
    ax.legend(loc="lower right")
    fig_roc.tight_layout()

    save_path = None
    save_cm_path = None
    if out_dir is not None:
        out_dir = Path(out_dir).expanduser()
        out_dir.mkdir(parents=True, exist_ok=True)
        save_path = out_dir / "ROC_AUC.png"
        save_cm_path = out_dir / "confusion_matrix.png"
        fig_roc.savefig(save_path, dpi=150, bbox_inches="tight")
        fig_cm.savefig(save_cm_path, dpi=150, bbox_inches="tight")

    return {
        "plot": fig_roc,
        "confmat_plot": fig_cm,
        "auc": float(auc),
        "roc_save_path": str(save_path) if save_path else None,
        "cm_save_path": str(save_cm_path) if save_cm_path else None
    }
    
def myROC_AUC_v2(
    df,
    responders,
    non_responders,
    score_col="score",
    threshold_method="youden",   # "youden", "closest_topleft", "f1", "fixed"
    threshold_conf_mat=0.5,      # used only if threshold_method="fixed"
    plot_criterion=False,
    figsize_roc=(5, 4),
    figsize_cm=(5, 4),
):
    """
    Compute ROC-AUC, automatically select a threshold from precomputed scores,
    and plot ROC curve + confusion matrix.

    Parameters
    ----------
    df : pd.DataFrame or pd.Series
        Either:
        - a DataFrame containing a score column
        - or a Series of scores
        Index must contain patient IDs.
    responders : list
        List of responder patient IDs.
    non_responders : list
        List of non-responder patient IDs.
    score_col : str
        Name of the score column if df is a DataFrame.
    threshold_method : str
        Threshold selection method:
        - "youden"          : maximize TPR - FPR
        - "closest_topleft" : minimize distance to (0,1) on ROC
        - "f1"              : maximize F1 score
        - "fixed"           : use threshold_conf_mat directly
    threshold_conf_mat : float
        Used only if threshold_method == "fixed".
    plot_criterion : bool
        If True, plot the threshold selection criterion.
    figsize_roc : tuple
        Figure size for ROC plot.
    figsize_cm : tuple
        Figure size for confusion matrix plot.

    Returns
    -------
    dict
        Contains y_true, y_score, auc, selected_threshold, confusion_matrix, etc.
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    from sklearn.metrics import (
        roc_auc_score,
        roc_curve,
        confusion_matrix,
        ConfusionMatrixDisplay,
        precision_recall_curve,
        accuracy_score,
        precision_score,
        recall_score,
        f1_score,
    )

    # ------------------------------------------------------------
    # Helper to choose threshold
    # ------------------------------------------------------------
    def _select_threshold(y_true, y_score, method="youden", fixed_threshold=0.5):
        y_true = np.asarray(y_true).astype(int)
        y_score = np.asarray(y_score).astype(float)

        if method == "fixed":
            return float(fixed_threshold), {"method": "fixed"}

        if method in ("youden", "closest_topleft"):
            fpr, tpr, thresholds = roc_curve(y_true, y_score)

            # sklearn may include inf as first threshold
            valid = np.isfinite(thresholds)
            fpr = fpr[valid]
            tpr = tpr[valid]
            thresholds = thresholds[valid]

            if len(thresholds) == 0:
                raise ValueError("No valid threshold available from roc_curve().")

            if method == "youden":
                crit = tpr - fpr
                best_idx = np.argmax(crit)
                info = {
                    "method": "youden",
                    "thresholds": thresholds,
                    "criterion_values": crit,
                    "best_idx": int(best_idx),
                    "fpr": fpr,
                    "tpr": tpr,
                }
                return float(thresholds[best_idx]), info

            if method == "closest_topleft":
                crit = np.sqrt((fpr - 0.0) ** 2 + (tpr - 1.0) ** 2)
                best_idx = np.argmin(crit)
                info = {
                    "method": "closest_topleft",
                    "thresholds": thresholds,
                    "criterion_values": crit,
                    "best_idx": int(best_idx),
                    "fpr": fpr,
                    "tpr": tpr,
                }
                return float(thresholds[best_idx]), info

        if method == "f1":
            precision, recall, thresholds = precision_recall_curve(y_true, y_score)

            # thresholds has length len(precision)-1
            denom = precision[:-1] + recall[:-1]
            f1_vals = np.where(denom == 0, 0, 2 * precision[:-1] * recall[:-1] / denom)

            valid = np.isfinite(thresholds)
            thresholds = thresholds[valid]
            f1_vals = f1_vals[valid]

            if len(thresholds) == 0:
                raise ValueError("No valid threshold available from precision_recall_curve().")

            best_idx = np.argmax(f1_vals)
            info = {
                "method": "f1",
                "thresholds": thresholds,
                "criterion_values": f1_vals,
                "best_idx": int(best_idx),
                "precision": precision,
                "recall": recall,
            }
            return float(thresholds[best_idx]), info

        raise ValueError(
            "threshold_method must be one of: "
            "'youden', 'closest_topleft', 'f1', 'fixed'"
        )

    # ------------------------------------------------------------
    # Extract scores
    # ------------------------------------------------------------
    scores = df[score_col] if isinstance(df, pd.DataFrame) else df
    scores = scores.copy()

    if not isinstance(scores, pd.Series):
        scores = pd.Series(scores)

    scores.index = scores.index.astype(str)

    responders_ids = [str(x) for x in responders]
    non_responders_ids = [str(x) for x in non_responders]

    overlap = set(responders_ids) & set(non_responders_ids)
    if overlap:
        raise ValueError(
            f"{len(overlap)} IDs are in BOTH lists "
            f"(example: {list(overlap)[:5]})"
        )

    # ------------------------------------------------------------
    # Build labels
    # ------------------------------------------------------------
    y = pd.Series(index=scores.index, dtype="float64")
    y.loc[list(set(responders_ids) & set(scores.index))] = 1
    y.loc[list(set(non_responders_ids) & set(scores.index))] = 0

    # Keep only labeled samples with non-missing scores
    mask = y.notna() & scores.notna()
    y_true = y.loc[mask].astype(int).to_numpy()
    y_score = scores.loc[mask].astype(float).to_numpy()
    used_ids = scores.loc[mask].index.to_list()

    if len(y_true) == 0:
        raise ValueError("No usable patients after matching IDs and removing missing scores.")

    if len(np.unique(y_true)) < 2:
        raise ValueError("Need at least two classes to compute ROC-AUC.")

    print(
        f"Using N={len(y_true)} patients "
        f"(responders={int((y_true == 1).sum())}, "
        f"non-responders={int((y_true == 0).sum())})"
    )

    # ------------------------------------------------------------
    # ROC / AUC
    # ------------------------------------------------------------
    fpr, tpr, roc_thresholds = roc_curve(y_true, y_score)
    auc = roc_auc_score(y_true, y_score)

    # ------------------------------------------------------------
    # Threshold selection
    # ------------------------------------------------------------
    selected_threshold, threshold_info = _select_threshold(
        y_true=y_true,
        y_score=y_score,
        method=threshold_method,
        fixed_threshold=threshold_conf_mat,
    )

    y_pred = (y_score >= selected_threshold).astype(int)

    # ------------------------------------------------------------
    # Confusion matrix + metrics
    # ------------------------------------------------------------
    conf_mat = confusion_matrix(y_true, y_pred, labels=[0, 1])

    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)

    # ------------------------------------------------------------
    # ROC plot
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize_roc)
    ax.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC curve")
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------------------
    # Confusion matrix plot
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize_cm)
    disp = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat,
        display_labels=["Non-Responder", "Responder"]
    )
    disp.plot(ax=ax, colorbar=False)
    ax.set_title(
        f"Confusion matrix\n"
        f"threshold_method = {threshold_method}, threshold = {selected_threshold:.3f}"
    )
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------------------
    # Optional criterion plot
    # ------------------------------------------------------------
    if plot_criterion and threshold_method != "fixed":
        thr = threshold_info["thresholds"]
        crit = threshold_info["criterion_values"]
        best_idx = threshold_info["best_idx"]

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot(thr, crit)
        ax.axvline(thr[best_idx], linestyle="--", label=f"Best threshold = {thr[best_idx]:.3f}")
        ax.set_xlabel("Threshold")
        if threshold_method == "youden":
            ax.set_ylabel("Youden index (TPR - FPR)")
            ax.set_title("Threshold selection by Youden index")
        elif threshold_method == "closest_topleft":
            ax.set_ylabel("Distance to top-left")
            ax.set_title("Threshold selection by closest top-left")
        elif threshold_method == "f1":
            ax.set_ylabel("F1 score")
            ax.set_title("Threshold selection by F1")
        ax.legend()
        plt.tight_layout()
        plt.show()

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------
    print(f"ROC-AUC = {auc:.4f}")
    print(f"Selected threshold method = {threshold_method}")
    print(f"Selected threshold = {selected_threshold:.6f}")
    print("Confusion matrix:")
    print(conf_mat)
    print(f"ACC={acc:.4f} | PREC={prec:.4f} | REC={rec:.4f} | F1={f1:.4f}")

    return {
        "ids": used_ids,
        "y_true": y_true,
        "y_score": y_score,
        "y_pred": y_pred,
        "auc": auc,
        "fpr": fpr,
        "tpr": tpr,
        "roc_thresholds": roc_thresholds,
        "selected_threshold": selected_threshold,
        "threshold_method": threshold_method,
        "threshold_info": threshold_info,
        "confusion_matrix": conf_mat,
        "accuracy": acc,
        "precision": prec,
        "recall": rec,
        "f1": f1,
    }