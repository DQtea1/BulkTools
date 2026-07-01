### Independent Component Analysis (ICA) of bulk transcriptomic data ###
# Two entry points, called from R (R_analysis/ICA.R via reticulate):
#   compute_mstd()          : Most Stable Transcriptomic Dimension (part 1)
#   run_stabilized_ica()    : stabilized ICA -> S (metagenes) + A (activities) (part 2)
# Plotting helpers (part 2), each saves a PNG and returns its path:
#   scatter_with_marginals(), plot_A_heatmap(), plot_source_distributions(),
#   plot_stability_index(), plot_component_correlation(), plot_clinical_association()
#
# All heavy numeric objects are exchanged with R as plain numpy matrices plus the
# corresponding name vectors, so that reticulate round-trips do not depend on how
# pandas indexes are preserved.

import os

import numpy as np
import pandas as pd
import matplotlib

# Headless rendering (the app runs without an X server); safe to call before pyplot.
try:
    matplotlib.use("Agg")
except Exception:
    pass
import matplotlib.pyplot as plt
from scipy import stats


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _as_matrix(x):
    """Coerce an object coming from R (matrix / data.frame) into a float ndarray."""
    if isinstance(x, pd.DataFrame):
        return x.to_numpy(dtype=float)
    return np.asarray(x, dtype=float)


def _expr_fit_matrix(expr, gene_names, sample_names):
    """Return X of shape (n_samples, n_genes) for sica.

    `expr` is the bulk matrix genes x samples (as passed from R). sica expects
    observations (samples) in rows and features (genes) in columns.
    """
    mat = _as_matrix(expr)  # genes x samples
    if mat.shape[0] != len(gene_names) or mat.shape[1] != len(sample_names):
        # Be forgiving if R passed it already oriented samples x genes.
        if mat.shape[0] == len(sample_names) and mat.shape[1] == len(gene_names):
            mat = mat.T
        else:
            raise ValueError(
                "Expression matrix shape %s does not match %d genes x %d samples."
                % (mat.shape, len(gene_names), len(sample_names))
            )
    return mat.T  # samples x genes


def _map_algorithm(algorithm):
    alg_map = {
        "parallel": "fastica_par",
        "deflation": "fastica_def",
        "fastica_par": "fastica_par",
        "fastica_def": "fastica_def",
    }
    return alg_map.get(str(algorithm), str(algorithm))


def _pval_to_stars(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def _empty_panel(save_path, message):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.set_axis_off()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


# --------------------------------------------------------------------------- #
# Part 1 : Most Stable Transcriptomic Dimension
# --------------------------------------------------------------------------- #
def compute_mstd(expr, gene_names, sample_names, m, M, step, n_runs,
                 out_dir, max_iter=2000, n_jobs=-1):
    """Run the MSTD diagnostic (stability index vs number of components).

    Returns {"save_path": <png path>}.
    """
    from sica.mstd import MSTD

    X = _expr_fit_matrix(expr, gene_names, sample_names)
    m, M, step, n_runs = int(m), int(M), int(step), int(n_runs)
    os.makedirs(out_dir, exist_ok=True)

    # MSTD draws into the current figure (it creates its own axes when ax=None).
    MSTD(X, m, M, step, n_runs, max_iter=int(max_iter), n_jobs=int(n_jobs))
    fig = plt.gcf()
    fig.set_size_inches(13, 6)
    save_path = os.path.join(out_dir, "MSTD_plot.png")
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return {"save_path": save_path}


# --------------------------------------------------------------------------- #
# Part 2 : stabilized ICA
# --------------------------------------------------------------------------- #
def run_stabilized_ica(expr, gene_names, sample_names, n_components, n_runs,
                       algorithm="fastica_par", max_iter=2000, n_jobs=-1):
    """Fit a StabilizedICA and return S (metagenes), A (activities) and stability.

    Returned matrices are plain lists/arrays plus names so R can rebuild
    data.frames deterministically:
      - S : (n_components, n_genes)   rows = components, cols = genes
      - A : (n_samples, n_components) rows = samples,    cols = components
    """
    from sica.base import StabilizedICA

    X = _expr_fit_matrix(expr, gene_names, sample_names)
    n_components = int(n_components)
    n_runs = int(n_runs)
    algorithm = _map_algorithm(algorithm)
    gene_names = list(gene_names)
    sample_names = list(sample_names)

    # sica 2.x takes n_runs in the constructor; older releases take it in fit().
    try:
        sICA = StabilizedICA(
            n_components=n_components, n_runs=n_runs, algorithm=algorithm,
            max_iter=int(max_iter), n_jobs=int(n_jobs),
        )
        sICA.fit(X)
    except TypeError:
        sICA = StabilizedICA(
            n_components=n_components, max_iter=int(max_iter), n_jobs=int(n_jobs),
        )
        sICA.fit(X, n_runs=n_runs, fun="logcosh", algorithm=algorithm)

    S_mat = np.asarray(sICA.S_, dtype=float)          # components x genes
    A_mat = getattr(sICA, "A_", None)
    if A_mat is None:
        # A = X . pinv(S) : (samples x genes) . (genes x components)
        A_mat = X @ np.linalg.pinv(S_mat)
    A_mat = np.asarray(A_mat, dtype=float)            # samples x components

    stab = None
    for attr in ("stability_indexes_", "stability_index_", "stability_"):
        if hasattr(sICA, attr):
            stab = np.asarray(getattr(sICA, attr), dtype=float).ravel()
            break
    if stab is None or stab.size != n_components:
        stab = np.full(n_components, np.nan)

    comp_names = ["IC%d" % (i + 1) for i in range(n_components)]
    # Return numpy arrays so reticulate hands R plain matrices / vectors.
    return {
        "S": S_mat,                     # components x genes
        "A": A_mat,                     # samples x components
        "stability": stab,              # length n_components
        "comp_names": comp_names,
        "gene_names": gene_names,
        "sample_names": sample_names,
    }


# --------------------------------------------------------------------------- #
# Part 2 : plots
# --------------------------------------------------------------------------- #
def scatter_with_marginals(A, comp_names, sample_names, comp_x, comp_y, out_dir,
                           clinic=None, color_col=None):
    """Scatter of two components with marginal histograms.

    Colours the points by `color_col` of `clinic` (aligned on sample IDs) when
    both are provided; otherwise a single colour is used.
    """
    os.makedirs(out_dir, exist_ok=True)
    A = _as_matrix(A)  # samples x components
    comp_names = list(comp_names)
    sample_names = list(sample_names)
    save_path = os.path.join(out_dir, "scatter_%s_vs_%s.png" % (str(comp_x), str(comp_y)))

    if comp_x not in comp_names or comp_y not in comp_names:
        return _empty_panel(save_path, "Component(s) not found:\n%s / %s" % (comp_x, comp_y))

    xi = comp_names.index(comp_x)
    yi = comp_names.index(comp_y)
    x = A[:, xi]
    y = A[:, yi]

    # Optional colouring from the clinical table.
    color_values = None
    if clinic is not None and color_col is not None and str(color_col) != "":
        color_values = _align_clinic_column(clinic, sample_names, color_col)

    fig = plt.figure(figsize=(9, 8))
    gs = fig.add_gridspec(4, 4, hspace=0.05, wspace=0.05)
    ax = fig.add_subplot(gs[1:4, 0:3])
    ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax)
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax)

    if color_values is None:
        ax.scatter(x, y, s=35, alpha=0.75, color="#1F30C9", edgecolor="white", linewidth=0.4)
        ax_top.hist(x, bins=30, color="#1F30C9", alpha=0.7)
        ax_right.hist(y, bins=30, orientation="horizontal", color="#1F30C9", alpha=0.7)
    else:
        vals = color_values
        numeric = pd.api.types.is_numeric_dtype(vals) and vals.nunique(dropna=True) > 6
        if numeric:
            sc = ax.scatter(x, y, s=35, alpha=0.8, c=vals.astype(float).values,
                            cmap="viridis", edgecolor="white", linewidth=0.4)
            cbar = fig.colorbar(sc, ax=ax_right, shrink=0.7)
            cbar.set_label(str(color_col))
            ax_top.hist(x, bins=30, color="grey", alpha=0.6)
            ax_right.hist(y, bins=30, orientation="horizontal", color="grey", alpha=0.6)
        else:
            cats = vals.astype("object").where(vals.notna(), "NA").astype(str)
            levels = sorted(cats.unique())
            cmap = plt.get_cmap("tab10" if len(levels) <= 10 else "tab20")
            for k, lev in enumerate(levels):
                mask = (cats.values == lev)
                col = cmap(k % cmap.N)
                ax.scatter(x[mask], y[mask], s=35, alpha=0.8, color=col,
                           edgecolor="white", linewidth=0.4, label=lev)
                ax_top.hist(x[mask], bins=20, color=col, alpha=0.5)
                ax_right.hist(y[mask], bins=20, orientation="horizontal", color=col, alpha=0.5)
            ax.legend(title=str(color_col), fontsize=8, title_fontsize=9,
                      loc="best", frameon=True)

    ax.set_xlabel(str(comp_x))
    ax.set_ylabel(str(comp_y))
    ax.axhline(0, color="grey", lw=0.6, ls="--")
    ax.axvline(0, color="grey", lw=0.6, ls="--")
    plt.setp(ax_top.get_xticklabels(), visible=False)
    plt.setp(ax_right.get_yticklabels(), visible=False)
    ax_top.set_yticks([]); ax_right.set_xticks([])
    fig.suptitle("Sample activities: %s vs %s" % (comp_x, comp_y), fontsize=13)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_A_heatmap(A, comp_names, sample_names, out_dir):
    """Heatmap of the mixing matrix A (components x samples)."""
    os.makedirs(out_dir, exist_ok=True)
    A = _as_matrix(A)  # samples x components
    comp_names = list(comp_names)
    sample_names = list(sample_names)
    save_path = os.path.join(out_dir, "A_heatmap.png")

    M = A.T  # components x samples
    n_comp, n_samp = M.shape
    vmax = np.nanmax(np.abs(M)) if M.size else 1.0
    fig, ax = plt.subplots(figsize=(max(6, n_samp * 0.18 + 3), max(4, n_comp * 0.4 + 2)))
    im = ax.imshow(M, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(n_comp)); ax.set_yticklabels(comp_names, fontsize=8)
    if n_samp <= 60:
        ax.set_xticks(range(n_samp)); ax.set_xticklabels(sample_names, rotation=90, fontsize=6)
    else:
        ax.set_xticks([]); ax.set_xlabel("%d samples" % n_samp)
    ax.set_title("Mixing matrix A (component activity per sample)", fontsize=12)
    fig.colorbar(im, ax=ax, shrink=0.7, label="activity")
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_source_distributions(S, comp_names, gene_names, out_dir, top_n=20):
    """For each metagene (source), show its weight distribution and top genes."""
    os.makedirs(out_dir, exist_ok=True)
    S = _as_matrix(S)  # components x genes
    comp_names = list(comp_names)
    gene_names = list(gene_names)
    top_n = int(top_n)
    save_path = os.path.join(out_dir, "source_distributions.png")

    n_comp = S.shape[0]
    ncol = 2 if n_comp > 1 else 1
    nrow = int(np.ceil(n_comp / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7 * ncol, 2.6 * nrow), squeeze=False)
    for idx in range(nrow * ncol):
        ax = axes[idx // ncol][idx % ncol]
        if idx >= n_comp:
            ax.set_axis_off()
            continue
        weights = S[idx, :]
        ax.hist(weights, bins=60, color="#4C72B0", alpha=0.8)
        order = np.argsort(-np.abs(weights))[:top_n]
        top_labels = ", ".join(str(gene_names[j]) for j in order[:min(top_n, 6)])
        ax.set_title("%s (top: %s ...)" % (comp_names[idx], top_labels), fontsize=8)
        ax.axvline(0, color="grey", lw=0.6, ls="--")
        ax.set_ylabel("genes")
    fig.suptitle("Metagene (source) weight distributions", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_stability_index(stability, comp_names, out_dir):
    """Barplot of the per-component stability index."""
    os.makedirs(out_dir, exist_ok=True)
    stability = np.asarray(stability, dtype=float).ravel()
    comp_names = list(comp_names)
    save_path = os.path.join(out_dir, "stability_index.png")

    if np.all(np.isnan(stability)):
        return _empty_panel(save_path, "Stability index not available\nfor this ICA backend.")

    n = len(comp_names)
    fig, ax = plt.subplots(figsize=(max(6, n * 0.35 + 2), 4.5))
    colors = plt.get_cmap("viridis")(np.clip(stability, 0, 1))
    ax.bar(range(n), stability, color=colors)
    ax.set_xticks(range(n)); ax.set_xticklabels(comp_names, rotation=90, fontsize=7)
    ax.set_ylim(0, 1)
    ax.axhline(0.5, color="red", lw=0.8, ls="--", label="0.5")
    ax.set_ylabel("stability index")
    ax.set_title("Component reproducibility (stability index)", fontsize=12)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_component_correlation(A, comp_names, out_dir):
    """Heatmap of pairwise correlation between component activities."""
    os.makedirs(out_dir, exist_ok=True)
    A = _as_matrix(A)  # samples x components
    comp_names = list(comp_names)
    save_path = os.path.join(out_dir, "component_correlation.png")

    n = A.shape[1]
    if n < 2 or A.shape[0] < 3:
        return _empty_panel(save_path, "Need >= 2 components and >= 3 samples\nfor correlation.")

    corr = np.corrcoef(A, rowvar=False)
    fig, ax = plt.subplots(figsize=(max(6, n * 0.6 + 2), max(5, n * 0.6 + 1)))
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(range(n)); ax.set_xticklabels(comp_names, rotation=90, fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels(comp_names, fontsize=7)
    if n <= 20:
        for i in range(n):
            for j in range(n):
                ax.text(j, i, "%.2f" % corr[i, j], ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, shrink=0.8, label="Pearson r")
    ax.set_title("Correlation between component activities", fontsize=12)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


# --------------------------------------------------------------------------- #
# Clinical association (one heatmap per variable type)
# --------------------------------------------------------------------------- #
def _clinic_frame(clinic):
    """Return the clinical table as a DataFrame with sample IDs as index."""
    if not isinstance(clinic, pd.DataFrame):
        clinic = pd.DataFrame(clinic)
    clinic = clinic.copy()
    if "ID_Patient" in clinic.columns:
        clinic = clinic.set_index(clinic["ID_Patient"].astype(str))
    else:
        clinic.index = clinic.index.astype(str)
    return clinic


def _align_clinic_column(clinic, sample_names, col):
    """Return a pandas Series of `col` aligned to sample_names (NaN when missing)."""
    clinic = _clinic_frame(clinic)
    if col not in clinic.columns:
        return None
    series = clinic[col]
    series = series[~series.index.duplicated(keep="first")]
    return series.reindex([str(s) for s in sample_names])


def _classify_columns(clinic, max_cat_levels_for_numeric=6):
    """Split clinical columns into ('continuous', 'categorical') name lists."""
    continuous, categorical = [], []
    for col in clinic.columns:
        if col == "ID_Patient":
            continue
        s = clinic[col]
        numeric = pd.to_numeric(s, errors="coerce")
        n_valid = numeric.notna().sum()
        if n_valid >= 0.8 * s.notna().sum() and numeric.nunique(dropna=True) > max_cat_levels_for_numeric:
            continuous.append(col)
        else:
            nlev = s.astype(str).replace("nan", np.nan).nunique(dropna=True)
            if 2 <= nlev <= 30:
                categorical.append(col)
    return continuous, categorical


def _draw_assoc_heatmap(values, pvals, comp_names, var_names, save_path, title,
                        cbar_label, vmin, vmax, cmap):
    n_comp, n_var = values.shape
    fig, ax = plt.subplots(figsize=(max(6, n_var * 0.9 + 3), max(4, n_comp * 0.4 + 2)))
    im = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yticks(range(n_comp)); ax.set_yticklabels(comp_names, fontsize=8)
    ax.set_xticks(range(n_var)); ax.set_xticklabels(var_names, rotation=45, ha="right", fontsize=8)
    for i in range(n_comp):
        for j in range(n_var):
            stars = _pval_to_stars(pvals[i, j])
            if stars:
                ax.text(j, i, stars, ha="center", va="center", fontsize=8, color="black")
    fig.colorbar(im, ax=ax, shrink=0.7, label=cbar_label)
    ax.set_title(title, fontsize=12)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_clinical_association(A, comp_names, sample_names, clinic, out_dir):
    """Association between component activities and clinical variables.

    Continuous variables -> signed Spearman correlation heatmap.
    Categorical variables -> -log10(p) heatmap (Mann-Whitney if binary,
    Kruskal-Wallis if >2 modalities). Stars mark significance in both.
    Returns {"continuous": path|None, "categorical": path|None}.
    """
    os.makedirs(out_dir, exist_ok=True)
    out = {"continuous": None, "categorical": None}
    if clinic is None:
        return out

    A = _as_matrix(A)  # samples x components
    comp_names = list(comp_names)
    sample_names = [str(s) for s in sample_names]
    clinic = _clinic_frame(clinic)

    A_df = pd.DataFrame(A, index=sample_names, columns=comp_names)
    common = [s for s in sample_names if s in set(clinic.index)]
    if len(common) < 3:
        out["continuous"] = _empty_panel(
            os.path.join(out_dir, "clinical_assoc_continuous.png"),
            "Fewer than 3 samples shared between\nclinic and expression matrix.")
        return out
    A_df = A_df.loc[common]
    clinic = clinic.loc[~clinic.index.duplicated(keep="first")].loc[common]

    continuous, categorical = _classify_columns(clinic)
    n_comp = len(comp_names)

    # --- continuous : Spearman ---
    if continuous:
        vals = np.full((n_comp, len(continuous)), np.nan)
        pvs = np.full((n_comp, len(continuous)), np.nan)
        for j, col in enumerate(continuous):
            xv = pd.to_numeric(clinic[col], errors="coerce").values
            for i, comp in enumerate(comp_names):
                yv = A_df[comp].values
                mask = ~np.isnan(xv) & ~np.isnan(yv)
                if mask.sum() >= 3 and np.unique(xv[mask]).size > 1:
                    r, p = stats.spearmanr(xv[mask], yv[mask])
                    vals[i, j] = r
                    pvs[i, j] = p
        out["continuous"] = _draw_assoc_heatmap(
            vals, pvs, comp_names, continuous,
            os.path.join(out_dir, "clinical_assoc_continuous.png"),
            "Component vs continuous clinical variables (Spearman)",
            "Spearman r", -1, 1, "coolwarm")

    # --- categorical : Mann-Whitney / Kruskal-Wallis -> -log10(p) ---
    if categorical:
        vals = np.full((n_comp, len(categorical)), np.nan)
        pvs = np.full((n_comp, len(categorical)), np.nan)
        for j, col in enumerate(categorical):
            g = clinic[col].astype(str)
            g = g.replace("nan", np.nan)
            for i, comp in enumerate(comp_names):
                yv = A_df[comp]
                groups = [yv[g == lev].values for lev in g.dropna().unique()]
                groups = [grp for grp in groups if len(grp) >= 2]
                if len(groups) < 2:
                    continue
                try:
                    if len(groups) == 2:
                        _, p = stats.mannwhitneyu(groups[0], groups[1], alternative="two-sided")
                    else:
                        _, p = stats.kruskal(*groups)
                except ValueError:
                    p = np.nan
                pvs[i, j] = p
                vals[i, j] = -np.log10(p) if (p is not None and p > 0) else np.nan
        finite = vals[np.isfinite(vals)]
        vmax = float(np.nanmax(finite)) if finite.size else 1.0
        out["categorical"] = _draw_assoc_heatmap(
            vals, pvs, comp_names, categorical,
            os.path.join(out_dir, "clinical_assoc_categorical.png"),
            "Component vs categorical clinical variables (-log10 p)",
            "-log10(p)", 0, max(vmax, 1.3), "magma")

    return out
