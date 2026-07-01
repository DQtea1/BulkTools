### Independent Component Analysis (ICA) of bulk transcriptomic data ###
# Mirrors the workflow of the user's Jupyter notebooks
# (01_Predimel / 05_Projet_Sarah  ->  04_PCA_ICA/02_ICA.ipynb):
#   from sica.base import StabilizedICA, MSTD
#   MSTD(X.values, m, M, step, n_runs)                 -> part 1 diagnostic
#   sICA = StabilizedICA(n_components, n_runs).fit(X)   -> part 2
#   A = sICA.S_ (metagenes x genes) ; S = pinv(A).T @ X.T (metagenes x samples)
#   scatter_with_marginals(ICA_S.T, clinic, anno_col, x, y, ...)  (user function)
#
# X is oriented samples x genes. Matrices are exchanged with R as plain numpy
# arrays plus name vectors so reticulate round-trips stay deterministic.

import os
import sys
import types


def _silence_tqdm():
    """Neutralize tqdm progress bars used by stabilized-ica.

    Under reticulate there is no Jupyter frontend, so `tqdm.auto` selects the
    ipywidgets variant and raises "ImportError: IProgress not found". We replace
    tqdm everywhere with a disabled console tqdm. Because reticulate keeps a
    single persistent Python session, we also rebind it inside any already
    imported `sica` module (which may have bound the notebook tqdm on an earlier
    failed run before this patch existed).
    """
    try:
        from tqdm.std import tqdm as _std_tqdm
    except Exception:
        return

    class _SilentTqdm(_std_tqdm):
        def __init__(self, *args, **kwargs):
            kwargs["disable"] = True          # never render a bar
            super().__init__(*args, **kwargs)

    def _silent_trange(*args, **kwargs):
        kwargs["disable"] = True
        return _SilentTqdm(range(*args), **kwargs)

    # Make every tqdm entry point resolve to the silent, console-only bar.
    for name in ("tqdm", "tqdm.auto", "tqdm.autonotebook", "tqdm.notebook", "tqdm.std"):
        mod = sys.modules.get(name)
        if mod is None:
            mod = types.ModuleType(name)
            sys.modules[name] = mod
        try:
            mod.tqdm = _SilentTqdm
            mod.trange = _silent_trange
        except Exception:
            pass

    # Rebind tqdm in any sica module that already captured it.
    for name, mod in list(sys.modules.items()):
        if (name == "sica" or name.startswith("sica.")) and mod is not None:
            for attr in ("tqdm", "trange"):
                if getattr(mod, attr, None) is not None:
                    try:
                        setattr(mod, attr, _SilentTqdm if attr == "tqdm" else _silent_trange)
                    except Exception:
                        pass


def _patch_sklearn_agglomerative():
    """Compatibility shim for stabilized-ica on recent scikit-learn.

    stabilized-ica 2.0 calls ``AgglomerativeClustering(affinity=...)``. scikit-learn
    renamed ``affinity`` to ``metric`` (deprecated in 1.2, removed in 1.4), so the
    call raises ``TypeError: ... unexpected keyword argument 'affinity'``. We wrap
    the estimator to translate ``affinity`` -> ``metric`` and rebind it inside any
    already-imported ``sica`` module (persistent reticulate session).
    """
    try:
        import inspect
        import sklearn.cluster as _cl
    except Exception:
        return

    orig = _cl.AgglomerativeClustering
    if "affinity" in inspect.signature(orig.__init__).parameters:
        return  # installed sklearn still accepts affinity; nothing to do

    base = getattr(orig, "_sica_orig", orig)

    class _AggCompat(base):
        def __init__(self, *args, affinity=None, **kwargs):
            if affinity is not None and "metric" not in kwargs:
                kwargs["metric"] = affinity
            super().__init__(*args, **kwargs)

    _AggCompat._sica_orig = base
    _cl.AgglomerativeClustering = _AggCompat
    for name, mod in list(sys.modules.items()):
        if (name == "sica" or name.startswith("sica.")) and mod is not None:
            if getattr(mod, "AgglomerativeClustering", None) is not None:
                try:
                    mod.AgglomerativeClustering = _AggCompat
                except Exception:
                    pass


_silence_tqdm()

import numpy as np
import pandas as pd
import matplotlib

# Headless rendering (the app runs without an X server).
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
    if isinstance(x, pd.DataFrame):
        return x.to_numpy(dtype=float)
    return np.asarray(x, dtype=float)


def _expr_df(expr, gene_names, sample_names):
    """Return X as a DataFrame of shape (samples, genes), like the notebooks."""
    mat = _as_matrix(expr)  # genes x samples (as passed from R)
    gene_names = list(gene_names)
    sample_names = list(sample_names)
    if mat.shape == (len(gene_names), len(sample_names)):
        mat = mat.T
    elif mat.shape != (len(sample_names), len(gene_names)):
        raise ValueError(
            "Expression matrix shape %s matches neither genes x samples (%d x %d) "
            "nor samples x genes." % (mat.shape, len(gene_names), len(sample_names)))
    return pd.DataFrame(mat, index=[str(s) for s in sample_names], columns=list(gene_names))


def _stability_indexes(sICA, n_components):
    for attr in ("stability_indexes_", "stability_index_", "stability_"):
        if hasattr(sICA, attr):
            stab = np.asarray(getattr(sICA, attr), dtype=float).ravel()
            if stab.size == n_components:
                return stab
    return np.full(int(n_components), np.nan)


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
# Part 1 : Most Stable Transcriptomic Dimension (sica.base.MSTD)
# --------------------------------------------------------------------------- #
def compute_mstd(expr, gene_names, sample_names, m, M, step, n_runs, out_dir):
    """Run MSTD exactly as in the notebook: MSTD(X.values, m, M, step, n_runs)."""
    from sica.base import MSTD
    _silence_tqdm()  # rebind tqdm inside sica now that it is imported
    _patch_sklearn_agglomerative()  # affinity -> metric shim for recent sklearn

    X = _expr_df(expr, gene_names, sample_names)  # samples x genes
    m, M, step, n_runs = int(m), int(M), int(step), int(n_runs)
    os.makedirs(out_dir, exist_ok=True)
    save_path = os.path.join(out_dir, "MSTD_plot.png")

    # The number of components cannot exceed the number of samples.
    M = min(M, X.shape[0] - 1, X.shape[1])
    if M < m:
        return _empty_panel(save_path,
                            "Invalid dimension range after clamping (min=%d, max=%d).\n"
                            "Max is limited by the number of samples." % (m, M))

    MSTD(X.values, m=m, M=M, step=step, n_runs=n_runs)
    fig = plt.gcf()
    fig.set_size_inches(13, 6)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return {"save_path": save_path}


# --------------------------------------------------------------------------- #
# Part 2 : stabilized ICA (same A / S convention as the notebook)
# --------------------------------------------------------------------------- #
def run_stabilized_ica(expr, gene_names, sample_names, n_components, n_runs, n_jobs=-1):
    """Fit StabilizedICA and reproduce the notebook's A / S matrices.

      A = sICA.S_                       (metagenes x genes)   -> "A_matrix"
      S = pinv(A).T @ X.T               (metagenes x samples) -> "S_matrix"
    """
    from sica.base import StabilizedICA
    _silence_tqdm()  # rebind tqdm inside sica now that it is imported
    _patch_sklearn_agglomerative()  # affinity -> metric shim for recent sklearn

    X = _expr_df(expr, gene_names, sample_names)  # samples x genes
    K = int(n_components)

    sICA = StabilizedICA(n_components=K, n_runs=int(n_runs), plot=False, n_jobs=int(n_jobs))
    sICA.fit(X.values)

    # Notebook convention (cell "Projection de nos samples dans l'espace des IC"):
    A_np = np.asarray(sICA.S_, dtype=float)             # metagenes x genes
    S_np = np.linalg.pinv(A_np).T @ X.values.T          # metagenes x samples

    stab = _stability_indexes(sICA, A_np.shape[0])
    comp_names = ["Metagene%d" % i for i in range(A_np.shape[0])]
    return {
        "A": A_np,                     # metagenes x genes   (gene weights)
        "S": S_np,                     # metagenes x samples (sample activities)
        "stability": stab,
        "comp_names": comp_names,
        "gene_names": list(gene_names),
        "sample_names": [str(s) for s in sample_names],
    }


# --------------------------------------------------------------------------- #
# scatter_with_marginals : verbatim from the user's notebook (matplotlib)
# --------------------------------------------------------------------------- #
def scatter_with_marginals(
    df_xy: pd.DataFrame,
    df_anno: pd.DataFrame,
    anno_col: str,
    *,
    x: str = "x",
    y: str = "y",
    key: str = None,                 # if None, join on index; else merge on this column
    palette: dict = None,            # optional: {category: color}
    shape_list: list = None,         # optional: list of shapes to use for the annotations
    hue_order: list = None,          # optional: order of legend categories
    s: float = 36,
    alpha: float = 0.9,
    legend_title: str = None,
    legend_loc: str = "best",
    figsize: tuple = (7.5, 6.5),
    dpi: int = 120,
    # Marginal histogram controls
    bins=30,                         # or 'auto', 'fd', etc.
    bins_x=None,
    bins_y=None,
    norm_marginals: bool = False,    # if True, each bin sums to 1 (proportions)
    hist_alpha: float = 0.95,
    shape_by: str = None,            # column controlling marker shapes
):
    """
    Scatter of df_xy[x] vs df_xy[y] colored by df_anno[anno_col] (categorical),
    with marginal stacked histograms (top for X, right for Y).
    """
    # --- Align data ---
    if key is None:
        merged = df_xy.join(df_anno[[anno_col] + ([shape_by] if shape_by else [])], how="inner")
    else:
        cols = [key, anno_col] + ([shape_by] if shape_by else [])
        merged = (
            df_xy.reset_index(drop=True)
                .merge(df_anno[cols], on=key, how="inner")
        )
    if merged.empty:
        raise ValueError("No overlapping rows between df_xy and df_anno.")

    # --- Categories & colors ---
    cats = pd.Categorical(merged[anno_col].astype("category"))
    categories = list(cats.categories)
    if hue_order is not None:
        categories = [c for c in hue_order if c in categories]

    if palette is None:
        base_colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
        if len(base_colors) < len(categories):
            cmap = plt.get_cmap("tab20")
            extra = [cmap(i % cmap.N) for i in range(len(categories) - len(base_colors))]
            colors = base_colors + extra
        else:
            colors = base_colors[:len(categories)]
        color_map = dict(zip(categories, colors))
    else:
        color_map = {c: palette.get(c, "lightgray") for c in categories}

    # --- Shapes ---
    if shape_by is not None:
        merged["_shape_key"] = merged[shape_by].astype("string").fillna("NA").str.strip()
        shape_classes = pd.unique(merged["_shape_key"])
        if shape_list is None:
            _markers = ['*', 'X', 'o', 's', '^', 'D', 'v', '<', '>', 'P']
        else:
            _markers = shape_list
        shape_map = {c: _markers[i % len(_markers)] for i, c in enumerate(shape_classes)}
    else:
        merged["_shape_key"] = "__all__"
        shape_classes = ["__all__"]
        shape_map = {"__all__": "o"}

    # --- Layout ---
    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = fig.add_gridspec(
        2, 2, width_ratios=(4, 1.1), height_ratios=(1.1, 4),
        wspace=0.06, hspace=0.06
    )
    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scatter)

    # --- Scatter ---
    done = set()  # ensure one label per color category
    for cat in categories:
        sub_cat = merged[cats == cat]
        for j, (shp, sub) in enumerate(sub_cat.groupby("_shape_key", dropna=False)):
            label = str(cat) if (cat not in done and j == 0) else None
            ax_scatter.scatter(
                sub[x], sub[y],
                s=s, alpha=alpha, label=label,
                c=[color_map[cat]],
                marker=shape_map.get(shp, "o"),
                edgecolors="k", linewidths=0.3
            )
        done.add(cat)

    ax_scatter.set_xlabel(x)
    ax_scatter.set_ylabel(y)
    ax_scatter.set_title('Scatter colored by "%s"' % anno_col)
    ax_scatter.grid(True, linestyle=":", linewidth=0.6, alpha=0.6)

    leg1 = ax_scatter.legend(title=legend_title or anno_col, frameon=False, loc=legend_loc)

    if shape_by is not None:
        from matplotlib.lines import Line2D
        ax_scatter.add_artist(leg1)
        shape_handles = [
            Line2D([0], [0], marker=shape_map[c], linestyle="None", color="k", markersize=7)
            for c in shape_classes
        ]
        ax_scatter.legend(shape_handles, list(shape_classes), title=shape_by, frameon=False, loc="lower left")

    # --- Bins (shared across categories) ---
    x_vals = merged[x].to_numpy()
    y_vals = merged[y].to_numpy()
    edges_x = np.histogram_bin_edges(x_vals, bins=bins_x or bins)
    edges_y = np.histogram_bin_edges(y_vals, bins=bins_y or bins)
    nbx = len(edges_x) - 1
    nby = len(edges_y) - 1

    counts_x = np.zeros((len(categories), nbx), dtype=float)
    for i, cat in enumerate(categories):
        counts_x[i], _ = np.histogram(merged.loc[cats == cat, x], bins=edges_x)
    totals_x = counts_x.sum(axis=0)
    counts_x_plot = (counts_x / totals_x.clip(min=1)) if norm_marginals else counts_x

    bottom_x = np.zeros(nbx, dtype=float)
    for i, cat in enumerate(categories):
        heights = counts_x_plot[i]
        ax_histx.bar(
            edges_x[:-1], heights, width=np.diff(edges_x),
            align="edge", bottom=bottom_x, color=color_map[cat],
            edgecolor="none", alpha=hist_alpha
        )
        bottom_x += heights

    counts_y = np.zeros((len(categories), nby), dtype=float)
    for i, cat in enumerate(categories):
        counts_y[i], _ = np.histogram(merged.loc[cats == cat, y], bins=edges_y)
    totals_y = counts_y.sum(axis=0)
    counts_y_plot = (counts_y / totals_y.clip(min=1)) if norm_marginals else counts_y

    left_y = np.zeros(nby, dtype=float)
    for i, cat in enumerate(categories):
        widths = counts_y_plot[i]
        ax_histy.barh(
            edges_y[:-1], widths, height=np.diff(edges_y),
            align="edge", left=left_y, color=color_map[cat],
            edgecolor="none", alpha=hist_alpha
        )
        left_y += widths

    # --- Cosmetics for marginals ---
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histx.tick_params(axis="y", left=False, labelleft=False)
    ax_histy.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    for ax in (ax_histx,):
        for sp in ("right", "top", "left"):
            ax.spines[sp].set_visible(False)
    for sp in ("right", "top", "bottom"):
        ax_histy.spines[sp].set_visible(False)

    if norm_marginals:
        ax_histx.set_ylim(0, 1)
        ax_histy.set_xlim(0, 1)

    return fig, {"scatter": ax_scatter, "histx": ax_histx, "histy": ax_histy}, merged, color_map


# --------------------------------------------------------------------------- #
# Plot wrappers (save a PNG, return its path) driven from R
# --------------------------------------------------------------------------- #
def _plain_scatter_marginals(df_xy, x, y, save_path):
    """Uncolored scatter + marginal histograms (used when no clinic is given)."""
    xv = df_xy[x].to_numpy()
    yv = df_xy[y].to_numpy()
    fig = plt.figure(figsize=(7.5, 6.5))
    gs = fig.add_gridspec(2, 2, width_ratios=(4, 1.1), height_ratios=(1.1, 4),
                          wspace=0.06, hspace=0.06)
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)
    ax.scatter(xv, yv, s=32, alpha=0.85, color="#1F30C9", edgecolors="k", linewidths=0.3)
    ax.set_xlabel(x); ax.set_ylabel(y)
    ax.axhline(0, color="grey", lw=0.6, ls="--"); ax.axvline(0, color="grey", lw=0.6, ls="--")
    ax.set_title("Sample activities: %s vs %s" % (x, y))
    ax_top.hist(xv, bins=30, color="#1F30C9", alpha=0.7)
    ax_right.hist(yv, bins=30, orientation="horizontal", color="#1F30C9", alpha=0.7)
    ax_top.tick_params(labelbottom=False); ax_right.tick_params(labelleft=False)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_scatter_marginals(S, comp_names, sample_names, comp_x, comp_y, out_dir,
                           clinic=None, color_col=None, bins=15):
    """Scatter of two metagene activities with marginals (user's function).

    `S` is the metagenes x samples activity matrix; it is transposed to
    samples x metagenes and passed to `scatter_with_marginals` as in the notebook.
    """
    os.makedirs(out_dir, exist_ok=True)
    S = _as_matrix(S)  # metagenes x samples
    comp_names = list(comp_names)
    sample_names = [str(s) for s in sample_names]
    save_path = os.path.join(out_dir, "scatter_%s_vs_%s.png" % (str(comp_x), str(comp_y)))

    if comp_x not in comp_names or comp_y not in comp_names:
        return _empty_panel(save_path, "Component(s) not found:\n%s / %s" % (comp_x, comp_y))

    df_xy = pd.DataFrame(S.T, index=sample_names, columns=comp_names)  # samples x metagenes

    has_clinic = clinic is not None and color_col is not None and str(color_col) != ""
    if has_clinic:
        clinic_df = _clinic_frame(clinic)
        if str(color_col) not in clinic_df.columns:
            return _plain_scatter_marginals(df_xy, comp_x, comp_y, save_path)
        try:
            fig, _, _, _ = scatter_with_marginals(
                df_xy, clinic_df, anno_col=str(color_col),
                x=comp_x, y=comp_y, bins=int(bins), legend_title=str(color_col), s=30,
            )
        except ValueError as e:
            return _empty_panel(save_path, str(e))
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    return _plain_scatter_marginals(df_xy, comp_x, comp_y, save_path)


def plot_activity_heatmap(S, comp_names, sample_names, out_dir):
    """Heatmap of metagene activities (metagenes x samples)."""
    os.makedirs(out_dir, exist_ok=True)
    S = _as_matrix(S)  # metagenes x samples
    comp_names = list(comp_names)
    sample_names = [str(s) for s in sample_names]
    save_path = os.path.join(out_dir, "activity_heatmap.png")

    n_comp, n_samp = S.shape
    vmax = np.nanmax(np.abs(S)) if S.size else 1.0
    fig, ax = plt.subplots(figsize=(max(6, n_samp * 0.18 + 3), max(4, n_comp * 0.4 + 2)))
    im = ax.imshow(S, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(n_comp)); ax.set_yticklabels(comp_names, fontsize=8)
    if n_samp <= 60:
        ax.set_xticks(range(n_samp)); ax.set_xticklabels(sample_names, rotation=90, fontsize=6)
    else:
        ax.set_xticks([]); ax.set_xlabel("%d samples" % n_samp)
    ax.set_title("Metagene activity per sample", fontsize=12)
    fig.colorbar(im, ax=ax, shrink=0.7, label="activity")
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_source_distributions(A, comp_names, gene_names, out_dir, top_n=20):
    """Gene-weight distribution and top genes for each metagene (A = metagenes x genes)."""
    os.makedirs(out_dir, exist_ok=True)
    A = _as_matrix(A)  # metagenes x genes
    comp_names = list(comp_names)
    gene_names = list(gene_names)
    top_n = int(top_n)
    save_path = os.path.join(out_dir, "metagene_weights.png")

    n_comp = A.shape[0]
    ncol = 2 if n_comp > 1 else 1
    nrow = int(np.ceil(n_comp / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7 * ncol, 2.6 * nrow), squeeze=False)
    for idx in range(nrow * ncol):
        ax = axes[idx // ncol][idx % ncol]
        if idx >= n_comp:
            ax.set_axis_off()
            continue
        weights = A[idx, :]
        ax.hist(weights, bins=60, color="#4C72B0", alpha=0.8)
        order = np.argsort(-np.abs(weights))
        top_labels = ", ".join(str(gene_names[j]) for j in order[:min(top_n, 6)])
        ax.set_title("%s (top: %s ...)" % (comp_names[idx], top_labels), fontsize=8)
        ax.axvline(0, color="grey", lw=0.6, ls="--")
        ax.set_ylabel("genes")
    fig.suptitle("Metagene gene-weight distributions", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_stability_index(stability, comp_names, out_dir):
    """Barplot of the per-metagene stability index."""
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
    ax.set_title("Metagene reproducibility (stability index)", fontsize=12)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_component_correlation(S, comp_names, out_dir):
    """Heatmap of pairwise correlation between metagene activities (S = metagenes x samples)."""
    os.makedirs(out_dir, exist_ok=True)
    S = _as_matrix(S)  # metagenes x samples
    comp_names = list(comp_names)
    save_path = os.path.join(out_dir, "metagene_correlation.png")

    n = S.shape[0]
    if n < 2 or S.shape[1] < 3:
        return _empty_panel(save_path, "Need >= 2 metagenes and >= 3 samples\nfor correlation.")

    corr = np.corrcoef(S)  # metagenes x metagenes
    fig, ax = plt.subplots(figsize=(max(6, n * 0.6 + 2), max(5, n * 0.6 + 1)))
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(range(n)); ax.set_xticklabels(comp_names, rotation=90, fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels(comp_names, fontsize=7)
    if n <= 20:
        for i in range(n):
            for j in range(n):
                ax.text(j, i, "%.2f" % corr[i, j], ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, shrink=0.8, label="Pearson r")
    ax.set_title("Correlation between metagene activities", fontsize=12)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return save_path


# --------------------------------------------------------------------------- #
# Clinical association (one heatmap per variable type)
# --------------------------------------------------------------------------- #
def _clinic_frame(clinic):
    if not isinstance(clinic, pd.DataFrame):
        clinic = pd.DataFrame(clinic)
    clinic = clinic.copy()
    if "ID_Patient" in clinic.columns:
        clinic = clinic.set_index(clinic["ID_Patient"].astype(str))
    else:
        clinic.index = clinic.index.astype(str)
    return clinic


def _classify_columns(clinic, max_cat_levels_for_numeric=6):
    continuous, categorical = [], []
    for col in clinic.columns:
        if col == "ID_Patient":
            continue
        s = clinic[col]
        numeric = pd.to_numeric(s, errors="coerce")
        n_valid = numeric.notna().sum()
        if n_valid >= 0.8 * max(s.notna().sum(), 1) and numeric.nunique(dropna=True) > max_cat_levels_for_numeric:
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


def plot_clinical_association(S, comp_names, sample_names, clinic, out_dir):
    """Association between metagene activities and clinical variables.

    `S` is metagenes x samples; activities are transposed to samples x metagenes.
    Continuous variables -> signed Spearman correlation heatmap.
    Categorical variables -> -log10(p) heatmap (Mann-Whitney if binary,
    Kruskal-Wallis if >2 modalities). Returns {"continuous": path|None, "categorical": path|None}.
    """
    os.makedirs(out_dir, exist_ok=True)
    out = {"continuous": None, "categorical": None}
    if clinic is None:
        return out

    S = _as_matrix(S)  # metagenes x samples
    comp_names = list(comp_names)
    sample_names = [str(s) for s in sample_names]
    clinic = _clinic_frame(clinic)

    A_df = pd.DataFrame(S.T, index=sample_names, columns=comp_names)  # samples x metagenes
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
            "Metagene vs continuous clinical variables (Spearman)",
            "Spearman r", -1, 1, "coolwarm")

    if categorical:
        vals = np.full((n_comp, len(categorical)), np.nan)
        pvs = np.full((n_comp, len(categorical)), np.nan)
        for j, col in enumerate(categorical):
            g = clinic[col].astype(str).replace("nan", np.nan)
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
            "Metagene vs categorical clinical variables (-log10 p)",
            "-log10(p)", 0, max(vmax, 1.3), "magma")

    return out
