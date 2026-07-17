"""Shared utilities for survival analysis plotting — JAMA Surgery style.
Auto-adds median survival and log-rank P to any KM plot. Colors from _jama_palette.

Usage:
    from _surv_utils import km_plot_with_stats, style_ax, JAMA_RC
    from _jama_palette import SURG  # JAMA surgery group colors

    plt.rcParams.update(JAMA_RC)
    fig, ax = plt.subplots(figsize=(6.85, 3.5))
    km_plot_with_stats(ax, data, group_col, time_col, event_col, SURG)
"""
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

def km_plot_with_stats(ax, df, group_col, time_col, event_col, palette,
                       max_t=120, risk_step=24, ylabel='Survival Probability',
                       title='', fontsize=1.0):
    """Kaplan-Meier plot with median survival + log-rank P built in."""
    groups = df[group_col].unique()
    groups = [g for g in groups if len(df[df[group_col]==g]) >= 30]
    groups.sort()
    dur_all, ev_all, grp_all = [], [], []
    medians = {}

    for g in groups:
        sub = df[df[group_col]==g]
        c = palette.get(g, '#333')
        kmf = KaplanMeierFitter()
        kmf.fit(sub[time_col], sub[event_col])
        m = kmf.median_survival_time_
        medians[g] = m
        label = f'{g} (n={len(sub)})'
        label += f', {int(m)}mo' if np.isfinite(m) else ', NR'
        kmf.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.12, lw=1.3, color=c, label=label)
        dur_all.extend(sub[time_col].tolist())
        ev_all.extend(sub[event_col].tolist())
        grp_all.extend([g]*len(sub))

    ax.axhline(0.5, color='#999', ls=':', lw=0.4, alpha=0.5, zorder=0)
    for g in groups:
        m = medians.get(g)
        if m is not None and np.isfinite(m):
            c = palette.get(g, '#333')
            ax.plot([m, m], [0, 0.5], '--', lw=0.6, color=c, alpha=0.5, zorder=0)

    if len(groups) >= 2:
        lr = multivariate_logrank_test(dur_all, ev_all, grp_all)
        p = lr.p_value
        txt = f'Log-rank P {"<0.001" if p < 0.001 else f"={p:.3f}"}'
        ax.text(0.98, 0.98, txt, transform=ax.transAxes,
                fontsize=6.5*fontsize, va='top', ha='right', color='#333', style='italic')

    ax.set_title(title, fontweight='bold', fontsize=8.5*fontsize, color='#1a1a1a', loc='left')
    ax.set_xlabel('Months', fontsize=7.5*fontsize)
    ax.set_ylabel(ylabel, fontsize=7.5*fontsize)
    ax.set_xlim(0, max_t)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(np.arange(0, max_t+1, risk_step))
    ax.legend(frameon=False, fontsize=6.5*fontsize, loc='lower left')

def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.4)
    ax.spines['bottom'].set_linewidth(0.4)
    ax.tick_params(width=0.4)
    ax.set_facecolor('white')

JAMA_RC = {
    'font.family': 'sans-serif', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 7, 'figure.dpi': 300,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.4, 'ytick.major.width': 0.4,
    'axes.edgecolor': '#313131',
    'text.color': '#313131',
    'axes.labelcolor': '#313131',
    'xtick.color': '#374E55',
    'ytick.color': '#374E55',
}

ASO_RC = JAMA_RC  # backward compat
