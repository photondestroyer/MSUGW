"""Figures for the Papers 1 and 2 detailed plan. Numbers come from checks/*.json (made by checks/check_*.py)."""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = f"{HERE}/fig"
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SEQ = ["#86b6ef", "#2a78d6", "#104281"]          # one hue, light to dark (ordered rule strictness / depth class)
plt.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white", "axes.facecolor": SURF, "savefig.dpi": 230})


def box(ax, xy, w, h, text, fc="#f3f2ee", ec=AXIS, lw=0.9, fs=8.5, bold=False, tc=INK):
    x, y = xy
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", color=tc,
            zorder=3, linespacing=1.25)


def arrow(ax, a, b, color=MUTED, lw=1.0, rad=0.0, ls="-"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=9, color=color, lw=lw, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2, zorder=1))


def fig_dag():
    fig, ax = plt.subplots(figsize=(7.1, 3.55))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.06); ax.axis("off"); ax.set_facecolor("white")
    W, H = 0.185, 0.115
    pos = dict(clim=(0.11, 0.86), terr=(0.11, 0.50), geol=(0.11, 0.14), pump=(0.43, 0.88), depth=(0.43, 0.56),
               runon=(0.43, 0.14), veg=(0.72, 0.84), fun=(0.89, 0.44))
    box(ax, pos["clim"], W, H, "Climate\n(rain, demand)")
    box(ax, pos["terr"], W, H, "Terrain position\n(valley or upland)")
    box(ax, pos["geol"], W, H, "Geology and soil")
    box(ax, pos["pump"], W + 0.03, H, "Pumping, canal seepage,\nstream connection", fc="#e6f0fc", ec=BLUE, lw=1.3)
    box(ax, pos["depth"], W, H, "Depth to water\n(exposure)", fc="#e6f0fc", ec=BLUE, lw=1.6, bold=True)
    box(ax, pos["runon"], W + 0.02, H, "Run-on, soil water,\nnutrients", fc="#fdebe3", ec=ORANGE, lw=1.3)
    box(ax, pos["veg"], W, H, "Vegetation type", fc="#e3f6ee", ec=AQUA, lw=1.3)
    box(ax, pos["fun"], W, H + 0.02, "Function held\nthrough drought\n(outcome)", fc="#f3f2ee", ec=INK2, lw=1.6, bold=True)

    def e(k, dx=0.0, dy=0.0):
        return (pos[k][0] + dx, pos[k][1] + dy)

    g = 0.022                                    # gap so that arrow heads stay visible outside the boxes
    arrow(ax, e("clim", W / 2 + g, -0.03), e("depth", -W / 2 - g, 0.035))
    arrow(ax, e("terr", W / 2 + g), e("depth", -W / 2 - g))
    arrow(ax, e("geol", W / 2 + g, 0.03), e("depth", -W / 2 - g, -0.035))
    arrow(ax, e("pump", 0, -H / 2 - g), e("depth", 0, H / 2 + g), color=BLUE, lw=1.6)
    arrow(ax, e("terr", W / 2 + g, -0.04), e("runon", -W / 2 - 0.01 - g, 0.04), color=ORANGE, lw=1.4)
    arrow(ax, e("geol", W / 2 + g), e("runon", -W / 2 - 0.01 - g), color=ORANGE, lw=1.4)
    arrow(ax, e("runon", W / 2 + 0.01 + g, 0.02), e("fun", -W / 2 - g, -0.05), color=ORANGE, lw=1.4)
    arrow(ax, e("depth", W / 2 + g), e("fun", -W / 2 - g, 0.02), color=BLUE, lw=2.2)
    arrow(ax, e("depth", W / 2 + g, 0.045), e("veg", -W / 2 - g, -0.03), color=AQUA, lw=1.4)
    arrow(ax, e("veg", 0.05, -H / 2 - g), e("fun", -0.02, H / 2 + 0.01 + g), color=AQUA, lw=1.4)
    ax.text(0.665, 0.565, "effect wanted", color="#1c5cab", fontsize=8, ha="center", fontweight="bold")
    ax.text(0.80, 0.135, "confounding path:\nlow ground is wetter with\nor without groundwater", color="#b84a1f",
            fontsize=7.6, ha="center", va="center", linespacing=1.2)
    ax.text(0.705, 0.665, "path through\nvegetation type", color="#12805a", fontsize=7.6, ha="center", va="center", linespacing=1.2)
    ax.text(0.43, 1.005, "variation in depth that terrain does not cause", color="#1c5cab", fontsize=7.6, ha="center",
            va="center")
    fig.subplots_adjust(0.01, 0.01, 0.99, 0.99)
    fig.savefig(f"{FIG}/fig1_dag.png"); plt.close(fig)


def fig_support():
    S = json.load(open(f"{HERE}/checks/support.json"))["rows"]
    d = np.arange(1, 11)
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.1, 2.85), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.22))
    wbar = 0.26
    for k, (key, lab) in enumerate((("ref15", "deeper than 15 m"), ("ref30", "deeper than 30 m"), ("ref60", "deeper than 60 m"))):
        v = np.array([100 * r[key] for r in S])
        a.bar(d + (k - 1) * wbar, v, wbar - 0.03, color=SEQ[k], label=lab, zorder=3)
    for i, ytxt in ((9, 36.0), (10, 24.0)):      # the two wettest deciles under the 30 m rule
        v = 100 * S[i - 1]["ref30"]
        a.annotate(f"{v:.1f} %", xy=(i, v + 0.6), xytext=(i, ytxt), ha="center", va="bottom", fontsize=7.3, color=INK2,
                   arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkA=1, shrinkB=0))
    a.set_xticks(d); a.set_xlabel("Aridity decile of the aquifer (1 = driest, 10 = wettest)")
    a.set_ylabel("% of decile area"); a.set_ylim(0, 85)
    a.yaxis.grid(True, color=GRID, lw=0.6, zorder=0); a.set_axisbelow(True)
    a.legend(frameon=False, fontsize=7.6, loc="upper right", title="Mapped water table", title_fontsize=7.6,
             handlelength=1.1, borderaxespad=0.1)
    a.set_title("a  Upper bound on rain-only reference land\n    (never irrigated in LANID, water table deeper than the rule)",
                loc="left", fontsize=8.4, color=INK)
    bot = np.zeros(10)
    for k, (key, lab) in enumerate((("nat_5_10", "5 to 10 m"), ("nat_2_5", "2 to 5 m"), ("nat_le2", "2 m or less"))):
        v = np.array([100 * r[key] for r in S])
        b.bar(d, v, 0.62, bottom=bot, color=SEQ[k], label=lab, zorder=3, edgecolor=SURF, linewidth=0.8)
        bot += v
    b.set_xticks(d); b.set_xlabel("Aridity decile"); b.set_ylim(0, 34)
    b.yaxis.grid(True, color=GRID, lw=0.6, zorder=0); b.set_axisbelow(True)
    h, l = b.get_legend_handles_labels()
    b.legend(h[::-1], l[::-1], frameon=False, fontsize=7.6, loc="upper left", title="Mapped water table",
             title_fontsize=7.6, handlelength=1.1)
    b.set_title("b  Never-irrigated land over shallow water\n    (where natural groundwater use can occur)", loc="left",
                fontsize=8.4, color=INK)
    fig.subplots_adjust(0.065, 0.155, 0.995, 0.845)
    fig.savefig(f"{FIG}/fig2_support.png"); plt.close(fig)


def fig_pipeline():
    fig, ax = plt.subplots(figsize=(7.1, 2.75))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_facecolor("white")
    steps = [("P1.0", "Block\naggregates", "histograms,\nobservation support"),
             ("P1.1", "Reference\npopulation", "depth rule per class,\nring test"),
             ("P1.2", "Rain-only\nexpectation", "bias inside\ndrought months"),
             ("P1.3", "Transfer", "applicability,\nweighted calibration"),
             ("P1.4", "Excess", "bound, model,\nlocal contrast"),
             ("P1.5", "Monthly\nstates", "excess and\ncoupling"),
             ("P1.6", "Source", "irrigation, ground-\nwater, surface water"),
             ("P1.7-8", "Validation", "pumping, pivots,\nother sensors")]
    n = len(steps); w = 0.080; xs = np.linspace(0.062, 0.938, n)
    for (sid, name, chk), x in zip(steps, xs):
        box(ax, (x, 0.62), w, 0.25, f"{sid}\n{name}", fc="#e6f0fc", ec=BLUE, lw=1.0, fs=7.6)
        ax.text(x, 0.42, chk, ha="center", va="top", fontsize=6.6, color=INK2, linespacing=1.2)
    for x0, x1 in zip(xs[:-1], xs[1:]):
        ax.plot([x0 + w / 2 + 0.014, x1 - w / 2 - 0.014], [0.62, 0.62], color=MUTED, lw=1.0, zorder=1)
    gates = [("G1a", (xs[1] + xs[2]) / 2, "reference is clean\nand large enough"),
             ("G1b", (xs[3] + xs[4]) / 2, "unbiased inside\ndrought months"),
             ("G1c", (xs[4] + xs[5]) / 2, "zero volume on\nrain-only land"),
             ("G1d", xs[7], "agrees with pumping\nand other sensors")]
    for g, x, t in gates:
        ax.plot([x, x], [0.79 if g == "G1d" else 0.50, 0.87], color=ORANGE, lw=1.2, ls=(0, (3, 2)), zorder=0)
        ax.scatter([x], [0.905], marker="D", s=44, color=ORANGE, zorder=4)
        ax.text(x - 0.014, 0.905, g, fontsize=7.8, fontweight="bold", color=INK, va="center", ha="right")
    gx = np.linspace(0.14, 0.86, 4)
    for (g, x, t), xx in zip(gates, gx):
        ax.text(xx, 0.085, f"{g}: {t}", ha="center", va="center", fontsize=7.0, color=INK2, linespacing=1.2)
    ax.text(0.5, 0.205, "Gates (orange diamonds): a failed gate stops the next step and triggers the fallback in Section 6",
            ha="center", fontsize=7.4, color=INK2)
    fig.subplots_adjust(0.005, 0.0, 0.995, 1.0)
    fig.savefig(f"{FIG}/fig3_pipeline.png"); plt.close(fig)


def _section(ax, wt_valley_depth, title, sub, raised=False, label_b=True):
    x = np.linspace(0, 10, 400)
    if raised:                                   # flat upland; water table mounded on the left, deep on the right
        surf = np.full_like(x, 6.0)
        wt = np.clip(5.6 - 0.075 * np.clip(x - 2.0, 0, None) ** 2, 1.4, None)
    else:
        surf = 6.0 - 3.2 * np.exp(-((x - 5) / 1.25) ** 2)
        base = surf.min() - wt_valley_depth
        wt = np.minimum(base + 0.10 * np.abs(x - 5) ** 1.1, surf - 0.03)
    ax.add_patch(Polygon(np.c_[np.r_[x, 10, 0], np.r_[surf, 0, 0]], closed=True, fc="#efe9dc", ec="#b5ac98", lw=0.9, zorder=1))
    ax.add_patch(Polygon(np.c_[np.r_[x, 10, 0], np.r_[wt, 0, 0]], closed=True, fc="#cfe2f8", ec="none", zorder=2))
    ax.plot(x, wt, color=BLUE, lw=1.6, zorder=3)
    pts = [(1.5, "A")] + ([(5.0, "B")] if label_b else [(8.5, "C")])
    for xx, lab in pts:
        ys = np.interp(xx, x, surf)
        ax.plot([xx, xx], [ys, ys + 0.55], color="#12805a", lw=1.4, zorder=4)
        ax.scatter([xx], [ys + 0.62], s=42, color=AQUA, zorder=5)
        ax.plot([xx, xx - 0.14, xx, xx + 0.14, xx], [ys, ys - 0.45, ys - 0.8, ys - 0.45, ys], color="#12805a", lw=0.7, zorder=4)
        ax.text(xx + 0.45, ys + 0.62, lab, ha="center", va="center", fontsize=8.6, fontweight="bold", color=INK)
    ax.set_xlim(0, 10); ax.set_ylim(0, 7.6); ax.axis("off")
    ax.set_title(title, fontsize=8.5, loc="left", color=INK, pad=2)
    ax.text(0.0, -0.35, sub, fontsize=7.4, color=INK2, va="top", linespacing=1.25)


def fig_design():
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.45), gridspec_kw=dict(wspace=0.10))
    _section(axs[0], 0.25, "Valley that reaches groundwater", "B: low ground, groundwater in reach\nA: upland, rain only")
    _section(axs[1], 2.4, "Dry valley (water table far below)", "B: low ground, no groundwater\nA: upland, rain only")
    _section(axs[2], 0.0, "Upland with a raised water table",
             "A: upland, groundwater in reach (canal seepage,\ninterdune areas).  C: upland, rain only", raised=True, label_b=False)
    for ax in axs:
        ax.text(9.9, 0.45, "water table", color="#1c5cab", fontsize=7.2, ha="right")
    fig.subplots_adjust(0.01, 0.17, 0.99, 0.91)
    fig.savefig(f"{FIG}/fig4_design.png"); plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_dag(); fig_support(); fig_pipeline(); fig_design()
    print("figures written to", FIG)
