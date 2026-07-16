"""Article-style LaTeX/PDF report for dvr.py.

No physics here: consumes the fit params and eigenvalues dvr.py produces and
renders (1) three figures mirroring the reference article -- the potential energy
curve with its Rydberg-6 fit, the fit residual, and the vibrational spectrum --
and (2) a LaTeX document with three tables, compiled to PDF via pdflatex.
"""
import os
import shutil
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import dvr

BOHR_TO_ANG = 0.52917721067
HARTREE_TO_KCAL = 627.509474

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "dejavuserif",
    "axes.linewidth": 0.8,
    "figure.dpi": 120,
})


# ------------------------------------------------------------------- figures
def _pec_curve(p, rr):
    """Fitted Rydberg-6 relative to the well bottom (well bottom = 0)."""
    return dvr.rydberg6(rr, p["De"], p["Re"], 0.0, p["c1"], p["c2"], p["c3"],
                        p["c4"], p["c5"], p["c6"])


def make_figures(r, v, p, results, A, B, assets):
    os.makedirs(assets, exist_ok=True)
    paths = {}
    E0 = results[0]                      # hartree, relative to minimum
    vrel = (v - p["xe"]) * dvr.CM        # data relative to minimum, cm-1
    rr = np.linspace(A, B, 1000)
    fit_rel = _pec_curve(p, rr) * dvr.CM
    De_cm = p["De"] * dvr.CM

    # --- Fig 1: PEC + Ryd6 fit + vibrational levels ---
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.scatter(r, vrel, s=14, facecolors="none", edgecolors="tab:blue",
               linewidths=0.9, label="ab initio (CSV)", zorder=3)
    ax.plot(rr, fit_rel, color="tab:blue", lw=1.6,
            label="Rydberg-6 fit", zorder=2)
    ax.axhline(De_cm, ls=":", color="0.5", lw=0.9)
    ax.text(B, De_cm, r"  $D_e$", va="center", ha="left", color="0.4",
            fontsize=9)
    ax.axvline(p["Re"], ls=":", color="0.7", lw=0.8)
    # vibrational levels drawn across their classically allowed region
    for e in E0 * dvr.CM:
        inside = rr[fit_rel <= e]
        if inside.size:
            ax.hlines(e, inside.min(), inside.max(), color="tab:green",
                      lw=0.5, alpha=0.6, zorder=1)
    ax.set_xlabel(r"$r$ (bohr)")
    ax.set_ylabel(r"$V(r)-V(R_e)$ (cm$^{-1}$)")
    ax.set_xlim(A, B)
    ax.set_ylim(-0.05 * De_cm, 1.15 * De_cm)

    # lower right is the one region no curve, level or D_e line reaches
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.set_title("Potential energy curve and sixth-degree Rydberg fit")
    fig.tight_layout()
    paths["pec"] = os.path.join(assets, "pec.pdf")
    fig.savefig(paths["pec"]); fig.savefig(paths["pec"][:-4] + ".png", dpi=300)
    plt.close(fig)

    # --- Fig 2: fit error (own panel; a twin axis on Fig 1 buried it) ---
    resid = (_pec_curve(p, r) * dvr.CM) - vrel
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.plot(r, resid, "o-", color="tab:red", ms=3, lw=1.0)
    ax.set_xlabel(r"$r$ (bohr)")
    ax.set_ylabel(r"fit error (cm$^{-1}$)")
    ax.set_xlim(A, B)
    ax.set_title("Rydberg-6 fit error (fit $-$ ab initio), RMS = "
                 f"{p['rms'] * dvr.CM:.2f} " + r"cm$^{-1}$")
    ax.grid(True, ls=":", lw=0.5, alpha=0.6)
    fig.tight_layout()
    paths["err"] = os.path.join(assets, "fiterror.pdf")
    fig.savefig(paths["err"]); fig.savefig(paths["err"][:-4] + ".png", dpi=300)
    plt.close(fig)

    # --- Fig 3: vibrational spectrum (levels ladder + spacing roll-off) ---
    fig, (axl, axc) = plt.subplots(1, 2, figsize=(7.6, 4.4))
    Ecm = E0 * dvr.CM
    for k, e in enumerate(Ecm):
        axl.hlines(e, 0.1, 0.9, color="tab:blue", lw=1.0)
        if k % max(1, len(Ecm) // 15) == 0:
            axl.text(0.93, e, f"$v={k}$", va="center", fontsize=7, color="0.3")
    axl.set_xlim(0, 1.25); axl.set_xticks([])
    axl.set_ylabel(r"$E_v$ (cm$^{-1}$)")
    axl.set_title("Vibrational levels (J=0)")

    vv = np.arange(len(Ecm) - 1)
    dE = np.diff(Ecm)
    axc.plot(vv, dE, "o-", color="tab:purple", ms=4, lw=1.0)
    axc.set_xlabel(r"$v$"); axc.set_ylabel(r"$\Delta E_v=E_{v+1}-E_v$ (cm$^{-1}$)")
    axc.set_title("Level spacing (anharmonic roll-off)")
    axc.grid(True, ls=":", lw=0.5, alpha=0.6)
    fig.tight_layout()
    paths["spec"] = os.path.join(assets, "spectrum.pdf")
    fig.savefig(paths["spec"]); fig.savefig(paths["spec"][:-4] + ".png", dpi=300)
    plt.close(fig)
    return paths


# --------------------------------------------------------------------- LaTeX
def _esc(s):
    return str(s).replace("_", r"\_").replace("&", r"\&")


def _fnum(x, nd=6):
    return f"{x:.{nd}f}"


def _enum(x, nd=6):
    return f"{x:.{nd}e}".replace("e", r"\mathrm{e}")


def make_tex(base, p, results, figs, mu):
    c = dvr.all_constants(results)
    E0cm, E1cm = results[0] * dvr.CM, results[1] * dvr.CM
    n = max(len(E0cm), len(E1cm))
    De_cm = p["De"] * dvr.CM

    L = []
    A = L.append
    A(r"\documentclass[11pt]{article}")
    A(r"\usepackage[a4paper,margin=2.3cm]{geometry}")
    A(r"\usepackage{booktabs,graphicx,longtable,caption,amsmath}")
    A(r"\usepackage[dvipsnames]{xcolor}")
    A(r"\captionsetup{font=small,labelfont=bf}")
    A(r"\renewcommand{\arraystretch}{1.15}")
    A(r"\begin{document}")
    A(r"\begin{center}{\Large\bfseries Rovibrational analysis: %s}\\[2pt]"
      % _esc(base))
    A(r"{\small Sinc-DVR (Colbert--Miller) + sixth-degree Rydberg potential}"
      r"\end{center}")
    A(r"\noindent Reduced mass $\mu = %.5f\ m_e$ (%.5f amu). Grid: %d sinc-DVR "
      r"points on the CSV range. Level count is automatic: all bound states "
      r"below $D_e$." % (mu, mu / dvr.AMU_TO_ME, dvr.NPOINTS))
    A(r"\vspace{4pt}")

    # Table 1 -- Rydberg-6 parameters
    A(r"\begin{table}[h]\centering\caption{Optimized sixth-degree Rydberg "
      r"parameters (atomic units unless noted).}")
    A(r"\begin{tabular}{l r}\toprule Parameter & Value \\\midrule")
    A(r"$R_e$ (bohr) & %s \\" % _fnum(p["Re"], 6))
    A(r"$R_e$ (\AA) & %s \\" % _fnum(p["Re"] * BOHR_TO_ANG, 6))
    A(r"$D_e$ (hartree) & %s \\" % _fnum(p["De"], 8))
    A(r"$D_e$ (cm$^{-1}$) & %s \\" % _fnum(De_cm, 3))
    A(r"$D_e$ (kcal mol$^{-1}$) & %s \\" % _fnum(p["De"] * HARTREE_TO_KCAL, 4))
    A(r"$V(R_e)$ (hartree) & %s \\" % _fnum(p["xe"], 8))
    for k in range(1, 7):
        A(r"$c_%d$ & %s \\" % (k, _fnum(p["c%d" % k], 8)))
    A(r"fit RMS (hartree) & %.3e \\" % p["rms"])
    A(r"\bottomrule\end{tabular}\end{table}")

    # Table 2 -- levels / spectrum / differences
    A(r"\begin{longtable}{r r r r r}")
    A(r"\caption{Vibrational energy levels, spectrum $E_v-E_0$ and level "
      r"differences $\Delta E_v$ (cm$^{-1}$).}\\")
    A(r"\toprule $v$ & $E_v$ (J=0) & $E_v$ (J=1) & $E_v-E_0$ (J=0) & "
      r"$\Delta E_v$ (J=0) \\\midrule\endfirsthead")
    A(r"\toprule $v$ & $E_v$ (J=0) & $E_v$ (J=1) & $E_v-E_0$ (J=0) & "
      r"$\Delta E_v$ (J=0) \\\midrule\endhead")
    for i in range(n):
        e0 = _fnum(E0cm[i], 4) if i < len(E0cm) else "--"
        e1 = _fnum(E1cm[i], 4) if i < len(E1cm) else "--"
        sp = _fnum(E0cm[i] - E0cm[0], 4) if i < len(E0cm) else "--"
        df = _fnum(E0cm[i + 1] - E0cm[i], 4) if i + 1 < len(E0cm) else "--"
        A(r"%d & %s & %s & %s & %s \\" % (i, e0, e1, sp, df))
    A(r"\bottomrule\end{longtable}")

    # Table 3 -- spectroscopic constants
    A(r"\begin{table}[h]\centering\caption{Spectroscopic constants "
      r"(cm$^{-1}$), finite differences of the first levels (Eq.~7 of "
      r"Silva \emph{et al.}, J. Mol. Model. \textbf{24}, 235, 2018).}")
    A(r"\begin{tabular}{l r r}\toprule Constant & J=0 & J=1 \\\midrule")
    A(r"$\omega_e$ & %s & %s \\" % (_fnum(c["we"][0], 4), _fnum(c["we"][1], 4)))
    A(r"$\omega_e x_e$ & %s & %s \\"
      % (_fnum(c["wexe"][0], 5), _fnum(c["wexe"][1], 5)))
    A(r"$\omega_e y_e$ & %.4e & %.4e \\" % (c["weye"][0], c["weye"][1]))
    A(r"\midrule")
    A(r"$B_e$ & \multicolumn{2}{c}{%s} \\" % _fnum(c["Be"], 6))
    A(r"$\alpha_e$ & \multicolumn{2}{c}{%.4e} \\" % c["alfae"])
    A(r"$\gamma_e$ & \multicolumn{2}{c}{%.4e} \\" % c["gamae"])
    A(r"\bottomrule\end{tabular}\end{table}")

    # Figures
    A(r"\begin{figure}[h]\centering\includegraphics[width=0.82\textwidth]{%s}"
      % figs["pec"].replace("\\", "/"))
    A(r"\caption{Potential energy curve: ab initio points, sixth-degree "
      r"Rydberg fit, and the J=0 vibrational levels (green).}\end{figure}")
    A(r"\begin{figure}[h]\centering\includegraphics[width=0.82\textwidth]{%s}"
      % figs["err"].replace("\\", "/"))
    A(r"\caption{Residual of the sixth-degree Rydberg fit, fit $-$ ab initio, "
      r"over the whole grid.}\end{figure}")
    A(r"\begin{figure}[h]\centering\includegraphics[width=0.9\textwidth]{%s}"
      % figs["spec"].replace("\\", "/"))
    A(r"\caption{Vibrational spectrum: energy-level ladder (left) and the "
      r"anharmonic decrease of the level spacing $\Delta E_v$ (right).}"
      r"\end{figure}")
    A(r"\end{document}")

    tex = f"{base}_report.tex"
    with open(tex, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    return tex


# ------------------------------------------------------------------- compile
def _find_pdflatex():
    exe = shutil.which("pdflatex")
    if exe:
        return exe
    guess = os.path.expandvars(
        r"%LOCALAPPDATA%\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe")
    return guess if os.path.exists(guess) else None


def build_report(base, r, v, p, results, A, B, mu):
    """Figures -> LaTeX -> PDF. Degrades gracefully if pdflatex is absent."""
    try:
        figs = make_figures(r, v, p, results, A, B, f"{base}_assets")
    except Exception as e:                       # never let plotting kill a run
        print(f"[dvr] aviso: figuras falharam ({e}); PDF sem figuras.",
              file=sys.stderr)
        figs = {"pec": "", "err": "", "spec": ""}
    tex = make_tex(base, p, results, figs, mu)

    exe = _find_pdflatex()
    if not exe:
        print(f"[dvr] pdflatex nao encontrado; deixei {tex} + figuras para "
              "compilar manualmente.", file=sys.stderr)
        return None
    for _ in range(2):                           # 2 passes settle longtable
        proc = subprocess.run(
            [exe, "-interaction=nonstopmode", "-halt-on-error", tex],
            capture_output=True, text=True)
    pdf = f"{base}_report.pdf"
    if proc.returncode != 0 or not os.path.exists(pdf):
        tail = "\n".join(proc.stdout.splitlines()[-15:])
        print(f"[dvr] pdflatex falhou:\n{tail}", file=sys.stderr)
        return None
    print(f"[dvr] PDF: {pdf}")
    if sys.platform == "win32":
        try:
            os.startfile(pdf)                    # best-effort auto-open
        except OSError:
            pass
    return pdf
