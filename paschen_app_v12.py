# -*- coding: utf-8 -*-
# ## ##########################################################################
#
# Authors: Miguel Nahuatlato
# License: MIT
#
# ========================================================
# Calculadora Ley de Paschen + Simulador Re-Gassing RECI — v12
# ========================================================
#
# Dependencias:
#    pip install flet matplotlib numpy sympy
#
# Uso:
#    py paschen_app_v12.py
#========================================================
# ## ##########################################################################

import flet as ft
import math, os, tempfile, shutil, threading, time
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN 1 — FÍSICA Y MATEMÁTICA
# ═══════════════════════════════════════════════════════════════════════════════

IONIZATION_DATA: dict[str, dict] = {
    'H2':  {'A': 5,   'B': 130},
    'N2':  {'A': 12,  'B': 342},
    'CO2': {'A': 20,  'B': 466},
    'He':  {'A': 3,   'B': 34},
    'Hg':  {'A': 20,  'B': 370},
    'air': {'A': 15,  'B': 365},
}
GASES_PERSONALIZADOS: dict[str, dict] = {}

PESOS_MOLECULARES = {
    'CO2': 44.01, 'N2': 28.014, 'H2': 2.016, 'He': 4.003,
    'air': 28.97,  'Hg': 200.59,
}

CONVERSIONES_A_PA = {
    'Pa': 1.0, 'Torr': 133.322, 'atm': 101325.0,
    'mbar': 100.0, 'psi': 6894.76, 'mmHg': 133.322,
    'kg/cm²': 98066.5, 'inHg': 3386.39,
}

SIM_DEF = {
    'longitud_tubo': 140.0, 'diametro_tubo': 8.0,
    'vol_tubo_manual': 5500.0, 'dist_electrodos': 125.0,
    'mezcla': {'CO2': 10.0, 'N2': 10.0, 'H2': 5.0, 'He': 75.0},
    'p_obj': 15.0, 'p_atm': 760.0, 'p_vac': 0.001,
    'v_ign': 24.0, 'v_op': 18.0, 'i_max': 28.0,
    'tasa_ev': 25.0, 'tasa_ll': 0.15,
    'vol_cuerpo': 50.0,
    'lg_gas': 100.0, 'dg_gas': 0.635,
    'lg_bomba': 100.0, 'dg_bomba': 1.27,
    'lg_laser': 20.0,  'dg_laser': 0.635,
    'gamma': 0.02,
}

DESKTOP  = os.path.join(os.path.expanduser("~"), "Desktop")
TEMP_DIR = tempfile.mkdtemp(prefix="Paschen_")
EXPORT_DIR = DESKTOP if os.path.isdir(DESKTOP) else TEMP_DIR

# ── Utilidades de gas ──────────────────────────────────────────────────────────

def todos_los_gases() -> dict:
    return {**IONIZATION_DATA, **GASES_PERSONALIZADOS}

def norm_gas(n: str) -> str:
    m = {'HE':'He','H2':'H2','N2':'N2','CO2':'CO2','HG':'Hg','AIR':'air','AIRE':'air'}
    return m.get(n.upper(), n.strip())

def norm_unidad(u: str) -> str:
    m = {'PA':'Pa','TORR':'Torr','ATM':'atm','MBAR':'mbar','PSI':'psi',
         'MMHG':'mmHg','KG/CM2':'kg/cm²','KG/CM²':'kg/cm²','KGCM2':'kg/cm²',
         'AT':'kg/cm²','INHG':'inHg','IN HG':'inHg','IN_HG':'inHg'}
    r = m.get(u.upper())
    if r is None:
        raise ValueError(f"Unidad '{u}' no válida. Usa: Pa Torr atm mbar psi mmHg inHg kg/cm²")
    return r

def convertir_presion(valor: float, unidad: str) -> dict:
    u   = norm_unidad(unidad)
    vpa = valor * CONVERSIONES_A_PA[u]
    return {k: vpa / f for k, f in CONVERSIONES_A_PA.items()}

def calcular_AB_mezcla(pct: dict) -> tuple:
    td   = todos_los_gases()
    suma = sum(pct.values())
    fr   = {g: p/suma for g, p in pct.items()}
    Am   = sum(fr[g]*td[g]['A'] for g in fr)
    Bm   = sum(fr[g]*td[g]['B'] for g in fr)
    return Am, Bm, fr, pct

# ── Ley de Paschen ─────────────────────────────────────────────────────────────

def paschen_Vb(pd: float, A: float, B: float, gamma: float) -> float:
    """Vb = B·pd / [ln(A·pd) − ln(ln(1+1/γ))]  — condición: A·pd > ln(1+1/γ)"""
    umbral = math.log(1.0 + 1.0/gamma)
    if A*pd <= umbral:
        return float('nan')
    den = math.log(A*pd) - math.log(umbral)
    return (B*pd)/den if den > 0 else float('nan')

def paschen_pd_min(A: float, gamma: float) -> float:
    """pd_min = e·ln(1+1/γ)/A"""
    return math.e * math.log(1.0+1.0/gamma) / A

def paschen_Vb_min(A: float, B: float, gamma: float) -> float:
    return B * paschen_pd_min(A, gamma)

def solve_for_p(Vb: float, A: float, B: float, d: float, gamma: float,
                tol: float=1e-8, max_iter: int=3000) -> float:
    """Bisección numérica — rama derecha (p > p_min)."""
    pdm  = paschen_pd_min(A, gamma)
    Vbm  = paschen_Vb_min(A, B, gamma)
    if Vb < Vbm*(1-1e-6):
        raise ValueError(f"Vb={Vb:.4g} V < Vb_min={Vbm:.4g} V — imposible físicamente.")
    def f(p):
        v = paschen_Vb(p*d, A, B, gamma)
        return v-Vb if math.isfinite(v) else float('nan')
    pm = pdm/d
    lo, hi = pm, pm*1e5
    while not math.isfinite(f(hi)) or f(hi)<0:
        hi *= 10
        if hi > 1e12: raise ValueError("No se encontró solución en la rama derecha.")
    flo = f(lo)
    if not (flo * f(hi) < 0):
        lo2 = (math.log(1+1/gamma)/A)/d*1.001
        fl2, fh2 = f(lo2), f(pm)
        if math.isfinite(fl2) and math.isfinite(fh2) and fl2*fh2 < 0:
            lo, hi, flo = lo2, pm, fl2
        else:
            raise ValueError("Sin solución — verifica parámetros.")
    for _ in range(max_iter):
        pm2 = (lo+hi)/2.0
        fm  = f(pm2)
        if not math.isfinite(fm): hi = pm2; continue
        if abs(fm) < tol: return pm2
        if flo*fm < 0: hi = pm2
        else: lo, flo = pm2, fm
    return (lo+hi)/2.0

# ── Física del simulador ───────────────────────────────────────────────────────

def vol_manifold(vc, lg, dg, lb, db, ll, dl):
    def vt(d, l): return math.pi*(d/2)**2*l
    return vc + vt(dg,lg) + vt(db,lb) + vt(dl,ll)

def masa_gas(vol_cm3: float, p_torr: float, pct: dict):
    R, T = 62.364, 298.0
    n    = (p_torr * vol_cm3/1000) / (R*T)
    fr   = {g: v/sum(pct.values()) for g,v in pct.items()}
    pm   = sum(fr[g]*PESOS_MOLECULARES.get(g,28.0) for g in fr)
    return n*pm*1000, n, pm, fr          # mg, mol, g/mol, fracs

def sim_evacuacion(p0, pf, tasa=5.0):
    """
    Simulación COMPLETA de evacuación según procedimiento profesional:
    
    FASE A: PREPARACIÓN Y LIMPIEZA (Bake-out)
    - Vacío inicial hasta 10⁻³ Torr
    - Prueba de estanqueidad (60 segundos)
    
    FASE B: FLUSHING (Barrido)
    - 3 ciclos de: llenar a 50 Torr → evacuar
    - Elimina aire residual de recovecos
    
    Retorna arrays de tiempo y presión con las fases marcadas
    """
    tau = p0/tasa
    ts, ps, fases = [], [], []
    t = 0.0
    
    # ═════════════════════════════════════════════════════════════════════
    # FASE A: EVACUACIÓN INICIAL (Vacío hasta 10⁻³ Torr)
    # ═════════════════════════════════════════════════════════════════════
    while True:
        p_actual = p0 * math.exp(-t/tau)
        if p_actual <= pf or t >= 700:
            break
        ts.append(t)
        ps.append(p_actual)
        fases.append("FASE A: EVACUACIÓN")
        t += 1.0
    
    # Asegurar que llegamos al vacío objetivo
    ts.append(t)
    ps.append(pf)
    fases.append("FASE A: EVACUACIÓN")
    
    # ═════════════════════════════════════════════════════════════════════
    # PRUEBA DE ESTANQUEIDAD (60 segundos en vacío)
    # ═════════════════════════════════════════════════════════════════════
    for i in range(60):
        t += 1.0
        ts.append(t)
        ps.append(pf)
        fases.append("PRUEBA FUGAS")
    
    # ═════════════════════════════════════════════════════════════════════
    # FASE B: FLUSHING (3 ciclos de barrido)
    # Procedimiento: inyectar Premix XX hasta 50 Torr → evacuar completamente
    # Repetir 3 veces. Tiempo de llenado de flush ≈ 80 s (válvula de aguja)
    # Evacuación usa misma constante tau que la evacuación principal.
    # ═════════════════════════════════════════════════════════════════════
    p_flush = 50.0   # Torr
    N_fill  = 80     # segundos para subir a p_flush (válvula de aguja manual)

    for ciclo in range(3):
        # ── Llenado lineal hasta p_flush ──────────────────────────────────
        for i in range(N_fill):
            t += 1.0
            p_ciclo = pf + (p_flush - pf) * (i / float(N_fill - 1))
            ts.append(t); ps.append(p_ciclo)
            fases.append(f"FASE B: FLUSH #{ciclo+1} ↑")

        # ── Mantener 10 s a p_flush ───────────────────────────────────────
        for _ in range(10):
            t += 1.0
            ts.append(t); ps.append(p_flush)
            fases.append(f"FASE B: FLUSH #{ciclo+1} EN {p_flush:.0f} Torr")

        # ── Evacuación exponencial completa: p_flush → pf usando tau principal
        # Calculamos el tiempo exacto necesario y generamos siempre una curva suave
        # con un máximo de 1500 puntos para no ralentizar la UI.
        t_ev_flush = tau * math.log(max(p_flush / pf, 1.0))
        n_pts = max(2, min(int(t_ev_flush) + 1, 1500))
        step  = t_ev_flush / n_pts
        for k in range(1, n_pts + 1):
            t_elapsed = k * step
            p_curr = p_flush * math.exp(-t_elapsed / tau)
            if p_curr < pf:
                p_curr = pf
            t += step
            ts.append(t); ps.append(p_curr)
            fases.append(f"FASE B: FLUSH #{ciclo+1} ↓")

        # ── Estabilización en vacío antes del siguiente ciclo ─────────────
        for _ in range(20):
            t += 1.0
            ts.append(t); ps.append(pf)
            fases.append(f"FASE B: FLUSH #{ciclo+1} VACÍO ✓")

    # Vacío final estable después de los 3 ciclos
    for _ in range(30):
        t += 1.0
        ts.append(t); ps.append(pf)
        fases.append("FASE B: VACÍO FINAL")
    
    return np.array(ts), np.array(ps), tau, fases

def sim_llenado(p0, pf, tasa):
    """
    Simulación de llenado según procedimiento:
    
    FASE C: LLENADO FINAL
    - Abrir válvula de gas MUY LENTAMENTE
    - Llenar hasta presión calculada (ej. 14-20 Torr)
    
    FASE D: PRUEBA Y ESTABILIZACIÓN
    - Esperar 10 minutos para verificar estabilidad
    - Comprobar que no hay micro-fugas
    
    Retorna arrays con fases marcadas
    """
    ts, ps, fases = [], [], []
    t = 0.0
    
    # ═════════════════════════════════════════════════════════════════════
    # FASE C: LLENADO GRADUAL (lento y controlado)
    # ═════════════════════════════════════════════════════════════════════
    while True:
        p_actual = p0 + tasa * t
        if p_actual >= pf or t >= 3600:
            break
        ts.append(t)
        ps.append(p_actual)
        fases.append("FASE C: LLENADO")
        t += 1.0
    
    # Asegurar que llegamos a presión objetivo
    ts.append(t)
    ps.append(pf)
    fases.append("FASE C: LLENADO")
    
    # ═════════════════════════════════════════════════════════════════════
    # FASE D: ESTABILIZACIÓN Y PRUEBA (10 minutos = 600 segundos)
    # ═════════════════════════════════════════════════════════════════════
    # Según procedimiento: "Esperar 10 minutos para que la presión se 
    # estabilice y verificar que no haya micro-fugas"
    
    for i in range(600):  # 10 minutos
        t += 1.0
        ts.append(t)
        ps.append(pf)  # Presión constante (sistema sellado)
        fases.append("FASE D: ESTABILIZACIÓN")
    
    return np.array(ts), np.array(ps), fases

def sim_vi(Vig, Vop, Imax):
    I = np.linspace(0, Imax, 300)
    V = np.where(I<0.1, Vig,
        np.where(I<5.0, Vig-(Vig-Vop)*(I/5.0), Vop))
    return I, V

# ── Gráficas matplotlib ────────────────────────────────────────────────────────

def render_latex(formula: str, fontsize: int=15, dark: bool=True) -> str:
    """Renderiza LaTeX como PNG y devuelve la ruta."""
    try:
        fig = plt.figure(figsize=(0.01, 0.01))
        color = "#f1f5f9" if dark else "#0a0f1a"
        bg    = "#162032" if dark else "#dde3ea"
        fig.patch.set_facecolor(bg)
        txt = fig.text(0.5, 0.5, f"${formula}$",
                       fontsize=fontsize, ha="center", va="center", color=color)
        fig.canvas.draw()
        bbox = txt.get_window_extent(renderer=fig.canvas.get_renderer())
        pad  = 22
        w = (bbox.width  + pad*2)/fig.dpi
        h = (bbox.height + pad*2)/fig.dpi
        fig.set_size_inches(max(w,1), max(h,0.45))
        path = os.path.join(TEMP_DIR, f"ltx_{abs(hash(formula+str(dark)))}.png")
        fig.savefig(path, dpi=fig.dpi, bbox_inches="tight", facecolor=bg)
        plt.close(fig)
        return path
    except Exception:
        return ""

def graficar_paschen(A, B, gamma, p_op=None, d_op=None, dark=True) -> str:
    umbral   = math.log(1.0+1.0/gamma)
    pd_ini   = umbral/A*1.001
    pdm      = paschen_pd_min(A, gamma)
    Vbm      = paschen_Vb_min(A, B, gamma)
    # Rango log amplio que siempre incluye el punto operativo
    pd_op_v = Vb_op_v = None
    if p_op and d_op and p_op > 0 and d_op > 0:
        pd_op_v = p_op * d_op
        v = paschen_Vb(pd_op_v, A, B, gamma)
        if math.isfinite(v) and v > 0: Vb_op_v = v
    pd_lo_log = max(pd_ini * 0.8, pd_ini)
    pd_hi_log = max(pdm * 2000, pd_ini * 1000)
    if pd_op_v: pd_hi_log = max(pd_hi_log, pd_op_v * 5)
    pd_log  = np.logspace(math.log10(pd_lo_log), math.log10(pd_hi_log), 8000)
    Vb_log  = np.array([paschen_Vb(pd, A, B, gamma) for pd in pd_log])
    # Zoom lineal: centrado en pdm, siempre incluye punto operativo
    pd_lo_lin = pdm * 0.05
    pd_hi_lin = pdm * 20
    if pd_op_v:
        pd_lo_lin = min(pd_lo_lin, pd_op_v * 0.2)
        pd_hi_lin = max(pd_hi_lin, pd_op_v * 1.8)
    pd_lin  = np.linspace(pd_lo_lin, pd_hi_lin, 6000)
    Vb_lin  = np.array([paschen_Vb(pd, A, B, gamma) for pd in pd_lin])
    # Recortar Vb_lin para no mostrar asíntotas
    Vb_cap  = max(Vbm * 10, Vb_op_v * 3 if Vb_op_v else Vbm * 10)
    Vb_lin_c= np.where((Vb_lin > Vb_cap) | (Vb_lin <= 0), np.nan, Vb_lin)

    bg,fg,gc,tc,cv = ("#0d1117","#e2e8f0","#1e293b","#64748b","#38bdf8") if dark else \
                     ("#f8fafc","#0f172a","#cbd5e1","#374151","#0369a1")
    plt.rcParams.update({'font.size': 11})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 8), dpi=110)
    fig.patch.set_facecolor(bg)

    for ax, pd_d, Vb_d, xlog, ttl in [
            (ax1, pd_log, Vb_log, True,  "Curva completa (escala log)"),
            (ax2, pd_lin, Vb_lin_c, False, "Zoom zona útil (lineal)")]:
        ax.set_facecolor(bg)
        ok = np.isfinite(Vb_d) & (Vb_d > 0)
        ax.plot(pd_d[ok], Vb_d[ok], color=cv, lw=2.4, label="Vb(p·d)", zorder=4)
        if not xlog: ax.fill_between(pd_d[ok], Vb_d[ok], alpha=0.08, color=cv)
        ax.axvline(pdm, color="#f97316", ls="--", lw=1.4, alpha=0.85, zorder=4)
        ax.axhline(Vbm, color="#f97316", ls="--", lw=1.4, alpha=0.85, zorder=4)
        ax.scatter([pdm], [Vbm], color="#f97316", s=100, zorder=7,
                   label=f"Mín: pd={pdm:.4g}, Vb={Vbm:.4g} V")
        if pd_op_v and Vb_op_v:
            # Asegurar que el punto está dentro del rango del eje antes de graficarlo
            if xlog:
                in_range = pd_lo_log <= pd_op_v <= pd_hi_log
            else:
                # Para lineal, el rango ya fue ajustado para incluirlo
                in_range = True
            if in_range:
                # Verificar que el punto realmente cae SOBRE la curva trazada
                ax.scatter([pd_op_v], [Vb_op_v], color="#fbbf24", s=140, zorder=9,
                           marker="*", label=f"Op: pd={pd_op_v:.4g}, Vb={Vb_op_v:.4g} V")
                ax.axvline(pd_op_v, color="#fbbf24", ls=":", lw=1.4, alpha=0.7, zorder=5)
                ax.axhline(Vb_op_v, color="#fbbf24", ls=":", lw=1.4, alpha=0.7, zorder=5)
        if xlog:
            ax.set_xscale("log")
            # Ajustar ylim en log para no mostrar valores impossibles
            ok2 = ok & (Vb_d > 0)
            if ok2.any():
                ylo = max(Vbm * 0.5, np.nanmin(Vb_d[ok2]))
                yhi = min(Vbm * 100, np.nanmax(Vb_d[ok2]))
                ax.set_ylim(ylo * 0.9, yhi * 1.1)
        ax.set_xlabel("p·d  [Torr·cm]", color=fg, fontsize=12)
        ax.set_ylabel("Vb  [V]", color=fg, fontsize=12)
        ax.set_title(ttl, color=fg, fontsize=12, fontweight="bold", pad=8)
        ax.tick_params(colors=tc, labelsize=10)
        for sp in ax.spines.values(): sp.set_edgecolor(gc)
        ax.grid(True, color=gc, lw=0.7, ls="--", alpha=0.9)
        ax.legend(facecolor=bg, edgecolor=gc, labelcolor=fg, fontsize=9, loc="upper right")
    fig.suptitle("Curva de Paschen — Vb vs p·d", color=fg, fontsize=14,
                 fontweight="bold", y=1.01)
    fig.text(0.5, -0.03,
             f"A={A:.4g}  B={B:.4g}  γ={gamma:.4g}  pd_min={pdm:.4g} Torr·cm  Vb_min={Vbm:.4g} V",
             ha="center", color=tc, fontsize=9)
    plt.tight_layout(pad=1.8)
    p = os.path.join(TEMP_DIR, f"paschen_curve_{int(time.time()*1000)}.png")
    fig.savefig(p, dpi=110, bbox_inches="tight", facecolor=bg); plt.close(fig)
    return p

def graficar_sim(t_ev,p_ev,t_ll,p_ll,A,B,cfg,vol_total,dark=True) -> str:
    BG,AX,TX = ("#0d1117","#1c2333","#e6edf3") if dark else ("#f8fafc","#e2e8f0","#0f172a")
    CE,CL,CO  = "#58a6ff","#3fb950","#f78166"
    CV,CP,CVI = "#ff7b72","#d2a8ff","#ffa657"
    gc = '#30363d' if dark else '#94a3b8'
    gamma = cfg.get('gamma',0.01); d = cfg['dist_electrodos']
    p_vac = cfg['p_vac']; p_obj = cfg['p_obj']

    # ── Limpiar arrays para escala log (evitar líneas verticales imposibles) ──
    p_ev_c = np.clip(np.array(p_ev, dtype=float), p_vac * 0.5, cfg.get('p_atm', 760) * 1.5)
    p_ll_c = np.clip(np.array(p_ll, dtype=float), p_vac * 0.5, cfg.get('p_atm', 760) * 1.5)
    t_ev_c = np.array(t_ev, dtype=float)
    t_ll_c = np.array(t_ll, dtype=float)
    # Eliminar puntos que saltarían verticalmente (delta_p / delta_t demasiado grande)
    def _smooth_log(t_arr, p_arr):
        if len(t_arr) < 2: return t_arr, p_arr
        dt = np.diff(t_arr); dt = np.where(dt == 0, 1e-9, dt)
        dp_ratio = np.abs(np.diff(np.log(np.maximum(p_arr, 1e-10)))) / dt
        mask = np.concatenate([[True], dp_ratio < 5.0])  # max 5 décadas/s
        return t_arr[mask], p_arr[mask]
    t_ev_s, p_ev_s = _smooth_log(t_ev_c, p_ev_c)
    t_ll_s, p_ll_s = _smooth_log(t_ll_c, p_ll_c)

    fig = plt.figure(figsize=(19, 13), dpi=100)
    fig.patch.set_facecolor(BG)
    gs_ = gridspec.GridSpec(2, 3, figure=fig, hspace=0.46, wspace=0.34,
                            width_ratios=[1.4, 1.4, 0.7],
                            left=0.07, right=0.97, top=0.94, bottom=0.07)

    def est(ax, t, xl, yl):
        ax.set_facecolor(AX); ax.set_title(t, color=TX, fontsize=10, fontweight='bold', pad=8)
        ax.set_xlabel(xl, color=TX, fontsize=9); ax.set_ylabel(yl, color=TX, fontsize=9)
        ax.tick_params(colors=TX, labelsize=8)
        for sp in ax.spines.values(): sp.set_edgecolor(gc)
        ax.grid(True, ls='--', alpha=0.22, color='#8b949e' if dark else '#94a3b8')

    # ── Presión vs tiempo (toda la fila superior) ─────────────────────────────
    ax1 = fig.add_subplot(gs_[0, :])
    off = float(t_ev_s[-1]) if len(t_ev_s) else 0
    ax1.plot(t_ev_s, p_ev_s, color=CE, lw=2.5, label='Fases A+B: Evacuación + Flushing')
    ax1.plot(off + t_ll_s, p_ll_s, color=CL, lw=2.5, label='Fases C+D: Llenado + Estabilización')
    ax1.axhline(p_obj, color=CO, ls='--', lw=1.6, label=f"P objetivo = {p_obj:.2f} Torr")
    ax1.axhline(p_vac, color=CV, ls=':', lw=1.2, label=f"Vacío = {p_vac:.2e} Torr")
    ax1.axvline(off, color='#f0b429', ls=':', lw=1.4, alpha=0.8, label='Inicio Fase C (llenado)')
    ax1.set_yscale('log')
    # ylim seguro para log
    all_p = np.concatenate([p_ev_s, p_ll_s])
    valid_p = all_p[(all_p > 0) & np.isfinite(all_p)]
    if len(valid_p):
        ax1.set_ylim(max(p_vac * 0.3, valid_p.min() * 0.5),
                     min(cfg.get('p_atm', 760) * 3, valid_p.max() * 3))
    est(ax1, 'PRESIÓN vs TIEMPO — Simulación Completa de Re-Gassing',
        'Tiempo (s)', 'Presión (Torr) [escala log]')
    ax1.legend(fontsize=9, facecolor=BG, labelcolor=TX, framealpha=0.92, loc='upper right')

    # ── Curva de Paschen ──────────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs_[1, 0])
    # Calcular Vb en el punto operativo
    pd_op = p_obj * d
    Vbop  = paschen_Vb(pd_op, A, B, gamma)
    # Rango que siempre incluye el punto operativo
    pd_lo = max(0.01, paschen_pd_min(A, gamma) * 0.05)
    pd_hi = max(1000, pd_op * 5)
    pd_a  = np.logspace(math.log10(pd_lo), math.log10(pd_hi), 2000)
    Vb_a  = np.array([paschen_Vb(pd, A, B, gamma) for pd in pd_a])
    ok    = np.isfinite(Vb_a) & (Vb_a > 0)
    ax2.plot(pd_a[ok] / d, Vb_a[ok] / 1000, color=CP, lw=2.2, label='Curva Paschen', zorder=4)
    if math.isfinite(Vbop) and Vbop > 0:
        ax2.scatter([p_obj], [Vbop / 1000], color=CO, s=100, zorder=7,
                    label=f"Op: {p_obj:.1f} Torr → {Vbop/1000:.3f} kV")
        ax2.axvline(p_obj, color=CO, ls=':', lw=1.2, alpha=0.7)
        ax2.axhline(Vbop / 1000, color=CO, ls=':', lw=1.2, alpha=0.7)
    ax2.axhline(cfg['v_ign'], color='#f0b429', ls='--', lw=1.4,
                label=f"V fuente = {cfg['v_ign']} kV")
    ax2.set_xscale('log')
    # Ajustar ylim para que el punto siempre sea visible
    if math.isfinite(Vbop) and Vbop > 0:
        vb_vis = Vb_a[ok] / 1000
        ylo = min(vb_vis[vb_vis > 0].min() * 0.8, Vbop / 1000 * 0.8)
        yhi = max(vb_vis.max() * 1.1,  Vbop / 1000 * 1.3, cfg['v_ign'] * 1.1)
        yhi = min(yhi, Vbop / 1000 * 15)  # no dejar el punto perdido abajo
        ax2.set_ylim(max(0, ylo), yhi)
    est(ax2, 'Ley de Paschen — V ruptura vs Presión', 'Presión (Torr)', 'Vb (kV)')
    ax2.legend(fontsize=8, facecolor=BG, labelcolor=TX)

    # ── Característica V-I ────────────────────────────────────────────────────
    ax3 = fig.add_subplot(gs_[1, 1])
    Ia, Va = sim_vi(cfg['v_ign'], cfg['v_op'], cfg['i_max'])
    ax3.plot(Ia, Va, color=CVI, lw=2.5, label='V-I plasma CO₂/He')
    ax3.axhline(cfg['v_ign'], color='#f0b429', ls='--', lw=1.2, alpha=0.8,
                label=f"V ign = {cfg['v_ign']} kV")
    ax3.axhline(cfg['v_op'], color=CL, ls='--', lw=1.2, alpha=0.8,
                label=f"V op = {cfg['v_op']} kV")
    ax3.annotate('Resistencia\nnegativa\ndel plasma',
                 xy=(3.0, (cfg['v_ign'] + cfg['v_op']) / 2),
                 xytext=(cfg['i_max'] * 0.45, cfg['v_ign'] * 0.93),
                 arrowprops=dict(arrowstyle='->', color=TX, lw=1),
                 fontsize=7.5, color=TX,
                 bbox=dict(boxstyle='round,pad=0.3', facecolor=AX, alpha=0.85))
    est(ax3, 'Característica V-I del Plasma', 'Corriente (mA)', 'Voltaje (kV)')
    ax3.set_xlim(0, cfg['i_max'])
    ax3.legend(fontsize=8, facecolor=BG, labelcolor=TX)

    # ── Pastel de mezcla (compacto) ───────────────────────────────────────────
    ax4 = fig.add_subplot(gs_[1, 2])
    ax4.set_facecolor(AX)
    mp = cfg.get('mezcla', {})
    if mp:
        gl = list(mp.keys()); fl = [mp[g] for g in gl]
        wedges, texts, autotexts = ax4.pie(
            fl, labels=gl, autopct='%1.1f%%', startangle=90,
            colors=[CE, CL, CP, CVI][:len(gl)],
            textprops={'color': TX, 'fontsize': 8},
            pctdistance=0.76, radius=0.78)
        for at in autotexts: at.set_fontsize(7)
    ax4.set_title('Composición\nMezcla', color=TX, fontsize=9, fontweight='bold', pad=4)

    fig.suptitle(
        f"SIMULADOR RE-GASSING — RECI  |  P = {p_obj:.2f} Torr  |  V = {vol_total/1000:.3f} L",
        color=TX, fontsize=13, fontweight='bold', y=0.998)
    plt.tight_layout(pad=1.5)
    p = os.path.join(TEMP_DIR, f"sim_graficas_{int(time.time()*1000)}.png")
    fig.savefig(p, dpi=100, bbox_inches='tight', facecolor=BG); plt.close(fig)
    return p

# ── Exportar reportes ──────────────────────────────────────────────────────────

def exportar(tipo: str, datos: dict) -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    ruta = os.path.join(EXPORT_DIR, f"reporte_paschen_{ts}.txt")
    sep  = "="*70
    with open(ruta,"w",encoding="utf-8") as f:
        f.write(sep+"\nREPORTE — PASCHEN v12 + SIMULADOR RE-GASSING\n"+sep+"\n")
        f.write(f"Fecha: {datetime.now():%Y-%m-%d %H:%M:%S}\nTipo: {tipo}\n"+sep+"\n\n")
        td = todos_los_gases()
        if tipo=="presion":
            d_=datos
            f.write(f"Vb={d_['Vb']} V  A={d_['A']}  B={d_['B']}  d={d_['d']} cm  γ={d_['gamma']}\n\n")
            pdm=paschen_pd_min(d_['A'],d_['gamma']); Vbm=paschen_Vb_min(d_['A'],d_['B'],d_['gamma'])
            f.write(f"Mínimo Paschen: pd_min={pdm:.6f}  Vb_min={Vbm:.6f} V\n\n")
            f.write(f"Presión = {d_['p_torr']:.6f} Torr\n")
            for u,v in convertir_presion(d_['p_torr'],'Torr').items():
                f.write(f"  {u:10s} = {v:.6f}\n")
        elif tipo=="mezcla":
            for g,pct in datos['pct'].items():
                f.write(f"  {g}: {pct:.2f}%  A={td[g]['A']}  B={td[g]['B']}\n")
            f.write(f"\nA_mezcla={datos['Am']:.6f}  B_mezcla={datos['Bm']:.6f}\n")
        elif tipo=="simulacion":
            d_=datos
            f.write(f"Volumen tubo:     {d_['vt']:.2f} cm³\n")
            f.write(f"Volumen manifold: {d_['vm']:.2f} cm³\n")
            f.write(f"Volumen total:    {d_['vtot']:.2f} cm³  ({d_['vtot']/1000:.4f} L)\n\n")
            f.write(f"P objetivo:       {d_['pobj']:.4f} Torr\n")
            f.write(f"Vb calculado:     {d_['vb']:.2f} V  ({d_['vb']/1000:.4f} kV)\n")
            f.write(f"Masa gas:         {d_['masa']:.4f} mg\n\n")
            f.write("Composición:\n")
            for g,pct in d_['mezcla'].items():
                f.write(f"  {g}: {pct:.1f}%\n")
        f.write("\n"+sep+"\nFIN\n"+sep+"\n")
    return ruta

def exp_img(src: str, prefijo: str="grafica") -> str:
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    dest = os.path.join(EXPORT_DIR,f"{prefijo}_{ts}.png")
    shutil.copy(src,dest); return dest

# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN 2 — TEMAS
# ═══════════════════════════════════════════════════════════════════════════════

DARK = {
    "bg":"#090e18","surface":"#0f172a","surface2":"#162032",
    "border":"#1e3a5f","border_f":"#38bdf8","text":"#f1f5f9",
    "muted":"#94a3b8","accent":"#38bdf8","label":"#7dd3fc",
    "green":"#4ade80","yellow":"#fde047","red":"#fb7185",
    "b_blue":"#0284c7","b_green":"#15803d","b_red":"#9f1239",
    "b_violet":"#7c3aed","b_gray":"#334155","b_orange":"#c2410c",
    "nav_sel":"#162032","nav_bg":"#0f172a","is_dark":True,
}
LIGHT = {
    "bg":"#f0f4f8","surface":"#ffffff","surface2":"#dde3ea",
    "border":"#64748b","border_f":"#0369a1","text":"#0a0f1a",
    "muted":"#334155","accent":"#0369a1","label":"#0c4a6e",
    "green":"#14532d","yellow":"#78350f","red":"#9f1239",
    "b_blue":"#0369a1","b_green":"#14532d","b_red":"#9f1239",
    "b_violet":"#5b21b6","b_gray":"#374151","b_orange":"#c2410c",
    "nav_sel":"#bfdbfe","nav_bg":"#dde3ea","is_dark":False,
}
C: dict = dict(DARK)

def set_theme(dark: bool):
    C.clear(); C.update(DARK if dark else LIGHT)

# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN 3 — COMPONENTES UI
# ═══════════════════════════════════════════════════════════════════════════════

def tf(label, hint="", value="", width=210):
    return ft.TextField(
        label=label, hint_text=hint, value=value,
        width=width, height=58, text_size=14,
        label_style=ft.TextStyle(color=C["label"],size=12),
        hint_style=ft.TextStyle(color=C["muted"],size=11),
        text_style=ft.TextStyle(color=C["text"],size=14),
        border_color=C["border"], focused_border_color=C["border_f"],
        cursor_color=C["accent"], bgcolor=C["surface"], border_radius=10,
    )

def rtf(f):
    f.label_style=ft.TextStyle(color=C["label"],size=12)
    f.hint_style=ft.TextStyle(color=C["muted"],size=11)
    f.text_style=ft.TextStyle(color=C["text"],size=14)
    f.border_color=C["border"]; f.focused_border_color=C["border_f"]
    f.cursor_color=C["accent"]; f.bgcolor=C["surface"]

def btn(txt, handler, color=None, icon=None, width=None):
    """Botón compatible Flet 0.80+ — Container clickeable."""
    c   = color or C["b_blue"]
    row = [ft.Text(txt, color="#ffffff", size=13,
                   weight=ft.FontWeight.W_600,
                   text_align=ft.TextAlign.CENTER)]
    if icon:
        row = [ft.Text(icon, size=15), *row]
    kwargs = dict(
        content=ft.Row(row, alignment=ft.MainAxisAlignment.CENTER,
                       spacing=6, tight=True),
        bgcolor=c, border_radius=10,
        padding=ft.Padding.symmetric(horizontal=18, vertical=11),
        on_click=handler, ink=True,
    )
    if width:
        kwargs['width'] = width
    return ft.Container(**kwargs)

def titulo(txt):
    return ft.Container(
        content=ft.Text(txt, size=19, weight=ft.FontWeight.BOLD, color=C["accent"]),
        padding=ft.Padding.only(bottom=8, top=4))

def sep():
    return ft.Divider(color=C["border"], height=1, thickness=1)

def card(items, pad=18):
    return ft.Container(
        content=ft.Column(items, spacing=10),
        bgcolor=C["surface"], border_radius=14, padding=pad,
        border=ft.Border.all(1, C["border"]))

def rrow(label, valor, color=None):
    return ft.Row([
        ft.Text(label, color=C["muted"], size=13, width=215),
        ft.Text(valor, color=color or C["text"], size=14,
                weight=ft.FontWeight.W_600, expand=True),
    ])

def nota(txt, color=None):
    return ft.Text(txt, color=color or C["yellow"], size=12, italic=True)

def banner_ok(msg):
    return ft.Container(
        content=ft.Row([ft.Text("✅",size=18),
                        ft.Text(msg,color="#fff",size=14,
                                weight=ft.FontWeight.W_600,expand=True)],spacing=10),
        bgcolor=C["b_green"], border_radius=10,
        padding=ft.Padding.symmetric(horizontal=16, vertical=10))

def banner_warn(msg):
    return ft.Container(
        content=ft.Row([ft.Text("⚠️",size=18),
                        ft.Text(msg,color="#fff",size=14,
                                weight=ft.FontWeight.W_600,expand=True)],spacing=10),
        bgcolor=C["b_orange"], border_radius=10,
        padding=ft.Padding.symmetric(horizontal=16, vertical=10))

def latex_img(formula, fontsize=15, dark=None):
    """Widget imagen de fórmula LaTeX; si falla, texto monospace."""
    d   = C["is_dark"] if dark is None else dark
    pth = render_latex(formula, fontsize, d)
    if pth and os.path.exists(pth):
        # Compatibilidad con diferentes versiones de Flet
        try:
            return ft.Image(src=pth, width=700, fit=ft.ImageFit.CONTAIN)
        except AttributeError:
            # Versiones antiguas de Flet no tienen ImageFit, usar string
            return ft.Image(src=pth, width=700)
    return ft.Text(f"[{formula}]", color=C["yellow"], size=12, font_family="Courier New")

def chip_gas(g, pct, dat, on_remove):
    return ft.Container(
        content=ft.Row([
            ft.Text(g, color=C["accent"], size=14, weight=ft.FontWeight.BOLD),
            ft.Text(f"{pct:.1f}%   A={dat['A']}  B={dat['B']}",
                    color=C["muted"], size=13),
            ft.Container(
                content=ft.Text("✕", color=C["red"], size=14,
                                weight=ft.FontWeight.BOLD),
                on_click=on_remove, padding=ft.Padding.symmetric(horizontal=4,vertical=2),
                border_radius=4, ink=True),
        ], spacing=8),
        bgcolor=C["surface2"], border_radius=8,
        padding=ft.Padding.symmetric(horizontal=10, vertical=5),
        border=ft.Border.all(1, C["border"]))

def ind_card(titulo_txt, lbl_ref, color_borde, sub="", ancho=150):
    """Tarjeta métrica SCADA que usa objeto ft.Text live."""
    return ft.Container(width=ancho,
        content=ft.Column([
            ft.Text(titulo_txt, color=C["muted"], size=10,
                    weight=ft.FontWeight.W_600, text_align=ft.TextAlign.CENTER),
            lbl_ref,
            ft.Text(sub, color=C["muted"], size=9,
                    text_align=ft.TextAlign.CENTER, italic=True),
        ], spacing=2, horizontal_alignment=ft.CrossAxisAlignment.CENTER),
        bgcolor=C["surface2"], border_radius=10, padding=10,
        border=ft.Border.all(1, color_borde))

# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN 4 — VISTAS
# ═══════════════════════════════════════════════════════════════════════════════

# ── ConversionView ─────────────────────────────────────────────────────────────
class ConversionView:
    def __init__(self, snack):
        self.snack  = snack
        self.tf_val = tf("Valor","ej. 760",width=190)
        self.tf_uni = tf("Unidad","Pa, Torr, atm…",width=170)
        self.res    = ft.Column([], spacing=5)
    def refresh_theme(self):
        rtf(self.tf_val)
        rtf(self.tf_uni)
        
        # Actualizar colores de resultados
        for control in self.res.controls:
            if isinstance(control, ft.Text):
                control.color = C["text"]
            elif isinstance(control, ft.Row):
                for item in control.controls:
                    if isinstance(item, ft.Text):
                        item.color = C["text"]
    def calcular(self, e):
        self.res.controls.clear()
        try:
            v = float(self.tf_val.value)
            c = convertir_presion(v, self.tf_uni.value.strip())
            self.res.controls.append(banner_ok(f"Conversión exitosa — {v} {self.tf_uni.value.strip()}"))
            for lbl,k in [("Pascal (Pa):","Pa"),("Torr:","Torr"),("Atmósferas:","atm"),
                          ("Milibar:","mbar"),("PSI:","psi"),("mmHg:","mmHg"),
                          ("inHg:","inHg"),("kg/cm²:","kg/cm²")]:
                self.res.controls.append(rrow(lbl, f"{c[k]:.6g} {k}"))
            self.snack("✅ Conversión realizada", C["green"])
        except Exception as ex:
            self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    def exp_txt(self, e):
        try:
            v=float(self.tf_val.value)
            c=convertir_presion(v,self.tf_uni.value.strip())
            p=exportar("conversion",{'valor':v,'unidad':self.tf_uni.value.strip(),'conversiones':c})
            # reutilizamos el tipo "conversion" pero el exportar actual lo tiene como conversion
            # Ajuste inline para compatibilidad
            self.snack(f"📄 {p}", C["accent"])
        except Exception as ex:
            self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    def build(self):
        return ft.Column([
            titulo("🔄  Conversión de Unidades de Presión"),
            card([ft.Row([self.tf_val, self.tf_uni], spacing=12),
                  ft.Row([btn("Convertir", self.calcular),
                          btn("🔁 Repetir", self.calcular, color=C["b_blue"]),
                          btn("Exportar .txt", self.exp_txt, color=C["b_green"])], spacing=10)]),
            sep(),
            ft.Text("Resultados", color=C["muted"], size=13),
            card([self.res]),
        ], spacing=14)

# ── MezclaView ─────────────────────────────────────────────────────────────────
class MezclaView:
    def __init__(self, snack, on_paschen=None, on_sim=None):
        self.snack       = snack
        self.on_paschen  = on_paschen
        self.on_sim      = on_sim
        self.gases: dict[str,float] = {}
        self.chips  = ft.Column([], spacing=5)
        self.res    = ft.Column([], spacing=5)
        self.tf_gas = tf("Gas","ej. N2, He, air",width=130)
        self.tf_pct = tf("%","ej. 30",width=100)
        self.tf_cg  = tf("Nombre","ej. Ar",width=115)
        self.tf_cA  = tf("A [1/(Torr·cm)]","ej. 14",width=155)
        self.tf_cB  = tf("B [V/(Torr·cm)]","ej. 350",width=155)
        self.lbl_g  = ft.Text(self._lbl(), color=C["muted"], size=12, italic=True)
        self._A = self._B = None
    def _lbl(self): return "Disponibles: "+", ".join(sorted(todos_los_gases().keys()))
    def refresh_theme(self):
        for f in [self.tf_gas,self.tf_pct,self.tf_cg,self.tf_cA,self.tf_cB]: 
            rtf(f)
        self.lbl_g.color = C["muted"]
        
        # Actualizar colores de resultados
        for control in self.res.controls:
            if isinstance(control, ft.Text):
                control.color = C["text"]
            elif isinstance(control, ft.Row):
                for item in control.controls:
                    if isinstance(item, ft.Text):
                        item.color = C["text"]
    def agregar(self, e):
        g = norm_gas(self.tf_gas.value.strip())
        if g not in todos_los_gases():
            self.snack(f"❌ '{g}' no reconocido.", C["red"]); e.page.update(); return
        try: pct = float(self.tf_pct.value)
        except: self.snack("❌ Porcentaje inválido.", C["red"]); e.page.update(); return
        self.gases[g] = self.gases.get(g, 0.0)+pct
        self._rebuild(); self.tf_gas.value = self.tf_pct.value = ""
        self.snack(f"✅ {g} {pct:.1f}% añadido", C["green"]); e.page.update()
    def quitar(self, g, page):
        self.gases.pop(g,None); self._rebuild(); page.update()
    def _rebuild(self):
        td = todos_los_gases(); self.chips.controls.clear()
        for g,pct in self.gases.items():
            dat = td.get(g,{'A':'?','B':'?'}); _g = g
            self.chips.controls.append(
                chip_gas(_g, pct, dat,
                         on_remove=lambda e, gg=_g: self.quitar(gg, e.page)))
        suma = sum(self.gases.values())
        sc = C["green"] if math.isclose(suma,100,rel_tol=1e-2) else C["yellow"]
        self.chips.controls.append(ft.Text(f"Suma: {suma:.2f}%", color=sc, size=13, italic=True))
    def calcular(self, e):
        self.res.controls.clear()
        if not self.gases:
            self.snack("❌ Agrega al menos un gas.", C["red"]); e.page.update(); return
        try:
            Am,Bm,fr,_ = calcular_AB_mezcla(self.gases)
            self._A, self._B = Am, Bm
            td = todos_los_gases(); rows=[]
            for g,frac in fr.items():
                rows += [rrow(f"{g} — frac molar:", f"{frac:.6f}"),
                         rrow("  contrib. A:", f"{frac*td[g]['A']:.6f}"),
                         rrow("  contrib. B:", f"{frac*td[g]['B']:.6f}")]
            rows += [sep(),
                     rrow("A (mezcla):",f"{Am:.6f} 1/(Torr·cm)",C["accent"]),
                     rrow("B (mezcla):",f"{Bm:.6f} V/(Torr·cm)",C["accent"]),
                     nota("💡 Pulsa '→ Paschen ⚡' o '→ Simulador 🔬'")]
            self.res.controls = [banner_ok(f"A y B calculados — {len(self.gases)} gas(es)")]+rows
            self.snack("✅ A y B calculados", C["green"])
        except Exception as ex:
            self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    def transferir_p(self, e):
        if self._A is None:
            self.snack("⚠️ Primero calcula A y B.", C["yellow"]); e.page.update(); return
        if self.on_paschen: self.on_paschen(self._A, self._B)
        self.snack(f"✅ A={self._A:.4f}, B={self._B:.4f} → Paschen", C["green"]); e.page.update()
    def transferir_s(self, e):
        if self._A is None:
            self.snack("⚠️ Primero calcula A y B.", C["yellow"]); e.page.update(); return
        if self.on_sim: self.on_sim(self._A, self._B, dict(self.gases))
        self.snack("✅ Mezcla enviada al Simulador", C["green"]); e.page.update()
    def exp_txt(self, e):
        if not self.gases:
            self.snack("❌ Sin gases.", C["red"]); e.page.update(); return
        try:
            Am,Bm,fr,pcts = calcular_AB_mezcla(self.gases)
            p = exportar("mezcla",{'pct':pcts,'Am':Am,'Bm':Bm})
            self.snack(f"📄 {p}", C["accent"])
        except Exception as ex: self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    def agregar_custom(self, e):
        n = norm_gas(self.tf_cg.value.strip())
        if n in IONIZATION_DATA:
            self.snack(f"⚠️ '{n}' ya es predefinido.", C["yellow"]); e.page.update(); return
        try: A,B = float(self.tf_cA.value), float(self.tf_cB.value)
        except: self.snack("❌ A y B deben ser números.", C["red"]); e.page.update(); return
        GASES_PERSONALIZADOS[n] = {'A':A,'B':B}
        self.lbl_g.value = self._lbl()
        self.tf_cg.value = self.tf_cA.value = self.tf_cB.value = ""
        self.snack(f"✅ Gas '{n}' guardado", C["green"]); e.page.update()
    def limpiar(self, e):
        self.gases.clear(); self.chips.controls.clear()
        self.res.controls.clear(); self._A = self._B = None; e.page.update()
    def build(self):
        return ft.Column([
            titulo("⚗️  Mezcla de Gases — Constantes de Townsend"),
            self.lbl_g,
            card([ft.Text("Añadir gas",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                  ft.Row([self.tf_pct,self.tf_gas,btn("Agregar",self.agregar)],spacing=10),
                  ft.Text("↑ Escribe primero el % y luego el nombre del gas.",
                          color=C["muted"],size=11,italic=True)]),
            ft.Text("Composición actual:",color=C["muted"],size=13),
            card([self.chips]),
            ft.Row([btn("Calcular A y B", self.calcular),
                    btn("🔁 Repetir", self.calcular, color=C["b_blue"]),
                    btn("→ Paschen ⚡", self.transferir_p, color=C["b_violet"]),
                    btn("→ Simulador 🔬", self.transferir_s, color=C["b_orange"]),
                    btn("Exportar .txt", self.exp_txt, color=C["b_green"]),
                    btn("Limpiar", self.limpiar, color=C["b_red"])],
                   spacing=10, wrap=True),
            sep(),
            ft.Text("Resultados",color=C["muted"],size=13),
            card([self.res]),
            sep(),
            card([ft.Text("➕ Gas personalizado",color=C["label"],size=13,
                          weight=ft.FontWeight.W_600),
                  ft.Row([self.tf_cg,self.tf_cA,self.tf_cB,
                          btn("Guardar",self.agregar_custom,color=C["b_green"])],
                         spacing=10, wrap=True)]),
        ], spacing=14)

# ── PaschenView ────────────────────────────────────────────────────────────────
class PaschenView:
    def __init__(self, snack, on_sim=None):
        self.snack    = snack
        self.on_sim   = on_sim
        self.tf_Vb    = tf("Vb — Voltaje ruptura (V)","ej. 400",width=195)
        self.tf_A     = tf("A [1/(Torr·cm)]","ej. 12",width=175)
        self.tf_B     = tf("B [V/(Torr·cm)]","ej. 342",width=175)
        self.tf_d     = tf("d — distancia (cm)","ej. 125",width=175)
        self.tf_gamma = tf("γ — coef. Townsend","ej. 0.01",width=175)
        self.res      = ft.Column([], spacing=5)
        self.img      = ft.Image(src="/",visible=False,width=820)
        self._last: dict = {}
    def refresh_theme(self):
        for f in [self.tf_Vb, self.tf_A, self.tf_B, self.tf_d, self.tf_gamma]:
            rtf(f)
        # Recolorear resultados
        for control in self.res.controls:
            if isinstance(control, ft.Text):
                control.color = C["text"]
            elif isinstance(control, ft.Row):
                for item in control.controls:
                    if isinstance(item, ft.Text):
                        item.color = C["text"]
        # Regenerar gráfica con nuevo tema si estaba visible
        if self.img.visible and self._last:
            try:
                p = self._parse()
                path = graficar_paschen(p['A'], p['B'], p['gamma'],
                                        p_op=self._last.get('p_torr'),
                                        d_op=p.get('d'), dark=C["is_dark"])
                self.img.src = ""
                self.img.src = path
            except Exception:
                pass
    def set_AB(self, A, B):
        self.tf_A.value = f"{A:.6f}"; self.tf_B.value = f"{B:.6f}"
    def _parse(self):
        return {'Vb':float(self.tf_Vb.value),'A':float(self.tf_A.value),
                'B':float(self.tf_B.value),'d':float(self.tf_d.value),
                'gamma':float(self.tf_gamma.value)}
    def calcular(self, e):
        self.res.controls.clear(); self.img.visible = False
        try:
            p = self._parse()
            pr  = solve_for_p(p['Vb'],p['A'],p['B'],p['d'],p['gamma'])
            pdo = pr*p['d']
            pdm = paschen_pd_min(p['A'],p['gamma'])
            Vbm = paschen_Vb_min(p['A'],p['B'],p['gamma'])
            Vchk= paschen_Vb(pdo,p['A'],p['B'],p['gamma'])
            conv= convertir_presion(pr,'Torr')
            self._last = {**p,'p_torr':pr}
            pm   = pdm/p['d']
            rama = "derecha (alta presión)" if pr>pm else "izquierda (baja presión)"
            margen = abs(pr-pm)/pm*100
            self.res.controls += [
                banner_ok(f"Presión = {pr:.4f} Torr  |  Vb verificado = {Vchk:.2f} V"),
                rrow("Presión calculada:",f"{pr:.6f} Torr",C["accent"]),
                rrow("p·d operativo:",f"{pdo:.6f} Torr·cm"),
                rrow("Rama de operación:",rama),
                rrow("Margen al mínimo:",f"{margen:.1f}%"),
                sep(),
                ft.Text("Mínimo de Paschen (analítico):",color=C["label"],size=13,
                        weight=ft.FontWeight.W_600),
                rrow("  pd_min:",f"{pdm:.6f} Torr·cm"),
                rrow("  p_min:", f"{pm:.6f} Torr"),
                rrow("  Vb_min:",f"{Vbm:.4f} V  ({Vbm/1000:.4f} kV)"),
                sep(),
                ft.Text("Conversiones:",color=C["muted"],size=13),
                *[rrow(k+":",f"{conv[k]:.6g} {k}") for k in conv],
            ]
            self.snack("✅ Cálculo completado", C["green"])
        except Exception as ex:
            self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    def graficar(self, e):
        try: p = self._parse()
        except:
            try:
                p={'A':float(self.tf_A.value),'B':float(self.tf_B.value),
                   'gamma':float(self.tf_gamma.value),'d':float(self.tf_d.value or "1"),'Vb':None}
            except Exception as ex:
                self.snack(f"❌ {ex}", C["red"]); e.page.update(); return
        try:
            path = graficar_paschen(p['A'], p['B'], p['gamma'],
                                    p_op=self._last.get('p_torr'), d_op=p.get('d'),
                                    dark=C["is_dark"])
            self.img.src = ""         # forzar reset en Flet
            self.img.src = path
            self.img.visible = True
            self.snack("✅ Gráfica generada", C["green"])
        except Exception as ex: self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    def exp_txt(self, e):
        if not self._last:
            self.snack("⚠️ Primero calcula.", C["yellow"]); e.page.update(); return
        try: p=exportar("presion",self._last); self.snack(f"📄 {p}", C["accent"])
        except Exception as ex: self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    def exp_png(self, e):
        if not self.img.visible:
            self.snack("⚠️ Primero genera la gráfica.", C["yellow"]); e.page.update(); return
        try: p=exp_img(self.img.src,"paschen"); self.snack(f"🖼️ {p}", C["accent"])
        except Exception as ex: self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    
    def transferir_s(self, e):
        """Enviar A, B, gamma y presión calculada a Simulación."""
        try:
            p = self._parse()
            p_torr = self._last.get('p_torr') if self._last else None
            if self.on_sim:
                self.on_sim(p['A'], p['B'], p['gamma'], p_torr)
            extra = f"  P={p_torr:.4f} Torr → P objetivo" if p_torr else "  (sin presión calculada)"
            self.snack(f"✅ A={p['A']:.4f}, B={p['B']:.4f}, γ={p['gamma']:.4f}{extra}", C["green"])
        except Exception as ex:
            self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    
    def build(self):
        return ft.Column([
            titulo("⚡  Cálculo de Presión Óptima — Ley de Paschen"),
            card([ft.Text("Parámetros de entrada",color=C["label"],size=13,
                          weight=ft.FontWeight.W_600),
                  ft.Row([self.tf_Vb,self.tf_A,self.tf_B],spacing=12,wrap=True),
                  ft.Row([self.tf_d,self.tf_gamma],spacing=12),
                  nota("💡 Ve a ⚗️ Mezcla, calcula A y B y usa '→ Paschen ⚡' para llenar A y B."),
                  ft.Row([btn("🔄 Recalcular",self.calcular,color=C["b_blue"]),
                          btn("Graficar V vs p·d",self.graficar,color=C["b_gray"]),
                          btn("→ Simulador 🔬",self.transferir_s,color=C["b_orange"])],spacing=10,wrap=True),
                  sep(),
                  ft.Row([btn("Exportar .txt",self.exp_txt,color=C["b_green"]),
                          btn("Exportar .png",self.exp_png,color=C["b_green"])],spacing=10),
                  nota(f"📁 Archivos en: {EXPORT_DIR}",C["muted"])]),
            sep(),
            ft.Text("Resultados",color=C["muted"],size=13),
            card([self.res]),
            ft.Container(height=8),
            ft.Row([self.img],alignment=ft.MainAxisAlignment.CENTER),
        ], spacing=14)

# ── GasesView ──────────────────────────────────────────────────────────────────
class GasesView:
    def __init__(self, snack):
        self.snack = snack; self.tabla = ft.Column([], spacing=5)
        self.tf_n = tf("Nombre","ej. Ar",width=120)
        self.tf_A = tf("A [1/(Torr·cm)]","ej. 14",width=155)
        self.tf_B = tf("B [V/(Torr·cm)]","ej. 350",width=155)
    def refresh_theme(self):
        for f in [self.tf_n,self.tf_A,self.tf_B]: rtf(f)
    def refrescar(self, page):
        self.tabla.controls.clear()
        self.tabla.controls.append(ft.Row([
            ft.Text("Gas",color=C["label"],size=14,width=80,weight=ft.FontWeight.BOLD),
            ft.Text("A  [1/(Torr·cm)]",color=C["label"],size=14,width=145,weight=ft.FontWeight.BOLD),
            ft.Text("B  [V/(Torr·cm)]",color=C["label"],size=14,width=145,weight=ft.FontWeight.BOLD),
            ft.Text("Tipo",color=C["label"],size=14,width=110,weight=ft.FontWeight.BOLD),
        ])); self.tabla.controls.append(sep())
        for g,d in sorted(todos_los_gases().items()):
            tipo = "Predefinido" if g in IONIZATION_DATA else "Personalizado"
            self.tabla.controls.append(ft.Row([
                ft.Text(g,color=C["text"],size=14,width=80),
                ft.Text(str(d['A']),color=C["accent"],size=14,width=145),
                ft.Text(str(d['B']),color=C["accent"],size=14,width=145),
                ft.Text(tipo,color=C["muted"] if tipo=="Predefinido" else C["yellow"],
                        size=13,width=110),
            ]))
        if page: page.update()
    def guardar(self, e):
        n = norm_gas(self.tf_n.value.strip())
        if n in IONIZATION_DATA:
            self.snack(f"⚠️ '{n}' ya predefinido.", C["yellow"]); e.page.update(); return
        try: A,B = float(self.tf_A.value), float(self.tf_B.value)
        except: self.snack("❌ A y B numéricos.", C["red"]); e.page.update(); return
        GASES_PERSONALIZADOS[n] = {'A':A,'B':B}
        self.tf_n.value = self.tf_A.value = self.tf_B.value = ""
        self.refrescar(e.page); self.snack(f"✅ '{n}' guardado", C["green"]); e.page.update()
    def build(self):
        self.refrescar(None)
        return ft.Column([
            titulo("🧪  Gases y Constantes de Townsend"),
            card([self.tabla]),
            sep(),
            card([ft.Text("Agregar gas personalizado",color=C["label"],size=13,
                          weight=ft.FontWeight.W_600),
                  ft.Row([self.tf_n,self.tf_A,self.tf_B,
                          btn("Guardar",self.guardar,color=C["b_green"])],
                         spacing=10, wrap=True)]),
        ], spacing=14)

# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN 5 — SIMULADOR RE-GASSING (vista principal nueva)
# ═══════════════════════════════════════════════════════════════════════════════

class SimuladorView:
    """
    Simulador completo de re-gassing para tubo láser RECI.
    Pestañas internas: Configuración → Resultados → Simulación → Guías
    """

    # ── Sub-pestañas ──────────────────────────────────────────────────────────
    TAB_CONFIG  = 0
    TAB_RES     = 1
    TAB_SIM     = 2
    TAB_GUIAS   = 3

    def __init__(self, snack):
        self.snack = snack

        # ── Campos de configuración ───────────────────────────────────────────
        self.tf_long   = tf("Longitud tubo (cm)",   "140",   str(SIM_DEF['longitud_tubo']),   width=155)
        self.tf_diam   = tf("Diámetro interno (cm)","8",     str(SIM_DEF['diametro_tubo']),   width=155)
        self.tf_volm   = tf("Volumen manual (cm³)", "5500",  str(SIM_DEF['vol_tubo_manual']), width=155)
        self.chk_volm  = ft.Checkbox(label="Usar volumen manual",value=False,
                                     active_color=C["accent"],
                                     label_style=ft.TextStyle(color=C["muted"],size=12))
        self.tf_vcuerpo= tf("Vol. cuerpo manifold (cm³)","50",str(SIM_DEF['vol_cuerpo']),width=175)
        self.tf_lgg    = tf("L manguera gas (cm)",  "100",  str(SIM_DEF['lg_gas']),  width=150)
        self.tf_dgg    = tf("D manguera gas (cm)",  "0.635",str(SIM_DEF['dg_gas']),  width=150)
        self.tf_lgb    = tf("L manguera bomba (cm)","100",  str(SIM_DEF['lg_bomba']),width=150)
        self.tf_dgb    = tf("D manguera bomba (cm)","1.27", str(SIM_DEF['dg_bomba']),width=150)
        self.tf_lgl    = tf("L manguera láser (cm)","20",   str(SIM_DEF['lg_laser']),width=150)
        self.tf_dgl    = tf("D manguera láser (cm)","0.635",str(SIM_DEF['dg_laser']),width=150)
        self.tf_pobj   = tf("P objetivo (Torr)",    "14.4", str(SIM_DEF['p_obj']),   width=155)
        self.tf_patm   = tf("P atm local (Torr)",   "585",  str(SIM_DEF['p_atm']),   width=155)
        self.tf_pvac   = tf("P vacío obj. (Torr)",  "0.001",str(SIM_DEF['p_vac']),   width=155)
        self.tf_vign   = tf("V ignición (kV)",       "24",  str(SIM_DEF['v_ign']),   width=155)
        self.tf_vop    = tf("V operación (kV)",      "18",  str(SIM_DEF['v_op']),    width=155)
        self.tf_imax   = tf("I máx (mA)",             "28", str(SIM_DEF['i_max']),   width=155)
        self.tf_delec  = tf("d electrodos (cm)",     "125", str(SIM_DEF['dist_electrodos']),width=155)
        self.tf_tev    = tf("Tasa evacuación (Torr/s)","5", str(SIM_DEF['tasa_ev']),  width=165)
        self.tf_tll    = tf("Tasa llenado (Torr/s)", "0.05",str(SIM_DEF['tasa_ll']),  width=165)
        self.tf_gamma  = tf("γ (coef. Townsend 2°)", "0.01",str(SIM_DEF['gamma']),    width=155)
        self.tf_co2    = tf("CO₂ (%)", "8",  str(SIM_DEF['mezcla']['CO2']),width=110)
        self.tf_n2     = tf("N₂ (%)",  "15", str(SIM_DEF['mezcla']['N2']), width=110)
        self.tf_h2     = tf("H₂ (%)",  "5",  str(SIM_DEF['mezcla']['H2']), width=110)
        self.tf_he     = tf("He (%)",  "72", str(SIM_DEF['mezcla']['He']), width=110)

        # ── Estado interno ────────────────────────────────────────────────────
        self._cfg: dict = {}
        self._res_col   = ft.Column([], spacing=5)
        self._img_sim   = ft.Image(src="/", visible=False, width=1100)
        self._tab_idx   = [self.TAB_CONFIG]

        # Indicadores SCADA (objetos vivos)
        self._lbl_p    = ft.Text("—",     color=C["accent"], size=20, weight=ft.FontWeight.BOLD,
                                  font_family="Courier New", text_align=ft.TextAlign.CENTER)
        self._lbl_fase = ft.Text("STANDBY",color=C["yellow"], size=13, weight=ft.FontWeight.BOLD,
                                  text_align=ft.TextAlign.CENTER)
        self._lbl_t    = ft.Text("0 s",   color=C["muted"], size=13, text_align=ft.TextAlign.CENTER)
        self._lbl_masa = ft.Text("—",     color=C["green"], size=14, text_align=ft.TextAlign.CENTER)
        self._lbl_vb   = ft.Text("—",     color=C["accent"],size=14, text_align=ft.TextAlign.CENTER)
        self._lbl_vtot = ft.Text("—",     color=C["muted"], size=13, text_align=ft.TextAlign.CENTER)
        self._prog_bar = ft.ProgressBar(value=0, color=C["accent"],
                                         bgcolor=C["surface2"], width=480)
        self._prog_lbl = ft.Text("0 %", color=C["muted"], size=11)

        # Estado válvulas
        self._vs  = {'V1':False,'V2':False,'V3':False}
        self._vb  = {}   # btn containers
        # Estado equipos
        self._eq  = {'bomba':False,'fuente':False,'tubo':False,'tanque':False}
        self._eqb = {}

        # Control simulación
        self._sim_running  = False
        self._sim_speed    = [1.0]   # multiplicador de velocidad
        self._sim_paused   = [False]

        # Tab containers (se rellenan en build)
        self._tab_content  = ft.Column([], spacing=0)
        self._tab_btns: list[ft.Container] = []

    # ── Refresh ───────────────────────────────────────────────────────────────
    def refresh_theme(self):
        for f in [self.tf_long,self.tf_diam,self.tf_volm,self.tf_vcuerpo,
                  self.tf_lgg,self.tf_dgg,self.tf_lgb,self.tf_dgb,self.tf_lgl,self.tf_dgl,
                  self.tf_pobj,self.tf_patm,self.tf_pvac,self.tf_vign,self.tf_vop,
                  self.tf_imax,self.tf_delec,self.tf_tev,self.tf_tll,self.tf_gamma,
                  self.tf_co2,self.tf_n2,self.tf_h2,self.tf_he]:
            rtf(f)
        for control in self._res_col.controls:
            if isinstance(control, ft.Text):
                control.color = C["text"]
            elif isinstance(control, ft.Row):
                for item in control.controls:
                    if isinstance(item, ft.Text):
                        item.color = C["text"]
        self._lbl_p.color    = C["accent"]
        self._lbl_fase.color = C["yellow"]
        self._lbl_t.color    = C["muted"]
        self._lbl_masa.color = C["green"]
        self._lbl_vb.color   = C["accent"]
        self._lbl_vtot.color = C["muted"]
        self._prog_bar.color  = C["accent"]
        self._prog_bar.bgcolor= C["surface2"]
        self._prog_lbl.color  = C["muted"]
        # Regenerar gráfica estática con nuevo tema si estaba visible
        if self._img_sim.visible and 't_ev' in self._cfg:
            try:
                c = self._cfg
                path = graficar_sim(c['t_ev'], c['p_ev'], c['t_ll'], c['p_ll'],
                                    c['Am'], c['Bm'],
                                    {'p_obj': c['p_obj'], 'p_vac': c['p_vac'],
                                     'p_atm': c.get('p_atm', 760),
                                     'v_ign': c['v_ign'], 'v_op': c['v_op'],
                                     'i_max': c['i_max'],
                                     'dist_electrodos': c['dist_electrodos'],
                                     'gamma': c['gamma'], 'mezcla': c['mezcla']},
                                    c['vtot'], dark=C["is_dark"])
                self._img_sim.src = ""
                self._img_sim.src = path
            except Exception:
                pass

    # ── Inyección desde Mezcla ────────────────────────────────────────────────
    def recibir_mezcla(self, A, B, pct):
        for g,tf_ in [('CO2',self.tf_co2),('N2',self.tf_n2),
                       ('H2',self.tf_h2),('He',self.tf_he)]:
            if g in pct: tf_.value = f"{pct[g]:.2f}"
    
    def recibir_AB_gamma(self, A, B, gamma, p_torr=None):
        """Recibe A, B, gamma (y opcionalmente la presión calculada) desde Paschen."""
        self.tf_gamma.value = f"{gamma:.6f}"
        if p_torr is not None and p_torr > 0:
            self.tf_pobj.value = f"{p_torr:.4f}"

    # ── Parse parámetros ──────────────────────────────────────────────────────
    def _parse(self):
        if self.chk_volm.value:
            vt = float(self.tf_volm.value)
        else:
            l  = float(self.tf_long.value); d = float(self.tf_diam.value)
            vt = math.pi*(d/2)**2*l
        vm  = vol_manifold(float(self.tf_vcuerpo.value),
                           float(self.tf_lgg.value),  float(self.tf_dgg.value),
                           float(self.tf_lgb.value),  float(self.tf_dgb.value),
                           float(self.tf_lgl.value),  float(self.tf_dgl.value))
        vtot = vt+vm
        pct  = {'CO2':float(self.tf_co2.value),'N2':float(self.tf_n2.value),
                'H2':float(self.tf_h2.value),  'He':float(self.tf_he.value)}
        sp   = sum(pct.values())
        if abs(sp-100)>5: raise ValueError(f"Mezcla suma {sp:.1f}% — debe ser ~100%")
        Am,Bm,fr,_ = calcular_AB_mezcla(pct)
        return {
            'vt':vt,'vm':vm,'vtot':vtot,'mezcla':pct,'fracs':fr,
            'Am':Am,'Bm':Bm,'gamma':float(self.tf_gamma.value),
            'p_obj':float(self.tf_pobj.value), 'p_atm':float(self.tf_patm.value),
            'p_vac':float(self.tf_pvac.value), 'v_ign':float(self.tf_vign.value),
            'v_op':float(self.tf_vop.value),   'i_max':float(self.tf_imax.value),
            'dist_electrodos':float(self.tf_delec.value),
            'tasa_ev':float(self.tf_tev.value),'tasa_ll':float(self.tf_tll.value),
        }

    # ── Calcular ──────────────────────────────────────────────────────────────
    def calcular(self, e):
        self._res_col.controls.clear(); self._img_sim.visible = False
        try:
            c = self._parse(); self._cfg = c
            mg,n,pm,fr = masa_gas(c['vtot'],c['p_obj'],c['mezcla'])
            Vb = paschen_Vb(c['p_obj']*c['dist_electrodos'],c['Am'],c['Bm'],c['gamma'])
            t_ev,p_ev,tau,fases_ev = sim_evacuacion(c['p_atm'],c['p_vac'],c['tasa_ev'])
            t_ll,p_ll,fases_ll     = sim_llenado(c['p_vac'],c['p_obj'],c['tasa_ll'])
            pdm = paschen_pd_min(c['Am'],c['gamma'])
            Vbm = paschen_Vb_min(c['Am'],c['Bm'],c['gamma'])
            self._cfg.update({'mg':mg,'n':n,'pm':pm,'Vb':Vb,
                              't_ev':t_ev,'p_ev':p_ev,'t_ll':t_ll,'p_ll':p_ll,
                              'tau':tau,'pdm':pdm,'Vbm':Vbm,
                              'fases_ev':fases_ev,'fases_ll':fases_ll})
            ok_v = math.isfinite(Vb) and Vb/1000<=c['v_ign']
            vb_s = f"{Vb/1000:.4f} kV" if math.isfinite(Vb) else "N/A"
            col  = C["green"] if ok_v else C["red"]
            ok_s = "✓ Fuente suficiente" if ok_v else "⚠ Voltaje de fuente INSUFICIENTE"
            # Actualizar indicadores SCADA
            self._lbl_p.value    = f"{c['p_atm']:.1f} Torr"
            self._lbl_masa.value = f"{mg:.4f} mg"
            self._lbl_vb.value   = vb_s
            self._lbl_vtot.value = f"{c['vtot']:.1f} cm³"
            self._lbl_fase.value = "CONFIGURADO ✓"
            rows = [
                banner_ok("Sistema calculado correctamente") if ok_v else banner_warn("Revisar voltaje de fuente"),
                ft.Text("Volúmenes:",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                rrow("  Tubo RECI:",        f"{c['vt']:.2f} cm³  ({c['vt']/1000:.4f} L)"),
                rrow("  Manifold muerto:",  f"{c['vm']:.2f} cm³  ({c['vm']/1000:.4f} L)"),
                rrow("  TOTAL:",            f"{c['vtot']:.2f} cm³  ({c['vtot']/1000:.4f} L)",C["accent"]),
                sep(),
                ft.Text("Mezcla SparkLaser:",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                *[rrow(f"  {g}:",f"{c['mezcla'][g]:.1f}%  (frac={fr[g]:.4f})")for g in c['mezcla']],
                rrow("  PM mezcla:",f"{pm:.4f} g/mol"),
                rrow("  A (Townsend):",f"{c['Am']:.4f} 1/(Torr·cm)",C["accent"]),
                rrow("  B (Townsend):",f"{c['Bm']:.4f} V/(Torr·cm)",C["accent"]),
                sep(),
                ft.Text("Masa de gas — PV=nRT (R=62.364 L·Torr/mol·K, T=298 K):",
                        color=C["label"],size=13,weight=ft.FontWeight.W_600),
                rrow("  P objetivo:",   f"{c['p_obj']:.4f} Torr"),
                rrow("  V total:",      f"{c['vtot']/1000:.6f} L"),
                rrow("  n (moles):",    f"{n:.6e} mol"),
                rrow("  Masa total:",   f"{mg:.4f} mg  ({mg*1000:.2f} µg)",C["green"]),
                *[rrow(f"    {g}:",f"{fr[g]*mg:.4f} mg")for g in c['mezcla']],
                sep(),
                ft.Text("Ley de Paschen:",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                rrow("  pd_min (analítico):", f"{pdm:.6f} Torr·cm"),
                rrow("  p_min:",              f"{pdm/c['dist_electrodos']:.6f} Torr"),
                rrow("  Vb_min:",             f"{Vbm:.4f} V  ({Vbm/1000:.4f} kV)"),
                rrow("  pd operativo:",       f"{c['p_obj']*c['dist_electrodos']:.4f} Torr·cm"),
                rrow("  Vb calculado:",       vb_s, col),
                rrow("  Voltaje fuente:",     f"{c['v_ign']:.1f} kV"),
                ft.Text(ok_s,color=col,size=13,italic=True),
                sep(),
                ft.Text("Proceso simulado:",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                rrow("  τ evacuación:",    f"{tau:.1f} s"),
                rrow("  t evacuación:",   f"~{float(t_ev[-1]):.0f} s  ({float(t_ev[-1])/60:.1f} min)"),
                rrow("  t llenado:",      f"~{float(t_ll[-1]):.0f} s  ({float(t_ll[-1])/60:.1f} min)"),
                rrow("  Tasa llenado:",   f"{c['tasa_ll']:.4f} Torr/s"),
            ]
            self._res_col.controls += rows
            self.snack("✅ Sistema calculado — ir a pestaña Resultados o Simulación", C["green"])
            # Cambiar a pestaña resultados automáticamente
            self._switch_tab(self.TAB_RES, e.page)
        except Exception as ex:
            self.snack(f"❌ {ex}", C["red"])
        e.page.update()

    def graficar(self, e):
        if 't_ev' not in self._cfg:
            self.snack("⚠️ Primero calcula el sistema.", C["yellow"]); e.page.update(); return
        try:
            c = self._cfg
            path = graficar_sim(c['t_ev'], c['p_ev'], c['t_ll'], c['p_ll'],
                                 c['Am'], c['Bm'],
                                 {'p_obj': c['p_obj'], 'p_vac': c['p_vac'],
                                  'p_atm': c.get('p_atm', 760),
                                  'v_ign': c['v_ign'], 'v_op': c['v_op'],
                                  'i_max': c['i_max'],
                                  'dist_electrodos': c['dist_electrodos'],
                                  'gamma': c['gamma'], 'mezcla': c['mezcla']},
                                 c['vtot'], dark=C["is_dark"])
            self._img_sim.src = ""    # forzar reset en Flet
            self._img_sim.src = path
            self._img_sim.visible = True
            self.snack("✅ Gráficas generadas", C["green"])
        except Exception as ex: self.snack(f"❌ {ex}", C["red"])
        e.page.update()

    # ── Simulación en ventana matplotlib ──────────────────────────────────────
    def iniciar_sim(self, e):
        if not self._cfg.get('p_obj'):
            self.snack("⚠️ Primero calcula el sistema (pestaña Configuración).", C["yellow"])
            e.page.update(); return
        if self._sim_running:
            self.snack("⚠️ Simulación ya en curso.", C["yellow"]); e.page.update(); return
        cfg      = dict(self._cfg)
        dark     = C["is_dark"]
        page_ref = e.page
        self._sim_paused[0]  = False
        self._sim_speed[0]   = 1.0
        def _run():
            self._sim_running = True
            self._lbl_fase.value = "INICIANDO..."
            try: page_ref.update()
            except: pass
            try:
                _sim_ventana(cfg, dark, self._sim_speed, self._sim_paused,
                    on_upd=lambda f,p,t,pr: self._upd(f,p,t,pr,page_ref),
                    on_end=lambda: self._end(page_ref))
            except Exception:
                pass
            finally:
                self._sim_running = False
        threading.Thread(target=_run, daemon=True).start()
        self.snack("▶️ Simulación iniciada — ventana matplotlib", C["accent"])
        e.page.update()

    def set_speed(self, factor: float, e):
        self._sim_speed[0] = factor
        self.snack(f"⚡ Velocidad ×{factor:.0f}", C["accent"]); e.page.update()

    def toggle_pausa(self, e):
        self._sim_paused[0] = not self._sim_paused[0]
        estado = "⏸ PAUSADO" if self._sim_paused[0] else "▶ REANUDADO"
        self._lbl_fase.value = estado
        self.snack(estado, C["yellow"]); e.page.update()

    def _upd(self, fase, presion, t, prog, page):
        try:
            self._lbl_p.value    = f"{presion:.3f} Torr"
            self._lbl_fase.value = fase
            self._lbl_t.value    = f"{t:.0f} s"
            self._prog_bar.value = prog
            self._prog_lbl.value = f"{prog*100:.0f}%"
            page.update()
        except: pass

    def _end(self, page):
        try:
            self._lbl_fase.value = "COMPLETADO ✓"
            self._prog_bar.value = 1.0; self._prog_lbl.value = "100%"
            page.update()
        except: pass

    # ── Toggle válvulas ───────────────────────────────────────────────────────
    def _toggle_v(self, nombre, page):
        self._vs[nombre] = not self._vs[nombre]
        est = self._vs[nombre]
        b = self._vb.get(nombre)
        if b:
            b._lbl.value  = f"{nombre}  {'ABIERTA ●' if est else 'CERRADA ○'}"
            b.bgcolor     = C["b_green"] if est else C["b_red"]
        fases_ok = self._validar_secuencia()
        self.snack(f"{'🟢' if est else '🔴'} {nombre} → {'ABIERTA' if est else 'CERRADA'}  |  {fases_ok}",
                   C["green"] if est else C["red"])
        page.update()

    def _validar_secuencia(self) -> str:
        v1,v2,v3 = self._vs['V1'],self._vs['V2'],self._vs['V3']
        if not v1 and not v2 and not v3: return "Estado: REPOSO"
        if not v1 and v2 and not v3:     return "Estado: PURGANDO manifold"
        if v1 and v2 and not v3:          return "Estado: EVACUANDO sistema"
        if v1 and not v2 and not v3:      return "Estado: ESTANQUEIDAD / LISTO"
        if v1 and not v2 and v3:          return "Estado: LLENANDO gas ✅"
        if v3 and v2:                     return "⚠ PELIGRO: vacío + gas simultáneos"
        return "Estado: combinación no estándar"

    def _make_vbtn(self, nombre):
        lbl = ft.Text(f"{nombre}  CERRADA ○", color="#fff", size=12,
                      weight=ft.FontWeight.W_600, text_align=ft.TextAlign.CENTER)
        b = ft.Container(content=lbl, bgcolor=C["b_red"], border_radius=8,
                         padding=ft.Padding.symmetric(horizontal=14,vertical=10),
                         on_click=lambda e,n=nombre: self._toggle_v(n,e.page), ink=True)
        b._lbl = lbl; self._vb[nombre] = b; return b

    # ── Toggle equipos ────────────────────────────────────────────────────────
    def _toggle_eq(self, nombre, page):
        self._eq[nombre] = not self._eq[nombre]
        est = self._eq[nombre]
        etiqs = {
            'bomba': ('💨 Bomba  ON','💨 Bomba  OFF'),
            'fuente':('⚡ Fuente ON','⚡ Fuente OFF'),
            'tubo':  ('🔬 Tubo  MONTADO','🔬 Tubo  DESMONTADO'),
            'tanque':('🛢 Tanque  CONECTADO','🛢 Tanque  DESCONECTADO'),
        }
        eon,eoff = etiqs.get(nombre,(nombre+' ON',nombre+' OFF'))
        b = self._eqb.get(nombre)
        if b:
            b._lbl.value = eon if est else eoff
            b.bgcolor    = C["b_green"] if est else C["b_gray"]
        self.snack(f"{'🟢' if est else '⚫'} {eon if est else eoff}", C["green"])
        page.update()

    def _make_eqbtn(self, nombre, etiq_off):
        lbl = ft.Text(etiq_off, color="#fff", size=12,
                      weight=ft.FontWeight.W_600, text_align=ft.TextAlign.CENTER)
        b = ft.Container(content=lbl, bgcolor=C["b_gray"], border_radius=8,
                         padding=ft.Padding.symmetric(horizontal=12,vertical=10),
                         on_click=lambda e,n=nombre: self._toggle_eq(n,e.page), ink=True)
        b._lbl = lbl; self._eqb[nombre] = b; return b

    # ── Guías ─────────────────────────────────────────────────────────────────
    def _guia_dlg(self, tipo, e):
        if not self._cfg.get('p_obj'):
            self.snack("⚠️ Primero calcula el sistema.", C["yellow"]); e.page.update(); return
        c   = self._cfg
        tau = c.get('tau', c['p_atm']/c['tasa_ev'])
        t_ev= float(c['t_ev'][-1]) if 't_ev' in c else 0
        t_ll= float(c['t_ll'][-1]) if 't_ll' in c else 0
        po  = c['p_obj']; pv=c['p_vac']; tll=c['tasa_ll']
        mg  = c.get('mg',0)
        if tipo=="evacuacion":
            lineas = [
                "═══ GUÍA DE EVACUACIÓN — PASO A PASO ═══",
                f"Meta: ≤ {pv:.2e} Torr  |  τ ≈ {tau:.0f} s  |  t_total ≈ {t_ev:.0f} s ({t_ev/60:.1f} min)",
                "","IDENTIFICACIÓN DE VÁLVULAS:",
                "  V1 — Maestra:    manifold ↔ tubo láser",
                "  V2 — Vacío:      manifold ↔ bomba de vacío rotativa",
                "  V3 — Inyección:  manifold ↔ cilindro de mezcla de gas",
                "  M1 — Vacuómetro Pirani  (760→10⁻³ Torr)",
                "  M2 — Capacitance Manometer Baratron (0–100 Torr, ±0.01 Torr)",
                "","──── PASO 1 [t=0s]  VERIFICACIÓN INICIAL ────",
                f"  • V1 CERRADA  V2 CERRADA  V3 CERRADA",
                f"  • Verifica aceite de la bomba rotativa",
                f"  • M1 debe leer ~{c['p_atm']:.0f} Torr (presión atmosférica local)",
                "","──── PASO 2 [t≈0–30s]  PURGADO DEL MANIFOLD ────",
                "  • Enciende la bomba de vacío",
                "  • Abre V2  →  V1 CERRADA  V2 ABIERTA  V3 CERRADA",
                f"  • M1: {c['p_atm']:.0f} → ~10 Torr en ~30 s (solo volumen del manifold)",
                "","──── PASO 3 [t≈30–60s]  CONECTAR EL TUBO ────",
                "  • Abre V1 (escucharás un soplo: es normal, equilibrio de presiones)",
                "  • V1 ABIERTA  V2 ABIERTA  V3 CERRADA",
                f"  • τ del sistema total ≈ {tau:.0f} s  ({tau/60:.1f} min)",
                "","──── PASO 4 [t≈1–5 min]  EVACUACIÓN VISCOSA ────",
                "  • Sin cambios de válvulas — deja trabajar la bomba",
                f"  • t=1 min → M1 ≈ 50–100 Torr",
                f"  • t=2 min → M1 ≈ 5–20 Torr",
                f"  • t=5 min → M1 ≈ 0.1–1 Torr",
                "  ⚠ Si se estanca >1 Torr por más de 2 min → busca fuga",
                "","──── PASO 5 [t≈5–{:.0f} min]  VACÍO PROFUNDO ────".format(t_ev/60),
                f"  • Meta: ≤ {pv:.2e} Torr",
                "  • Si tienes turbomolecular: enciéndela cuando p<0.1 Torr",
                "","──── PASO 6 [t≈{:.0f}s]  CONFIRMACIÓN Y PRUEBA ────".format(t_ev),
                "  • Cierra V2 (aísla la bomba)",
                "  • Observa M1 por 1 minuto sin tocar nada:",
                "    ✓ Sube < 0.01 Torr/min → SISTEMA ESTANCO → procede al llenado",
                "    ✗ Sube > 0.05 Torr/min → HAY FUGA → revisar conexiones y ferrules",
            ]
        else:
            lineas = [
                "═══ GUÍA DE LLENADO — PASO A PASO ═══",
                f"Meta: {po:.2f} ± 0.2 Torr  |  tasa: ~{tll:.4f} Torr/s  |  t_total ≈ {t_ll:.0f} s ({t_ll/60:.1f} min)",
                f"Masa necesaria: {mg:.4f} mg  |  Composición:"+
                "  ".join(f"{g}:{c['mezcla'][g]:.0f}%" for g in c['mezcla']),
                "","──── PASO 1 [t=0s]  ESTADO INICIAL ────",
                f"  • M1 debe leer ≤ {pv:.2e} Torr (evacuación completada)",
                "  • V1 ABIERTA  V2 CERRADA  V3 CERRADA",
                "","──── PASO 2 [t≈0–10s]  PURGA DE MANGUERA ────",
                "  • Presuriza la manguera de gas con el regulador (2–5 psi)",
                "  • Abre/cierra V3 rápido 3 veces para purgar el aire",
                "  • Deja V3 CERRADA al terminar",
                "","──── PASO 3 [t≈10–30s]  APERTURA CONTROLADA V3 ────",
                "  • Abre V3 muy suavemente (¼ vuelta máximo — válvula de aguja)",
                "  • M2 comenzará a subir desde 0 Torr",
                f"  • Si sube > {tll*2:.3f} Torr/s → cierra un poco V3",
                f"  • Si sube < {tll*0.5:.3f} Torr/s → abre un poco V3",
                "","──── PASO 4 [t≈30s–{:.0f}s]  LLENADO CONTINUO ────".format(t_ll*0.7),
                f"  • t=60s  → M2 ≈ {min(tll*60,po):.2f} Torr",
                f"  • t=120s → M2 ≈ {min(tll*120,po):.2f} Torr",
                f"  • t=180s → M2 ≈ {min(tll*180,po):.2f} Torr",
                "","──── PASO 5 [aproximación final] ────",
                f"  • Al llegar a {po*0.9:.2f} Torr (90%) → reduce apertura V3 a la mitad",
                f"  • Al llegar a {po*0.97:.2f} Torr (97%) → cierre casi total",
                "  ⚠ NO sobre-presiones — retirar gas es difícil",
                "","──── PASO 6  PRESIÓN OBJETIVO — CIERRE ────",
                f"  • Cuando M2 marque {po:.2f} ± 0.2 Torr: cierra V3 completamente",
                "  • V1 ABIERTA  V2 CERRADA  V3 CERRADA",
                "  • Espera 30 s y verifica estabilidad (estabilización térmica)",
                "","──── PASO 7  PRUEBA DE ENCENDIDO ────",
                "  • Enciende la fuente de alta tensión",
                "  • El plasma debe encender en ≤ 2 s",
                "  • Color esperado: azul-violeta (CO₂/He dominante)",
                "  ✓ Enciende suave → presión correcta",
                "  ✗ No enciende en 5s → verifica presión y voltaje de fuente",
                "  ✗ Arco ruidoso → presión demasiado alta → re-evacua y repite",
            ]
        dlg = ft.AlertDialog(
            modal=True,
            title=ft.Text(
                "📋 Guía de Evacuación" if tipo=="evacuacion" else "📋 Guía de Llenado",
                color=C["accent"],size=16,weight=ft.FontWeight.BOLD),
            content=ft.Container(width=600,height=500,
                content=ft.Column([
                    ft.Text(l,
                        color=C["accent"] if l.startswith("═") else
                              C["label"]  if l.startswith("────") else
                              C["text"],
                        size=12, font_family="Courier New",
                        weight=ft.FontWeight.BOLD if (l.startswith("PASO") or l.startswith("═") or l.startswith("────")) else ft.FontWeight.NORMAL)
                    for l in lineas], spacing=2, scroll=ft.ScrollMode.AUTO)),
            actions=[ft.TextButton("Cerrar",
                on_click=lambda ev: (setattr(dlg,"open",False), e.page.update()))],
            bgcolor=C["surface"])
        e.page.overlay.append(dlg); dlg.open=True; e.page.update()

    # ── Exportar ──────────────────────────────────────────────────────────────
    def exp_txt(self, e):
        if not self._cfg.get('p_obj'):
            self.snack("⚠️ Primero calcula.", C["yellow"]); e.page.update(); return
        try:
            c=self._cfg
            p=exportar("simulacion",{'vt':c['vt'],'vm':c['vm'],'vtot':c['vtot'],
                'pobj':c['p_obj'],'vb':c.get('Vb',0),
                'masa':c.get('mg',0),'mezcla':c['mezcla']})
            self.snack(f"📄 {p}", C["accent"])
        except Exception as ex: self.snack(f"❌ {ex}", C["red"])
        e.page.update()
    def exp_png(self, e):
        if not self._img_sim.visible:
            self.snack("⚠️ Primero genera las gráficas.", C["yellow"]); e.page.update(); return
        try: p=exp_img(self._img_sim.src,"simulacion"); self.snack(f"🖼️ {p}", C["accent"])
        except Exception as ex: self.snack(f"❌ {ex}", C["red"])
        e.page.update()

    # ── Navegación de sub-pestañas ─────────────────────────────────────────────
    def _switch_tab(self, idx, page):
        self._tab_idx[0] = idx
        for i,b in enumerate(self._tab_btns):
            b.bgcolor = C["b_blue"] if i==idx else C["surface2"]
            b.border  = ft.Border.all(1, C["border_f"] if i==idx else C["border"])
        self._tab_content.controls = [self._build_tab(idx)]
        if page: page.update()

    def _tab_btn(self, label, idx):
        b = ft.Container(
            content=ft.Text(label, color=C["text"], size=12,
                            weight=ft.FontWeight.W_600, text_align=ft.TextAlign.CENTER),
            bgcolor=C["surface2"], border_radius=8, padding=ft.Padding.symmetric(horizontal=14,vertical=8),
            border=ft.Border.all(1, C["border"]),
            on_click=lambda e, i=idx: self._switch_tab(i, e.page), ink=True)
        self._tab_btns.append(b); return b

    def _build_tab(self, idx):
        if idx == self.TAB_CONFIG:  return self._tab_config()
        if idx == self.TAB_RES:     return self._tab_resultados()
        if idx == self.TAB_SIM:     return self._tab_simulacion()
        if idx == self.TAB_GUIAS:   return self._tab_guias()
        return ft.Column([])

    # ── Tab 0: Configuración ──────────────────────────────────────────────────
    def _tab_config(self):
        # Buscar la imagen del manifold en varias ubicaciones posibles
        posibles_rutas = [
            "manifold.png",                          # Misma carpeta que el script
            os.path.join(os.path.dirname(__file__), "manifold.png"),  # Carpeta del script
            os.path.join(os.getcwd(), "manifold.png"),                # Carpeta actual
        ]
        
        ruta_imagen = None
        for ruta in posibles_rutas:
            if os.path.exists(ruta):
                ruta_imagen = ruta
                break
        
        # Si no encuentra la imagen, usar texto explicativo
        if ruta_imagen:
            try:
                img_manifold = ft.Image(
                    src=ruta_imagen,
                    width=600, fit=ft.ImageFit.CONTAIN,
                    error_content=ft.Text("(imagen manifold.png no encontrada)",
                                           color=C["muted"],size=11,italic=True))
            except AttributeError:
                # Versión antigua de Flet
                img_manifold = ft.Image(
                    src=ruta_imagen,
                    width=600,
                    error_content=ft.Text("(imagen manifold.png no encontrada)",
                                           color=C["muted"],size=11,italic=True))
        else:
            # Si no hay imagen, mostrar diagrama ASCII
            img_manifold = ft.Container(
                content=ft.Text(
                    """
    Bomba de                                            Mezcla
     Vacío                                             de Gases
       │                                                   │
       │         [M_V2]         [M_V3]                     │
       │           │              │                        │
       └───[V2]────┴──────────────┴────[V3]────────────────┘
                   │              │
                   │      [V1]    │
                   │       │      │
                   └───────┴──────┘
                           │
                      Tubo Láser
                    
Coloca la imagen 'manifold.png' en la misma carpeta que este script
                    """,
                    color=C["muted"], size=11, font_family="Courier New"
                ),
                bgcolor=C["surface2"],
                padding=10,
                border_radius=8
            )
        return ft.Column([
            ft.Container(
                content=ft.Column([
                    ft.Text("¿Qué es el re-gassing?", color=C["accent"],size=14,
                            weight=ft.FontWeight.BOLD),
                    ft.Text(
                        "El tubo láser RECI es un tubo de vidrio sellado con una mezcla "
                        "de gases (CO₂/N₂/H₂/He). Con el tiempo el gas se contamina o "
                        "se agota y el láser pierde potencia. El re-gassing consiste en: "
                        "1) evacuar el tubo con una bomba de vacío hasta ~10⁻³ Torr, "
                        "2) inyectar la mezcla fresca hasta la presión de operación "
                        "(~14.4 Torr para RECI W4 en CDMX). El manifold actúa como "
                        "estación de control con tres válvulas.",
                        color=C["muted"],size=13),
                    ft.Container(height=8),
                    ft.Row([img_manifold], alignment=ft.MainAxisAlignment.CENTER),
                    ft.Text("Diagrama del manifold: Bomba de vacío (izq) → V2 → colector central → V1 → Tubo de vidrio (abajo) / V3 → Mezcla de gases (der)",
                            color=C["muted"],size=11,italic=True,
                            text_align=ft.TextAlign.CENTER),
                ], spacing=6, horizontal_alignment=ft.CrossAxisAlignment.CENTER),
                bgcolor=C["surface2"], border_radius=12, padding=16,
                border=ft.Border.all(1, C["border"])),
            ft.Container(height=4),
            # Explicación detallada de Geometría
            ft.Container(
                content=ft.Column([
                    ft.Text("📏 Sobre la Geometría del Tubo", size=13, weight=ft.FontWeight.W_600, color=C["accent"]),
                    ft.Text(
                        "El volumen interno del tubo láser se puede calcular geométricamente "
                        "como un cilindro (π·r²·L), pero el volumen REAL suele ser menor debido a:\n"
                        "  • Electrodos internos que ocupan espacio\n"
                        "  • Capilares de menor diámetro en algunas secciones\n"
                        "  • Variaciones en el diámetro interno del vidrio\n\n"
                        "Para mayor precisión, activa 'Usar volumen manual' e ingresa el volumen "
                        "medido experimentalmente (ej: 5500 cm³ para RECI W4).",
                        size=11, color=C["muted"], italic=True
                    ),
                ]),
                bgcolor=C["surface2"], padding=12, border_radius=8,
                border=ft.border.all(1, C["border"]), margin=ft.margin.only(bottom=10)
            ),
            card([
                ft.Text("🔧 Geometría del Tubo RECI",color=C["label"],size=13,
                        weight=ft.FontWeight.W_600),
                ft.Row([self.chk_volm]),
                ft.Row([self.tf_long,self.tf_diam,self.tf_volm],spacing=10,wrap=True),
                nota("💡 Activa 'Usar volumen manual' para el volumen medido directamente."),
            ]),
            # Explicación detallada del Manifold
            ft.Container(
                content=ft.Column([
                    ft.Text("🔧 Sobre el Volumen Muerto del Manifold", size=13, weight=ft.FontWeight.W_600, color=C["accent"]),
                    ft.Text(
                        "El 'volumen muerto' es el volumen adicional de las mangueras y el cuerpo "
                        "del manifold que NO está dentro del tubo, pero SÍ se llena con gas.\n\n"
                        "Componentes del volumen muerto:\n"
                        "  • Cuerpo del manifold de latón (~50 cm³)\n"
                        "  • Manguera de gas (tanque → manifold)\n"
                        "  • Manguera de vacío (bomba → manifold)\n"
                        "  • Manguera al láser (manifold → tubo)\n\n"
                        "Este volumen AUMENTA la cantidad total de gas necesaria. Mangueras más "
                        "largas o de mayor diámetro requieren más gas para el mismo resultado.",
                        size=11, color=C["muted"], italic=True
                    ),
                ]),
                bgcolor=C["surface2"], padding=12, border_radius=8,
                border=ft.border.all(1, C["border"]), margin=ft.margin.only(bottom=10)
            ),
            card([
                ft.Text("⚙️ Volumen Muerto del Manifold",color=C["label"],size=13,
                        weight=ft.FontWeight.W_600),
                nota("Incluye: cuerpo de latón + mangueras de gas, bomba y láser."),
                ft.Row([self.tf_vcuerpo],spacing=10),
                ft.Row([self.tf_lgg,self.tf_dgg],spacing=10),
                ft.Row([self.tf_lgb,self.tf_dgb],spacing=10),
                ft.Row([self.tf_lgl,self.tf_dgl],spacing=10),
            ]),
            card([
                ft.Text("⚡ Parámetros de Operación",color=C["label"],size=13,
                        weight=ft.FontWeight.W_600),
                ft.Row([self.tf_pobj,self.tf_patm,self.tf_pvac],spacing=10,wrap=True),
                ft.Row([self.tf_vign,self.tf_vop,self.tf_imax],spacing=10,wrap=True),
                ft.Row([self.tf_delec,self.tf_tev,self.tf_tll],spacing=10,wrap=True),
                ft.Row([self.tf_gamma],spacing=10,wrap=True),
                nota("γ (gamma): Coeficiente de emisión secundaria de electrones. Típico: 0.01"),
            ]),
            card([
                ft.Text("⚗️ Composición de la Mezcla (%) — SparkLaser RECI",
                        color=C["label"],size=13,weight=ft.FontWeight.W_600),
                ft.Row([self.tf_co2,self.tf_n2,self.tf_h2,self.tf_he],spacing=10,wrap=True),
                nota("Suma debe ser ~100%. También puedes enviar desde ⚗️ Mezcla."),
            ]),
            ft.Row([
                btn("🔢 Calcular sistema", self.calcular, width=200),
            ], spacing=10),
            nota(f"📁 Exportaciones: {EXPORT_DIR}",C["muted"]),
        ], spacing=12, scroll=ft.ScrollMode.AUTO)

    # ── Tab 1: Resultados ─────────────────────────────────────────────────────
    def _tab_resultados(self):
        return ft.Column([
            ft.Row([
                btn("🔄 Recalcular",        self.calcular,  color=C["b_blue"]),
                btn("📊 Graficar resultados",self.graficar,  color=C["b_gray"]),
                btn("📄 Exportar .txt",      self.exp_txt,   color=C["b_green"]),
                btn("🖼️ Exportar .png",      self.exp_png,   color=C["b_green"]),
            ],spacing=10,wrap=True),
            sep(),
            card([self._res_col]),
            ft.Container(height=8),
            ft.Row([self._img_sim],alignment=ft.MainAxisAlignment.CENTER),
        ], spacing=12, scroll=ft.ScrollMode.AUTO)

    # ── Tab 2: Simulación ─────────────────────────────────────────────────────
    def _tab_simulacion(self):
        bv1 = self._make_vbtn("V1"); bv2 = self._make_vbtn("V2"); bv3 = self._make_vbtn("V3")
        bb  = self._make_eqbtn("bomba",  "💨 Bomba  OFF")
        bf  = self._make_eqbtn("fuente", "⚡ Fuente OFF")
        bt  = self._make_eqbtn("tubo",   "🔬 Tubo  DESMONTADO")
        bta = self._make_eqbtn("tanque", "🛢 Tanque  DESCONECTADO")

        return ft.Column([
            # Panel SCADA indicadores
            ft.Text("Panel SCADA — Estado en Tiempo Real",
                    color=C["label"],size=13,weight=ft.FontWeight.W_600),
            ft.Container(
                content=ft.Column([
                    ft.Row([
                        ind_card("PRESIÓN",    self._lbl_p,    C["accent"],"Torr actual",155),
                        ind_card("FASE",       self._lbl_fase, C["yellow"],"",155),
                        ind_card("TIEMPO",     self._lbl_t,    C["muted"], "simulado",130),
                        ind_card("MASA GAS",   self._lbl_masa, C["green"], "calculada",145),
                        ind_card("Vb PASCHEN", self._lbl_vb,   C["accent"],"ruptura",140),
                        ind_card("VOL. TOTAL", self._lbl_vtot, C["muted"], "cm³",140),
                    ],spacing=8,wrap=True),
                    ft.Container(height=6),
                    ft.Row([ft.Text("Progreso:",color=C["muted"],size=12),
                            self._prog_bar, self._prog_lbl],
                           spacing=8, vertical_alignment=ft.CrossAxisAlignment.CENTER),
                ],spacing=6),
                bgcolor=C["surface2"],border_radius=12,padding=14,
                border=ft.Border.all(1,C["border"])),

            # Control de la simulación
            ft.Text("Control de Simulación",color=C["label"],size=13,weight=ft.FontWeight.W_600),
            ft.Container(
                content=ft.Column([
                    ft.Row([
                        btn("▶️ INICIAR Simulación",self.iniciar_sim,color=C["b_green"],width=220),
                        btn("⏸ Pausar / Reanudar",self.toggle_pausa,color=C["b_gray"],width=210),
                    ],spacing=10,wrap=True),
                    ft.Container(height=4),
                    ft.Text("Velocidad de simulación:",color=C["muted"],size=12),
                    ft.Row([
                        btn("×1  Normal",    lambda e: self.set_speed(1.0, e), color=C["b_blue"], width=130),
                        btn("×3  Rápido",    lambda e: self.set_speed(3.0, e), color=C["b_violet"],width=130),
                        btn("×5  Muy rápido",lambda e: self.set_speed(5.0, e), color=C["b_orange"],width=150),
                        btn("×10 Turbo",     lambda e: self.set_speed(10.0,e), color=C["b_red"],   width=130),
                    ],spacing=10,wrap=True),
                    nota("La simulación se abre en una ventana matplotlib separada con animación en tiempo real."),
                    nota("Puedes cambiar la velocidad antes o durante la simulación.", C["muted"]),
                ],spacing=8),
                bgcolor=C["surface2"],border_radius=12,padding=14,
                border=ft.Border.all(1,C["border"])),

            # Control de válvulas
            ft.Text("Control Manual de Válvulas",color=C["label"],size=13,weight=ft.FontWeight.W_600),
            ft.Container(
                content=ft.Column([
                    nota("Alterna cada válvula — el sistema valida la secuencia de operación.",C["muted"]),
                    ft.Row([
                        ft.Column([
                            ft.Text("V1 — Maestra\n(manifold ↔ tubo)",
                                    color=C["muted"],size=10,text_align=ft.TextAlign.CENTER),
                            bv1,
                        ],spacing=4,horizontal_alignment=ft.CrossAxisAlignment.CENTER),
                        ft.Column([
                            ft.Text("V2 — Vacío\n(manifold ↔ bomba)",
                                    color=C["muted"],size=10,text_align=ft.TextAlign.CENTER),
                            bv2,
                        ],spacing=4,horizontal_alignment=ft.CrossAxisAlignment.CENTER),
                        ft.Column([
                            ft.Text("V3 — Inyección\n(manifold ↔ gas)",
                                    color=C["muted"],size=10,text_align=ft.TextAlign.CENTER),
                            bv3,
                        ],spacing=4,horizontal_alignment=ft.CrossAxisAlignment.CENTER),
                    ],spacing=20,wrap=True),
                    sep(),
                    ft.Text("Equipos del sistema:",color=C["muted"],size=12),
                    ft.Row([
                        ft.Column([ft.Text("Bomba de vacío",color=C["muted"],size=10,
                                           text_align=ft.TextAlign.CENTER),bb],
                                  spacing=4,horizontal_alignment=ft.CrossAxisAlignment.CENTER),
                        ft.Column([ft.Text("Fuente HV",color=C["muted"],size=10,
                                           text_align=ft.TextAlign.CENTER),bf],
                                  spacing=4,horizontal_alignment=ft.CrossAxisAlignment.CENTER),
                        ft.Column([ft.Text("Tubo RECI",color=C["muted"],size=10,
                                           text_align=ft.TextAlign.CENTER),bt],
                                  spacing=4,horizontal_alignment=ft.CrossAxisAlignment.CENTER),
                        ft.Column([ft.Text("Tanque mezcla",color=C["muted"],size=10,
                                           text_align=ft.TextAlign.CENTER),bta],
                                  spacing=4,horizontal_alignment=ft.CrossAxisAlignment.CENTER),
                    ],spacing=16,wrap=True),
                    nota("💡 Los botones de válvulas simulan el estado físico del manifold. "
                         "El sistema avisa si la combinación es insegura (ej. vacío + gas juntos).",
                         C["muted"]),
                ],spacing=8),
                bgcolor=C["surface2"],border_radius=12,padding=14,
                border=ft.Border.all(1,C["border"])),
        ], spacing=12, scroll=ft.ScrollMode.AUTO)

    # ── Tab 3: Guías ──────────────────────────────────────────────────────────
    def _tab_guias(self):
        return ft.Column([
            ft.Text("Guías Paso a Paso del Proceso",
                    color=C["label"],size=13,weight=ft.FontWeight.W_600),
            ft.Row([
                btn("📋 Guía de Evacuación",lambda e: self._guia_dlg("evacuacion",e),
                    color=C["b_orange"],width=230),
                btn("📋 Guía de Llenado",   lambda e: self._guia_dlg("llenado",e),
                    color=C["b_orange"],width=230),
            ],spacing=12,wrap=True),
            ft.Container(height=4),
            card([
                ft.Text("📐 Resumen del proceso de re-gassing",
                        color=C["accent"],size=14,weight=ft.FontWeight.BOLD),
                ft.Container(height=4),
                ft.Text("1. PREPARACIÓN",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                ft.Text("   Monta el tubo en el manifold (V1). Conecta la bomba (V2). Conecta el cilindro de gas (V3).\n"
                        "   Todas las válvulas CERRADAS al inicio.",color=C["muted"],size=12),
                ft.Text("2. EVACUACIÓN",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                ft.Text("   Enciende la bomba → abre V2 (purga manifold ~30 s) → abre V1 (evacúa el tubo).\n"
                        "   Meta: ≤ 10⁻³ Torr. Verifica estanqueidad cerrando V2 por 1 min.",
                        color=C["muted"],size=12),
                ft.Text("3. LLENADO",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                ft.Text("   Cierra V2 → purga manguera de gas 3 veces → abre V3 muy suavemente (válvula de aguja).\n"
                        "   Controla la tasa con M2 (Baratron). Objetivo ±0.2 Torr.",
                        color=C["muted"],size=12),
                ft.Text("4. VERIFICACIÓN",color=C["label"],size=13,weight=ft.FontWeight.W_600),
                ft.Text("   Cierra V3. Espera 30 s (estabilización térmica). Enciende fuente HV.\n"
                        "   El plasma debe encender en ≤2 s. Color azul-violeta = CO₂/He correcto.",
                        color=C["muted"],size=12),
                sep(),
                ft.Text("⚠️ Reglas de seguridad críticas",
                        color=C["red"],size=13,weight=ft.FontWeight.W_600),
                ft.Text("• NUNCA abras V2 y V3 al mismo tiempo (vacío + gas = contaminación).\n"
                        "• Verifica siempre la estanqueidad ANTES de introducir gas.\n"
                        "• No sobre-presiones: retirar gas del tubo es muy difícil.\n"
                        "• Usa siempre gafas al conectar/desconectar el cilindro de gas comprimido.",
                        color=C["muted"],size=12),
            ],pad=18),
        ], spacing=12, scroll=ft.ScrollMode.AUTO)

    # ── Build principal ───────────────────────────────────────────────────────
    def build(self):
        self._tab_btns.clear()
        tb = ft.Row([
            self._tab_btn("⚙️ Configuración", self.TAB_CONFIG),
            self._tab_btn("📊 Resultados",     self.TAB_RES),
            self._tab_btn("▶️ Simulación",     self.TAB_SIM),
            self._tab_btn("📋 Guías",          self.TAB_GUIAS),
        ], spacing=8, wrap=True)
        # marcar activa
        idx = self._tab_idx[0]
        if self._tab_btns:
            self._tab_btns[idx].bgcolor    = C["b_blue"]
            self._tab_btns[idx].border     = ft.Border.all(1, C["border_f"])
        self._tab_content.controls = [self._build_tab(idx)]
        return ft.Column([
            titulo("🔬  Simulador de Re-Gassing — Tubo Láser RECI"),
            tb,
            sep(),
            self._tab_content,
        ], spacing=10, expand=True)


# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN 6 — SIMULACIÓN EN VENTANA MATPLOTLIB
# ═══════════════════════════════════════════════════════════════════════════════

def _sim_ventana(cfg, dark, speed_ref, paused_ref, on_upd=None, on_end=None):
    """
    Ventana matplotlib con animación frame-a-frame del re-gassing.
    speed_ref: lista[float] compartida con la UI para cambiar velocidad en vivo.
    paused_ref: lista[bool] para pausar/reanudar.
    """
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.gridspec as gridspec

    p_atm, p_vac, p_obj = cfg['p_atm'], cfg['p_vac'], cfg['p_obj']
    tasa_ev, tasa_ll     = cfg['tasa_ev'], cfg['tasa_ll']
    tau  = p_atm / tasa_ev
    A, B, gamma  = cfg['Am'], cfg['Bm'], cfg.get('gamma', 0.01)
    d            = cfg['dist_electrodos']
    mezcla       = cfg.get('mezcla', {})
    v_ign        = cfg['v_ign']

    BG, AX, TX = ("#0d1117","#1c2333","#e6edf3") if dark else ("#f8fafc","#e2e8f0","#0f172a")
    CE, CL, CO, CV = "#58a6ff","#3fb950","#f78166","#ffa657"

    # ─── FASES según procedimiento real de re-gassing SparkLaser RECI ────────
    # Duración de evacuación calculada: tiempo para bajar de p_atm a p_vac
    _t_evac = max(120, int(tau * 8))   # ~8 tau para bajar ~99.97%
    _t_fill = max(60,  int(p_obj / max(tasa_ll, 1e-6)) + 20)
    _p_flush = 50.0  # Torr — presión de purga en cada ciclo de flushing

    # tipo: "idle"=presión congelada, "evac"=baja, "gas"=sube, "flush_up"=sube a 50T, "flush_dn"=baja de 50T
    FASES = [
        # ── Fase A: Preparación ──────────────────────────────────────────────
        ("A: VERIFICACIÓN",    6,   "idle",     False, False, False,
         "Revisar vidrio, aceite de bomba y conexiones"),
        ("A: INSTALAR TUBO",   5,   "idle",     False, False, False,
         "Acoplar tubo RECI al conector Ultra-Torr"),
        ("A: ENCENDER BOMBA",  10,  "idle",     False, False, False,
         "Bomba mecánica ON — esperar estabilización (~10 s)"),
        # ── Fase B: Evacuación ───────────────────────────────────────────────
        ("B: PURGAR MANIFOLD", 30,  "evac",     False, True,  False,
         "V2 abierta — evacuando manifold y mangueras"),
        ("B: EVACUAR SISTEMA", _t_evac, "evac", True,  True,  False,
         "V1+V2 abiertas — evacuando tubo RECI completo hasta vacío"),
        # Prueba de estanqueidad: solo 60 s con V2 cerrada
        ("B: PRUEBA ESTANQ.",  60,  "idle",     True,  False, False,
         "V2 cerrada — vigilar 60 s: si presión sube hay fuga"),
        # ── Fase B: Flushing 3 ciclos ────────────────────────────────────────
        ("B: FLUSH CICLO 1 ↑", 18, "flush_up",  True,  False, True,
         "Ciclo 1/3 — inyectar Premix hasta ~50 Torr (abrir V3)"),
        ("B: FLUSH CICLO 1 ↓", 25, "flush_dn",  True,  True,  False,
         "Ciclo 1/3 — evacuar gas de purga (cerrar V3, abrir V2)"),
        ("B: FLUSH CICLO 2 ↑", 18, "flush_up",  True,  False, True,
         "Ciclo 2/3 — inyectar Premix hasta ~50 Torr"),
        ("B: FLUSH CICLO 2 ↓", 25, "flush_dn",  True,  True,  False,
         "Ciclo 2/3 — evacuar gas de purga"),
        ("B: FLUSH CICLO 3 ↑", 18, "flush_up",  True,  False, True,
         "Ciclo 3/3 — último ciclo de purga del cátodo"),
        ("B: FLUSH CICLO 3 ↓", 25, "flush_dn",  True,  True,  False,
         "Ciclo 3/3 — re-evacuar hasta vacío profundo"),
        # ── Fase C: Llenado final ────────────────────────────────────────────
        ("C: LLENADO FINAL",   _t_fill, "gas",  True,  False, True,
         "V3 (aguja) abierta lentamente — inyectando mezcla final"),
        # ── Fase D: Diagnóstico y sellado ────────────────────────────────────
        ("D: ESTABILIZACIÓN",  30,  "done",     True,  False, False,
         "Esperar: verificar que presión NO suba (micro-fugas)"),
        ("D: DIAGNÓSTICO",     12,  "done",     True,  False, False,
         "Encender HV — color ROSA/PÚRPURA = mezcla correcta ✓"),
        ("D: SELLADO",         5,   "done",     False, False, False,
         "Cerrar V1 — retirar tubo RECI — proceso completado ✓"),
    ]
    t_total = sum(f[1] for f in FASES)

    # Presión acumulada al final de cada fase (para p_en_t correcto)
    def _p_end_of_phase(tipo, dur, p_in, tau_loc):
        if tipo == "evac":      return max(p_in * math.exp(-dur / tau_loc), p_vac)
        if tipo == "gas":       return min(p_in + tasa_ll * dur, p_obj)
        if tipo == "flush_up":  return _p_flush
        if tipo == "flush_dn":  return max(_p_flush * math.exp(-dur / tau_loc), p_vac)
        return p_in  # idle / done

    def p_en_t(t_glob):
        acum = 0; p = p_atm
        for _, dur, tipo, v1, v2, v3, _ in FASES:
            dt = t_glob - acum
            if dt <= 0: return p
            if dt <= dur:
                # Presión dentro de esta fase
                if tipo == "evac":
                    return max(p * math.exp(-dt / tau), p_vac)
                if tipo == "gas":
                    return min(p + tasa_ll * dt, p_obj)
                if tipo == "flush_up":
                    # Sube de p_vac hacia _p_flush en dur segundos
                    tasa_f = (_p_flush - p) / dur
                    return min(p + tasa_f * dt, _p_flush)
                if tipo == "flush_dn":
                    return max(p * math.exp(-dt / tau), p_vac)
                return p  # idle / done
            # Avanzar p al final de esta fase
            p = _p_end_of_phase(tipo, dur, p, tau)
            acum += dur
        return p_obj

    # ── Layout de la ventana ──────────────────────────────────────────────────
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=(18, 10), dpi=95)
    fig.patch.set_facecolor(BG)
    try: fig.canvas.manager.set_window_title("SCADA — Re-Gassing RECI v12")
    except: pass

    # Grilla: 2 filas × 3 cols
    # Fila 0: presión (cols 0-1) + manifold (col 2, más ancho)
    # Fila 1: Paschen + estado + pastel
    gs_ = gridspec.GridSpec(2, 3, figure=fig,
                             height_ratios=[1.55, 1.0],
                             width_ratios=[1.3, 1.0, 1.1],
                             hspace=0.50, wspace=0.38,
                             left=0.06, right=0.97, top=0.93, bottom=0.07)

    # ── Panel 1: Presión vs tiempo ────────────────────────────────────────────
    ax_p = fig.add_subplot(gs_[0, :2])
    ax_p.set_facecolor(AX); ax_p.set_yscale('log')
    ax_p.set_xlim(0, t_total)
    ax_p.set_ylim(p_vac * 0.3, p_atm * 2.5)
    ax_p.set_xlabel("Tiempo (s)", color=TX, fontsize=9)
    ax_p.set_ylabel("Presión (Torr) [log]", color=TX, fontsize=9)
    ax_p.set_title("PRESIÓN vs TIEMPO — Procedimiento completo de Re-Gassing (tiempo real)",
                    color=TX, fontsize=10, fontweight='bold', pad=6)
    ax_p.axhline(p_obj, color=CO,  ls='--', lw=1.6, label=f"P obj = {p_obj:.2f} Torr", zorder=3)
    ax_p.axhline(p_vac, color=CV,  ls=':',  lw=1.2, label=f"Vacío = {p_vac:.2e} Torr",  zorder=3)
    ax_p.axhline(50.0,  color=CL,  ls=':',  lw=0.8, alpha=0.5, label="P flush = 50 Torr", zorder=2)

    # Bandas de color por tipo de fase
    _fc_map = {'idle':'#1e293b','evac':'#1e3a8a','flush_up':'#4c1d95',
               'flush_dn':'#5b21b6','gas':'#14532d','done':'#78350f'}
    _t_ac = 0
    for _fn, _fd, _ft, *_ in FASES:
        _fc = _fc_map.get(_ft, '#1e293b')
        ax_p.axvspan(_t_ac, _t_ac + _fd, alpha=0.13, color=_fc, zorder=1)
        ax_p.axvline(_t_ac, color='#475569', lw=0.5, alpha=0.35, zorder=2)
        # Etiqueta solo si la fase es suficientemente ancha
        if _fd >= t_total * 0.04:
            _mid = _t_ac + _fd / 2
            _lbl = _fn.split(':')[-1].strip()[:10]
            ax_p.text(_mid, p_atm * 1.6, _lbl,
                      color=TX, fontsize=5.5, ha='center', va='bottom',
                      alpha=0.75, rotation=55, clip_on=True)
        _t_ac += _fd

    ax_p.tick_params(colors=TX, labelsize=8)
    for sp in ax_p.spines.values(): sp.set_edgecolor('#30363d' if dark else '#94a3b8')
    ax_p.grid(True, ls='--', alpha=0.18, color='#8b949e')
    ax_p.legend(fontsize=8, facecolor=AX, labelcolor=TX, loc='upper right', framealpha=0.9)
    line_p, = ax_p.plot([], [], color=CE, lw=2.2, zorder=5)
    dot_p,  = ax_p.plot([], [], 'o', color=CO, ms=8, zorder=8)
    txt_fase = ax_p.text(0.02, 0.97, '', transform=ax_p.transAxes,
                         color='#fbbf24', fontsize=9, va='top', fontweight='bold',
                         bbox=dict(facecolor=AX, edgecolor='#475569', alpha=0.92, pad=4))

    # ── Panel 2: Diagrama manifold (animado) ─────────────────────────────────
    ax_m = fig.add_subplot(gs_[0, 2])
    ax_m.set_facecolor(AX); ax_m.axis('off')

    def draw_manifold(v1, v2, v3, pr, fase_desc=""):
        ax_m.cla(); ax_m.set_facecolor(AX)
        ax_m.set_xlim(0, 10); ax_m.set_ylim(0, 12)
        ax_m.axis('off')
        ax_m.set_title("Manifold — Estado de Válvulas", color=TX,
                        fontsize=10, fontweight='bold', pad=6)

        # Color del tubo principal según modo activo
        pc = CE if (v2 and not v3) else (CL if v3 else '#334155')

        # ── Tubo principal horizontal ──────────────────────────────────────
        ax_m.plot([2.0, 8.0], [6.0, 6.0], color=pc, lw=6,
                  solid_capstyle='round', zorder=2)

        # ── Bomba (izquierda) ──────────────────────────────────────────────
        fc_b = CE if v2 else '#2d3748'
        ec_b = CE if v2 else '#475569'
        ax_m.add_patch(mpatches.FancyBboxPatch(
            (0.1, 5.0), 1.6, 2.0,
            boxstyle="round,pad=0.12", facecolor=fc_b, edgecolor=ec_b, lw=2.0))
        ax_m.text(0.9, 6.0, '⚙️\nBOMBA', color='#ffffff', fontsize=7,
                  ha='center', va='center', fontweight='bold')
        ax_m.plot([1.7, 2.0], [6.0, 6.0],
                  color=CE if v2 else '#334155', lw=4)
        ax_m.text(0.9, 4.7, 'V2', color=CE if v2 else '#e24b4a',
                  fontsize=7, ha='center', fontweight='bold')

        # ── Tanque gas (derecha) ───────────────────────────────────────────
        fc_g = CL if v3 else '#2d3748'
        ec_g = CL if v3 else '#475569'
        ax_m.add_patch(mpatches.FancyBboxPatch(
            (8.3, 4.2), 1.5, 3.6,
            boxstyle="round,pad=0.2", facecolor=fc_g, edgecolor=ec_g, lw=2.0))
        ax_m.text(9.05, 6.0, '🛢️\nGAS', color='#ffffff', fontsize=7,
                  ha='center', va='center', fontweight='bold')
        ax_m.plot([8.0, 8.3], [6.0, 6.0],
                  color=CL if v3 else '#334155', lw=4)
        ax_m.text(9.05, 4.0, 'V3', color=CL if v3 else '#e24b4a',
                  fontsize=7, ha='center', fontweight='bold')

        # ── Cuerpo del manifold ────────────────────────────────────────────
        ax_m.add_patch(mpatches.FancyBboxPatch(
            (3.6, 5.3), 2.8, 1.4,
            boxstyle="round,pad=0.08", facecolor='#1e3a5f',
            edgecolor='#38bdf8', lw=2.0))
        ax_m.text(5.0, 6.0, 'MANIFOLD', color='#38bdf8', fontsize=8,
                  ha='center', va='center', fontweight='bold')

        # ── Válvula V1 (abajo del manifold hacia el tubo) ──────────────────
        v1_col = '#4ade80' if v1 else '#e24b4a'
        ax_m.add_patch(mpatches.RegularPolygon(
            (5.0, 4.6), 4, radius=0.42, orientation=math.pi/4,
            facecolor=v1_col, edgecolor='white', lw=2.0, zorder=5))
        ax_m.text(5.0, 4.0, 'V1', color=v1_col, fontsize=7.5,
                  ha='center', va='top', fontweight='bold')
        ax_m.plot([5.0, 5.0], [5.3, 4.6 + 0.42],
                  color=v1_col, lw=4, zorder=3)
        ax_m.plot([5.0, 5.0], [4.6 - 0.42, 3.5],
                  color=v1_col if v1 else '#334155', lw=4, zorder=3)

        # ── Tubo RECI (abajo) ──────────────────────────────────────────────
        p_ratio = 0.0
        if p_obj > p_vac:
            p_ratio = min(1.0, max(0.0, (pr - p_vac) / (p_obj - p_vac)))
        fc_t = '#0d2040' if v1 else '#111827'
        ec_t = '#38bdf8' if v1 else '#475569'
        ax_m.add_patch(mpatches.FancyBboxPatch(
            (3.0, 0.6), 4.0, 2.7,
            boxstyle="round,pad=0.12", facecolor=fc_t, edgecolor=ec_t, lw=2.2))
        if p_ratio > 0.01:
            # barra de "llenado"
            ax_m.add_patch(mpatches.FancyBboxPatch(
                (3.1, 0.7), 3.8 * p_ratio, 2.5,
                boxstyle="square,pad=0",
                facecolor=CE + '55', edgecolor='none'))
        ax_m.text(5.0, 2.0, f'TUBO RECI\n{pr:.4g} Torr',
                  color='#38bdf8' if v1 else '#64748b',
                  fontsize=7.5, ha='center', va='center', fontweight='bold')

        # ── Manómetros ─────────────────────────────────────────────────────
        for mx, my, ml, mc in [(3.5, 8.0, 'M1\nPirani', '#38bdf8'),
                                (6.5, 8.0, 'M2\nBaratron', '#f0b429')]:
            ax_m.add_patch(mpatches.Circle(
                (mx, my), 0.62, facecolor='#0d1a2a', edgecolor=mc, lw=1.5))
            ax_m.text(mx, my, ml, color=mc, fontsize=6,
                      ha='center', va='center', fontweight='bold')
            ax_m.plot([mx, mx], [my - 0.62, 6.0 + 0.7 if abs(mx - 5) > 1 else 6.7],
                      color=mc, lw=1.2, ls='--', alpha=0.6)

        # ── Válvulas en el tubo horizontal (V2 y V3 como diamantes) ────────
        for vx, vn, vst, vc in [(2.0,'V2',v2,'#38bdf8'), (8.0,'V3',v3,'#f97316')]:
            fc_v = vc if vst else '#e24b4a'
            ax_m.add_patch(mpatches.RegularPolygon(
                (vx, 6.0), 4, radius=0.40, orientation=math.pi/4,
                facecolor=fc_v, edgecolor='white', lw=1.8, zorder=6))
            lbl = f'{vn}\n{"A" if vst else "C"}'
            ax_m.text(vx, 6.95, lbl, color=fc_v, fontsize=6,
                      ha='center', va='bottom', fontweight='bold')

        # ── Descripción de la fase actual ──────────────────────────────────
        ax_m.text(5.0, 11.4, fase_desc[:42],
                  color='#fbbf24', fontsize=6.5, ha='center', va='top',
                  style='italic',
                  bbox=dict(facecolor=AX, edgecolor='#475569', alpha=0.85, pad=3))

    draw_manifold(False, False, False, p_atm, "Sistema en estado inicial")

    # Panel 3: Curva Paschen con punto dinámico
    ax_pch = fig.add_subplot(gs_[1,0])
    ax_pch.set_facecolor(AX)
    pd_a = np.logspace(-1,3,800)
    Vb_a = np.array([paschen_Vb(pd,A,B,gamma) for pd in pd_a])
    ok   = np.isfinite(Vb_a)&(Vb_a>0)
    ax_pch.plot(pd_a[ok]/d, Vb_a[ok]/1000, color='#d2a8ff', lw=2)
    Vbop = paschen_Vb(p_obj*d,A,B,gamma)
    if math.isfinite(Vbop):
        ax_pch.plot(p_obj, Vbop/1000,'o',color=CO,ms=8,zorder=6,
                    label=f'Op: {Vbop/1000:.3f} kV')
    ax_pch.axhline(v_ign,color='#f0b429',ls='--',lw=1.3,label=f'V fuente={v_ign} kV')
    ax_pch.set_xscale('log')
    ax_pch.set_title('Curva de Paschen',color=TX,fontsize=9,fontweight='bold')
    ax_pch.set_xlabel('Presión (Torr)',color=TX,fontsize=8)
    ax_pch.set_ylabel('Vb (kV)',color=TX,fontsize=8)
    ax_pch.tick_params(colors=TX,labelsize=7)
    for sp in ax_pch.spines.values(): sp.set_edgecolor('#30363d')
    ax_pch.grid(True,ls='--',alpha=0.2,color='#8b949e')
    ax_pch.legend(fontsize=7.5,facecolor=AX,labelcolor=TX)
    dot_pch, = ax_pch.plot([],[],'*',color='#fbbf24',ms=13,zorder=9,label='p actual')

    # Panel 4: Indicadores numéricos
    ax_inf = fig.add_subplot(gs_[1,1])
    ax_inf.set_facecolor(AX); ax_inf.axis('off')
    ax_inf.set_title('Estado del Sistema',color=TX,fontsize=9,fontweight='bold')
    etiq = [('Presión:',0.90),('Fase:',0.77),('Tiempo:',0.64),
            ('V1 Maestra:',0.51),('V2 Vacío:',0.38),('V3 Gas:',0.25),('Progreso:',0.12)]
    tinf = {}
    for lbl,y in etiq:
        ax_inf.text(0.04,y,lbl,color=TX,fontsize=8.5,transform=ax_inf.transAxes)
        tinf[lbl] = ax_inf.text(0.52,y,'—',color='#38bdf8',fontsize=9,
                                  fontweight='bold',transform=ax_inf.transAxes)

    # Panel 5: Pie de mezcla (compacto)
    ax_pie = fig.add_subplot(gs_[1, 2])
    ax_pie.set_facecolor(AX)
    if mezcla:
        gl = list(mezcla.keys()); fl = [mezcla[g] for g in gl]
        wedges, texts, autotexts = ax_pie.pie(
            fl, labels=gl, autopct='%1.0f%%', startangle=90,
            colors=[CE, CL, '#d2a8ff', CV][:len(gl)],
            textprops={'color': TX, 'fontsize': 7},
            pctdistance=0.76, radius=0.72)
        for at in autotexts: at.set_fontsize(6)
    ax_pie.set_title('Composición\nMezcla', color=TX, fontsize=8,
                      fontweight='bold', pad=2)

    fig.suptitle('SCADA — SIMULACIÓN RE-GASSING  |  TUBO LÁSER RECI  v12',
                  color=TX, fontsize=11, fontweight='bold', y=0.995)

    plt.ion(); plt.show()

    # ── Loop de animación ─────────────────────────────────────────────────────
    BASE_DT  = 1.0 / 10.0   # 10 fps base
    t_global = 0.0
    t_hist   = []; p_hist = []
    fi = 0; ft_tick = 0.0

    while t_global <= t_total + 0.5:
        if not plt.fignum_exists(fig.number): break

        # Pausa
        while paused_ref[0]:
            if not plt.fignum_exists(fig.number): break
            try: fig.canvas.flush_events()
            except: break
            time.sleep(0.05)
        if not plt.fignum_exists(fig.number): break

        speed   = max(0.1, float(speed_ref[0]))
        dt      = BASE_DT * speed
        fi_real = min(fi, len(FASES)-1)
        fase_nom,dur,tipo,v1,v2,v3,desc = FASES[fi_real]
        pr  = p_en_t(t_global)
        prog= min(t_global/t_total, 1.0)

        t_hist.append(t_global); p_hist.append(max(pr, p_vac*0.9))

        # Actualizar gráficas
        line_p.set_data(t_hist, p_hist)
        dot_p.set_data([t_global],[pr])
        txt_fase.set_text(f"{fase_nom}: {desc}")
        Vb_now = paschen_Vb(pr*d, A, B, gamma)
        if math.isfinite(Vb_now) and Vb_now>0 and pr>0:
            dot_pch.set_data([pr],[Vb_now/1000])
        else:
            dot_pch.set_data([],[])
        draw_manifold(v1, v2, v3, pr, desc)
        # Indicadores
        tinf['Presión:'].set_text(f"{pr:.4f} Torr")
        tinf['Presión:'].set_color(CE if v2 else (CL if v3 else '#fbbf24'))
        tinf['Fase:'].set_text(fase_nom)
        _fase_col = {'evac':'#58a6ff','flush_up':'#d2a8ff','flush_dn':'#a78bfa',
                     'gas':'#3fb950','done':'#fbbf24','idle':'#94a3b8'}
        tinf['Fase:'].set_color(_fase_col.get(tipo, '#94a3b8'))
        tinf['Tiempo:'].set_text(f"{t_global:.1f} s  ({t_global/60:.1f} min)")
        for vn,vst in [('V1 Maestra:',v1),('V2 Vacío:',v2),('V3 Gas:',v3)]:
            tinf[vn].set_text('ABIERTA ●' if vst else 'CERRADA ○')
            tinf[vn].set_color('#4ade80' if vst else '#fb7185')
        tinf['Progreso:'].set_text(f"{prog*100:.1f} %")
        # Callback UI
        if on_upd: on_upd(fase_nom, pr, t_global, prog)
        try:
            fig.canvas.draw_idle(); fig.canvas.flush_events()
        except: break
        time.sleep(max(BASE_DT*0.3, 0.01))  # nunca duerme más de lo necesario
        # Avanzar tiempo y fase
        ft_tick += dt
        if ft_tick >= FASES[fi_real][1]:
            ft_tick = 0.0
            if fi < len(FASES)-1: fi += 1
        t_global += dt

    if on_end: on_end()
    plt.ioff()
    try: plt.show(block=False)
    except: pass


# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN 7 — AYUDA (con LaTeX y conceptos completos)
# ═══════════════════════════════════════════════════════════════════════════════

class AyudaView:
    def _bloque(self, emoji, tit, pasos, n=""):
        items = [ft.Row([
            ft.Container(
                content=ft.Text(str(i+1),color="#0ea5e9",size=12,weight=ft.FontWeight.BOLD),
                width=26,height=26,bgcolor=C["surface2"],border_radius=13,
                alignment=ft.Alignment(0,0)),
            ft.Text(p,color=C["text"],size=13,expand=True),
        ],spacing=10,vertical_alignment=ft.CrossAxisAlignment.START)
        for i,p in enumerate(pasos)]
        ch = [ft.Row([ft.Text(emoji,size=20),
                      ft.Text(tit,color=C["accent"],size=15,weight=ft.FontWeight.BOLD)],spacing=8),
              ft.Container(height=4),*items]
        if n: ch.append(ft.Container(content=nota(f"💡  {n}"),padding=ft.Padding.only(top=6)))
        return card(ch,pad=16)

    def _p(self, sym, name, unit, desc):
        return ft.Row([
            ft.Container(ft.Text(sym,color=C["accent"],size=14,
                                  weight=ft.FontWeight.BOLD,font_family="Courier New"),width=46),
            ft.Column([ft.Text(name,color=C["text"],size=13,weight=ft.FontWeight.W_600),
                       ft.Text(f"{unit}  —  {desc}",color=C["muted"],size=12)],
                      spacing=1,expand=True),
        ],spacing=8,vertical_alignment=ft.CrossAxisAlignment.START)

    def _sec(self, txt, color=None):
        return ft.Text(txt,color=color or C["accent"],size=15,weight=ft.FontWeight.BOLD)

    def _unidad_row(self, nombre, simb, equiv):
        return ft.Row([
            ft.Text(nombre,color=C["text"],size=13,width=140,weight=ft.FontWeight.W_600),
            ft.Text(simb,color=C["accent"],size=13,width=180,font_family="Courier New"),
            ft.Text(equiv,color=C["muted"],size=12,expand=True),
        ],spacing=8)

    def build(self):
        return ft.Column([
            titulo("📖  Guía de Uso — Paschen v12 + Simulador Re-Gassing"),
            ft.Container(
                content=ft.Column([
                    ft.Text("Calculadora Ley de Paschen + Simulador Re-Gassing RECI — v12",
                            size=17,weight=ft.FontWeight.BOLD,color=C["text"]),
                    ft.Container(height=4),
                    ft.Text(
                        "Esta herramienta combina la calculadora rigurosa de la Ley de Paschen "
                        "con un simulador físico completo del proceso de re-gassing para tubos "
                        "láser RECI CO₂. Puedes calcular presiones de operación óptimas, "
                        "simular el proceso evacuación→llenado, y seguir guías paso a paso.",
                        size=13,color=C["muted"]),
                    ft.Container(height=6),
                    ft.Text("Características principales:",size=13,color=C["label"],
                            weight=ft.FontWeight.W_600),
                    ft.Text(
                        "• Matemática corregida: A·pd > ln(1+1/γ) — condición física real.\n"
                        "• Mínimo analítico exacto: pd_min = e·ln(1+1/γ)/A.\n"
                        "• Bisección numérica robusta (error < 10⁻⁸ V).\n"
                        "• Mezcla de gases: A y B ponderados por fracción molar.\n"
                        "• Gas ideal PV=nRT para masa exacta de la mezcla.\n"
                        "• Simulación animada en ventana matplotlib (velocidad ×1/×3/×5/×10).\n"
                        "• Panel SCADA con indicadores en tiempo real.\n"
                        "• Control manual de válvulas con validación de secuencia.",
                        size=13,color=C["muted"]),
                ],spacing=4),
                bgcolor=C["surface2"],border_radius=12,padding=18,
                border=ft.Border.all(1,C["border"])),

            # ── Ecuaciones LaTeX ──────────────────────────────────────────────
            card([
                self._sec("⚛️  Ecuaciones Físicas (LaTeX)"),
                ft.Text("Ley de Paschen — Voltaje de ruptura:",
                        color=C["label"],size=13,weight=ft.FontWeight.W_600),
                latex_img(
                    r"V_b = \dfrac{B \cdot p \cdot d}{\ln(A \cdot p \cdot d)"
                    r" - \ln\!\left(\ln\!\left(1+\dfrac{1}{\gamma}\right)\right)}",
                    fontsize=16),
                ft.Text("Condición física del denominador:",color=C["muted"],size=12),
                latex_img(
                    r"A \cdot p \cdot d \;>\; \ln\!\left(1 + \dfrac{1}{\gamma}\right)",
                    fontsize=14),
                ft.Container(height=6),
                ft.Text("Mínimo de Paschen (analítico):",
                        color=C["label"],size=13,weight=ft.FontWeight.W_600),
                latex_img(
                    r"(p\cdot d)_{\min} = \dfrac{e\cdot\ln(1+1/\gamma)}{A}"
                    r"\qquad V_{b,\min} = B\cdot(p\cdot d)_{\min}",
                    fontsize=13),
                ft.Container(height=6),
                ft.Text("Constantes de Townsend para mezcla de gases:",
                        color=C["label"],size=13,weight=ft.FontWeight.W_600),
                latex_img(
                    r"A_{\mathrm{mezcla}} = \sum_i x_i A_i \qquad"
                    r"B_{\mathrm{mezcla}} = \sum_i x_i B_i",
                    fontsize=13),
                ft.Container(height=6),
                ft.Text("Masa de gas — Ley del Gas Ideal:",
                        color=C["label"],size=13,weight=ft.FontWeight.W_600),
                latex_img(r"PV = nRT \;\Rightarrow\; n = \dfrac{PV}{RT}", fontsize=14),
                ft.Text("R = 62.364 L·Torr/(mol·K)   |   T = 298 K (25 °C)",
                        color=C["muted"],size=12,font_family="Courier New"),
                ft.Container(height=6),
                ft.Text("Evacuación exponencial — bomba rotatoria:",
                        color=C["label"],size=13,weight=ft.FontWeight.W_600),
                latex_img(
                    r"p(t) = p_0 \cdot e^{-t/\tau} \qquad \tau = \dfrac{p_0}{\dot{p}_{ev}}",
                    fontsize=13),
                ft.Container(height=6),
                ft.Text("Llenado lineal — válvula de aguja:",
                        color=C["label"],size=13,weight=ft.FontWeight.W_600),
                latex_img(r"p(t) = p_{\mathrm{vac}} + \dot{p}_{ll} \cdot t", fontsize=13),
            ],pad=18),

            # ── Conceptos físicos ─────────────────────────────────────────────
            card([
                self._sec("🔬  Conceptos Físicos Clave"),
                ft.Container(height=4),
                ft.Text("Voltaje de ruptura (Vb)",color=C["label"],size=14,weight=ft.FontWeight.W_600),
                ft.Text(
                    "Es el voltaje mínimo al que un gas deja de ser aislante y se vuelve "
                    "conductor, permitiendo flujo de corriente entre los electrodos. "
                    "Por debajo de Vb el gas es dieléctrico. Por encima se produce "
                    "una descarga (arco, plasma o glow discharge). En el tubo láser CO₂ "
                    "esta descarga es el plasma que genera la inversión de población.",
                    color=C["muted"],size=13),
                ft.Container(height=6),
                ft.Text("El Mínimo de Paschen",color=C["label"],size=14,weight=ft.FontWeight.W_600),
                ft.Text(
                    "La curva Vb vs p·d tiene forma de U con mínimo bien definido. "
                    "A presiones muy bajas, los electrones no colisionan → más voltaje. "
                    "A presiones muy altas, muchas colisiones pero poca energía → más voltaje. "
                    "En el mínimo el balance es óptimo → voltaje mínimo de ruptura. "
                    "Operar cerca del mínimo maximiza la eficiencia del láser.",
                    color=C["muted"],size=13),
                ft.Container(height=6),
                ft.Text("Coeficiente γ de Townsend",color=C["label"],size=14,weight=ft.FontWeight.W_600),
                ft.Text(
                    "Describe cuántos electrones secundarios emite el cátodo por cada ion "
                    "positivo que lo impacta. Depende del material del cátodo y del gas. "
                    "Valores típicos: γ = 0.01–0.1 para metales comunes en gases nobles.",
                    color=C["muted"],size=13),
                ft.Container(height=6),
                ft.Text("Mezcla SparkLaser para RECI",color=C["label"],size=14,weight=ft.FontWeight.W_600),
                ft.Text(
                    "CO₂ 8% — Gas activo láser (transición vibracional 10.6 µm). "
                    "N₂ 15% — Molécula de transferencia; bombea CO₂ al nivel superior. "
                    "H₂ 5%  — Acelera la recuperación del nivel inferior (más velocidad de repetición). "
                    "He 72% — Gas de buffer; enfría el gas, ajusta Vb y mejora la conductividad térmica.",
                    color=C["muted"],size=13),
            ],pad=18),

            # ── Parámetros de la ecuación ─────────────────────────────────────
            card([
                self._sec("📐  Parámetros de la Ecuación"),
                ft.Container(height=4),
                self._p("Vb","Voltaje de ruptura","Voltios [V]",
                        "tensión mínima para ionizar el gas entre electrodos"),
                self._p("p", "Presión del gas","Torr",
                        "presión absoluta entre los electrodos"),
                self._p("d", "Distancia electrodos","cm",
                        "separación entre ánodo y cátodo"),
                self._p("A", "Constante de Townsend A","1/(Torr·cm)",
                        "coeficiente de ionización primaria por colisión; depende del gas"),
                self._p("B", "Constante de Townsend B","V/(Torr·cm)",
                        "relacionada con el potencial de ionización del gas"),
                self._p("γ", "Coeficiente de Townsend γ","adimensional",
                        "emisión secundaria del cátodo; típico 0.01–0.1"),
                self._p("τ", "Constante de tiempo evacuación","s",
                        "τ = p₀/ṗ_ev  (tiempo característico de la bomba)"),
            ],pad=18),

            # ── Unidades ──────────────────────────────────────────────────────
            card([
                self._sec("🔄  Unidades de Presión Disponibles"),
                ft.Text("Escribe cualquiera de estas (mayúsculas, minúsculas o mixtas):",
                        color=C["muted"],size=13),
                ft.Container(height=4),
                ft.Row([ft.Text("Nombre",color=C["label"],size=13,width=140,weight=ft.FontWeight.BOLD),
                        ft.Text("Escribir",color=C["label"],size=13,width=185,weight=ft.FontWeight.BOLD),
                        ft.Text("Equivalencia en Pa",color=C["label"],size=13,weight=ft.FontWeight.BOLD)],spacing=8),
                sep(),
                self._unidad_row("Pascal",          "Pa",          "1 Pa (base SI)"),
                self._unidad_row("Torr",            "Torr",        "133.322 Pa ≈ 1 mmHg"),
                self._unidad_row("Atmósfera",       "atm",         "101 325 Pa"),
                self._unidad_row("Milibar",         "mbar",        "100 Pa"),
                self._unidad_row("PSI",             "psi",         "6 894.76 Pa"),
                self._unidad_row("mm de mercurio",  "mmHg",        "133.322 Pa = 1 Torr"),
                self._unidad_row("Pulgadas de Hg",  "inHg",        "3 386.39 Pa"),
                self._unidad_row("Kg por cm²",      "kg/cm², at",  "98 066.5 Pa"),
                nota("💡 Ejemplo: 14.4 Torr = 1 919.8 Pa = 0.01895 atm = 0.2786 psi"),
            ],pad=18),

            # ── Instrucciones por módulo ──────────────────────────────────────
            self._bloque("🔄","Conversión de Presión",[
                "Ingresa el valor numérico en 'Valor'.",
                "Escribe la unidad de origen en 'Unidad' (ver tabla arriba).",
                "Pulsa 'Convertir' → aparece el banner verde y todos los resultados.",
                "Exporta con 'Exportar .txt' → se guarda en tu Escritorio."],
                "760 Torr = 101 325 Pa = 1 atm = 14.696 psi"),

            self._bloque("⚗️","Mezcla de Gases → A y B",[
                "Escribe primero el % y luego el nombre del gas; pulsa 'Agregar'.",
                "Repite para cada gas de la mezcla — se normalizan a 100% automáticamente.",
                "Pulsa 'Calcular A y B' → constantes de Townsend ponderadas.",
                "Usa '→ Paschen ⚡' para enviar A y B al calculador de Paschen.",
                "Usa '→ Simulador 🔬' para enviar la composición al simulador RECI.",
                "Agrega gases personalizados con sus propios A y B en la sección inferior."],
                "Gases predefinidos: H2, N2, CO2, He, Hg, air"),

            self._bloque("⚡","Ley de Paschen — Presión óptima",[
                "Completa Vb, A, B, d y γ (o usa '→ Paschen' desde Mezcla).",
                "'Calcular presión' → presión óptima con verificación cruzada.",
                "'Graficar V vs p·d' → curva de Paschen doble panel.",
                "Panel log: vista completa de toda la curva U.",
                "Panel lineal: zoom ±20× alrededor del mínimo con tu punto de operación.",
                "Exporta .txt (reporte) y .png (gráfica) al Escritorio."],
                "Ejemplo: Vb=400 V, A=12, B=342 (N₂), d=1 cm, γ=0.01 → p ≈ 1.81 Torr"),

            self._bloque("🔬","Simulador Re-Gassing RECI",[
                "Pestaña ⚙️ Configuración: ingresa geometría del tubo, manifold, operación y mezcla.",
                "Pulsa '🔢 Calcular sistema' → mueve a pestaña 📊 Resultados automáticamente.",
                "En 📊 Resultados: volúmenes, masa gas (PV=nRT), Paschen completo, tiempos.",
                "En ▶️ Simulación: panel SCADA, control de válvulas/equipos, botón INICIAR.",
                "Botones de velocidad ×1/×3/×5/×10 controlan la velocidad de la animación.",
                "Los botones de válvulas (V1/V2/V3) simulan el estado físico con validación.",
                "Los equipos (bomba, fuente HV, tubo, tanque) se encienden/montan aquí.",
                "En 📋 Guías: pasos detallados de evacuación y llenado + reglas de seguridad."],
                "La simulación se abre en ventana matplotlib separada con 5 paneles animados."),

            self._bloque("🧪","Gases Disponibles",[
                "Tabla de las constantes A y B de todos los gases (predefinidos + personalizados).",
                "Agrega gases nuevos con nombre, A y B — quedan disponibles en ⚗️ Mezcla."],
                "Gases predefinidos: H2, N2, CO2, He, Hg, air"),

            ft.Container(height=24),
        ],spacing=12)


# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN 8 — APP PRINCIPAL
# ═══════════════════════════════════════════════════════════════════════════════

def main(page: ft.Page):
    page.title         = "Paschen + Simulador Re-Gassing RECI — v12"
    page.bgcolor       = C["bg"]
    page.window_width  = 1140
    page.window_height = 940
    page.window_resizable = True
    page.padding       = 0

    snack = ft.SnackBar(content=ft.Text(""), duration=4000)
    page.overlay.append(snack)
    def show_snack(msg, color=None):
        snack.content = ft.Text(msg, color=color or C["green"], size=13)
        snack.bgcolor = C["surface"]; snack.open = True; page.update()

    # ── Vistas ────────────────────────────────────────────────────────────────
    v_sim  = SimuladorView(show_snack)
    v_pch  = PaschenView(show_snack,
                          on_sim=lambda A,B,gamma,p=None: v_sim.recibir_AB_gamma(A,B,gamma,p))
    v_mez  = MezclaView(show_snack,
                         on_paschen=lambda A,B: v_pch.set_AB(A,B),
                         on_sim=lambda A,B,g: v_sim.recibir_mezcla(A,B,g))
    v_conv = ConversionView(show_snack)
    v_gas  = GasesView(show_snack)
    v_ayu  = AyudaView()

    content = ft.ListView(controls=[], expand=True, spacing=0,
                           padding=ft.Padding.symmetric(horizontal=28,vertical=20))

    NAV = [
        ("📖","Ayuda",      v_ayu),
        ("🔄","Conversión", v_conv),
        ("⚗️", "Mezcla",   v_mez),
        ("⚡", "Paschen",   v_pch),
        ("🔬", "Simulador", v_sim),
        ("🧪", "Gases",     v_gas),
    ]
    nav_btns = []; sel=[0]; dark=[True]

    def switch(idx, pg):
        sel[0]=idx
        for i,b in enumerate(nav_btns):
            b.bgcolor = C["nav_sel"] if i==idx else "transparent"
            b.border  = ft.Border.all(1, C["border_f"] if i==idx else "transparent")
            for ctrl in b.content.controls:
                if isinstance(ctrl,ft.Text): ctrl.color=C["text"]
        _,_,vista = NAV[idx]
        for v in [v_conv,v_mez,v_pch,v_sim,v_gas]:
            if hasattr(v,"refresh_theme"): v.refresh_theme()
        content.controls = [vista.build()]
        if idx==5: v_gas.refrescar(pg)
        pg.update()

    for i,(ic,lb,_) in enumerate(NAV):
        b = ft.Container(
            content=ft.Column([ft.Text(ic,size=22),
                               ft.Text(lb,color=C["text"],size=11,
                                       text_align=ft.TextAlign.CENTER)],
                              horizontal_alignment=ft.CrossAxisAlignment.CENTER,spacing=3),
            width=112, padding=ft.Padding.symmetric(vertical=12,horizontal=4),
            border_radius=12, bgcolor="transparent",
            border=ft.Border.all(1,"transparent"),
            on_click=lambda e,idx=i: switch(idx,e.page), ink=True)
        nav_btns.append(b)

    tema_sw = ft.Switch(value=True, active_color="#38bdf8",
                        inactive_thumb_color="#94a3b8")
    lbl_tema= ft.Text("🌙 Oscuro",color=C["muted"],size=11,text_align=ft.TextAlign.CENTER)
    def on_tema(e):
        dark[0]=tema_sw.value; set_theme(dark[0])
        lbl_tema.value="🌙 Oscuro" if dark[0] else "☀️ Claro"
        lbl_tema.color=C["muted"]; page.bgcolor=C["bg"]
        sidebar.bgcolor=C["nav_bg"]
        sidebar.border=ft.Border.only(right=ft.BorderSide(1,C["border"]))
        for v in [v_conv,v_mez,v_pch,v_sim,v_gas]:
            if hasattr(v,"refresh_theme"): v.refresh_theme()
        switch(sel[0],page)
    tema_sw.on_change=on_tema

    sidebar = ft.Container(
        content=ft.Column([
            ft.Container(height=16),
            ft.Text("PASCHEN",color=C["accent"],size=14,
                    weight=ft.FontWeight.BOLD,text_align=ft.TextAlign.CENTER),
            ft.Text("v12",color=C["muted"],size=10,text_align=ft.TextAlign.CENTER),
            ft.Divider(color=C["border"],height=20),
            *nav_btns,
            ft.Container(expand=True),
            ft.Divider(color=C["border"],height=1),
            ft.Container(
                content=ft.Column([lbl_tema,
                    ft.Row([ft.Text("☀️",size=16),tema_sw,ft.Text("🌙",size=16)],
                           alignment=ft.MainAxisAlignment.CENTER,spacing=2)],
                    horizontal_alignment=ft.CrossAxisAlignment.CENTER,spacing=3),
                padding=ft.Padding.symmetric(vertical=10)),
            ft.Container(height=8),
        ],horizontal_alignment=ft.CrossAxisAlignment.CENTER,spacing=4),
        width=114, bgcolor=C["nav_bg"],
        border=ft.Border.only(right=ft.BorderSide(1,C["border"])))

    page.add(ft.Row([sidebar,content],expand=True,spacing=0))
    switch(0,page)

if __name__ == "__main__":
    ft.run(main)