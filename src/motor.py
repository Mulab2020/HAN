import numpy as np
import numba
from scipy.signal import fftconvolve
import brainpy as bp
import brainpy.math as bm

import matplotlib.pyplot as plt
import itertools
from scipy.stats import ttest_ind

def double_exponential_filter(x, dt, tau_rise, tau_decay):
    """
    Double-exponential filtering implemented via FFT convolution (causal, non-negative), with area normalization.
    """
    T_window = 5 * tau_decay
    L = max(1, int(T_window / dt))
    t_kernel = np.arange(L) * dt

    if tau_rise == tau_decay:
        w = (t_kernel / tau_decay) * np.exp(-t_kernel / tau_decay)
    else:
        norm = 1.0 / (tau_decay - tau_rise)
        w = norm * (np.exp(-t_kernel / tau_decay) - np.exp(-t_kernel / tau_rise))

    s = w.sum()
    if s > 0:
        w /= s

    y = fftconvolve(x, w, mode='full')[:len(x)]
    return y

@numba.jit(nopython=True)
def gated_trace_fast(x, theta, alpha, gain=1.0, max_val=None):
    g = np.zeros_like(x)
    rel_all = gain * np.maximum(x - theta, 0.0)

    current_g = 0.0
    for t in range(1, len(x)):
        current_g = (1.0 - alpha) * current_g + alpha * rel_all[t]
        if max_val is not None and current_g > max_val:
            current_g = max_val
        g[t] = current_g
    return g

@numba.jit(nopython=True)

def poisson_noise_rate_and_std_modulated(modulator, dt,
                                         base_rate_hz=1.0,
                                         rate_gain_hz=5.0,
                                         rate_clip_hz=200.0,
                                         base_std=0.05,
                                         std_gain=0.5,
                                         std_clip=5.0):
    """
    modulator: 1D array, used to modulate intensity (recommended: pass signal = log1p(...))
    dt: ms

    rate(t) = base_rate_hz + rate_gain_hz * modulator(t)
    std(t)  = base_std     + std_gain     * modulator(t)

    At each time step, an event is triggered with prob=rate(t)*dt/1000; on trigger, add a N(0, std(t)) pulse.
    Return abs(noise) to keep it non-negative (matches previous logic).
    """
    n = modulator.shape[0]
    noise = np.zeros(n)

    for i in range(n):
        # --- rate modulation ---
        rate_hz = base_rate_hz + rate_gain_hz * modulator[i]
        if rate_hz < 0.0:
            rate_hz = 0.0
        if rate_hz > rate_clip_hz:
            rate_hz = rate_clip_hz

        prob = rate_hz * dt / 1000.0
        if np.random.random() < prob:
            # --- std modulation ---
            std_t = base_std + std_gain * modulator[i]
            if std_t < 0.0:
                std_t = 0.0
            if std_t > std_clip:
                std_t = std_clip

            noise[i] = np.random.normal(0.0, std_t)

    return np.abs(noise)


def process_motor_readout(runner, params=None):
    p = {
        'inp_thresh': 1.0,
        'inp_gain': 1.0,
        'tau_gate': 500.0,
        'max_gate': 500.0,
        'tau_rise': 100.0,   # rise time constant
        'tau_decay': 3000.0, # previous tau_hb becomes the decay term

        # rate modulation
        'noise_base_hz': 0.5,
        'noise_rate_gain_hz': 1.0,
        'noise_rate_clip_hz': 3.0,

        # std modulation
        'noise_std_base': 0.1,
        'noise_std_gain': 0.1,
        'noise_std_clip': 3.0,
        }
    if params:
        p.update(params)

    ts = runner.mon['ts']
    dt = float(np.mean(np.diff(ts)))  # ms

    # firing rate
    Inp_A_rate = bp.measure.firing_rate(runner.mon['Inp_A.spike'], width=30.)
    Inp_B_rate = bp.measure.firing_rate(runner.mon['Inp_B.spike'], width=30.)
    hb_A_rate  = bp.measure.firing_rate(runner.mon['hb_A.spike'],  width=30.)
    hb_B_rate  = bp.measure.firing_rate(runner.mon['hb_B.spike'],  width=30.)

    # A) gate
    alpha_gate = dt / p['tau_gate']
    gate_A = gated_trace_fast(Inp_B_rate, p['inp_thresh'], alpha_gate, p['inp_gain'], p['max_gate'])
    gate_B = gated_trace_fast(Inp_A_rate, p['inp_thresh'], alpha_gate, p['inp_gain'], p['max_gate'])

    # B) hb integration
    hb_A_int = double_exponential_filter(hb_A_rate, dt, p['tau_rise'], p['tau_decay'])
    hb_B_int = double_exponential_filter(hb_B_rate, dt, p['tau_rise'], p['tau_decay'])

    # C) deterministic signal (this will control poisson strength)
    signal_A = np.log1p(gate_A * (1.0 + 1.0 * hb_A_int))
    signal_B = np.log1p(gate_B * (1.0 + 1.0 * hb_B_int))

    # deterministic signal
    signal_A = np.log1p(gate_A * (1.0 + 1.0 * hb_A_int))
    signal_B = np.log1p(gate_B * (1.0 + 1. * hb_B_int))

    # poisson noise with std controlled by signal
    noise_A = poisson_noise_rate_and_std_modulated(
        signal_A, dt,
        base_rate_hz=p['noise_base_hz'],
        rate_gain_hz=p['noise_rate_gain_hz'],
        rate_clip_hz=p['noise_rate_clip_hz'],
        base_std=p['noise_std_base'],
        std_gain=p['noise_std_gain'],
        std_clip=p['noise_std_clip'],
    )

    noise_B = poisson_noise_rate_and_std_modulated(
        signal_B, dt,
        base_rate_hz=p['noise_base_hz'],
        rate_gain_hz=p['noise_rate_gain_hz'],
        rate_clip_hz=p['noise_rate_clip_hz'],
        base_std=p['noise_std_base'],
        std_gain=p['noise_std_gain'],
        std_clip=p['noise_std_clip'],
    )

    motor_A = noise_A
    motor_B = noise_B


    return ts, motor_A, motor_B


def plot_final_results_separate(ts, motor_A, motor_B,
                                savepath="motor_output_example.pdf",
                                show=False):
    import matplotlib.pyplot as plt

    ts = np.asarray(ts)
    motor_A = np.asarray(motor_A)
    motor_B = np.asarray(motor_B)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True)

    # --- Motor A ---
    axes[0].plot(ts, motor_A, alpha=0.8)
    axes[0].set_ylabel('Motor Left')
    axes[0].grid(True, alpha=0.3)

    # --- Motor B ---
    axes[1].plot(ts, motor_B, alpha=0.8)
    axes[1].set_ylabel('Motor Right')
    axes[1].grid(True, alpha=0.3)

    # --- Decision variable ---
    diff = motor_A - motor_B
    decision = np.cumsum(diff)
    axes[2].plot(ts, decision, label='Cumulative Decision', alpha=0.9)
    axes[2].axhline(0, linestyle='--', linewidth=1)
    axes[2].set_xlabel('Time (ms)')
    axes[2].set_ylabel('Avoiding Distance')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(savepath, format="pdf", bbox_inches="tight", transparent=True)
    if show:
        plt.show()
    plt.close(fig)


# load saved average rates and generate motor outputs
# ----------------- motor readout with signal+noise options -----------------

# ---------- motor readout ----------
def process_motor_readout_final(average_rates, ts, params=None, output_mode="noise"):
    """
    output_mode:
      - "noise": motor_A/B = noise_A/B  (current usage)
      - "signal+noise": motor_A/B = signal_A/B + noise_A/B
      - "signal": motor_A/B = signal_A/B
    """
    p = {
        'inp_thresh': 1.0,
        'inp_gain': 1.0,
        'tau_gate': 500.0,
        'max_gate': 500.0,
        'tau_rise': 100.0,
        'tau_decay': 3000.0,

        # rate modulation
        'noise_base_hz': 0.5,
        'noise_rate_gain_hz': 1.0,
        'noise_rate_clip_hz': 3.0,

        # std modulation
        'noise_std_base': 0.1,
        'noise_std_gain': 0.1,
        'noise_std_clip': 3.0,

        # signal
        'hb_gain': 1.0,     # the 1.0 in your (1.0 + 1.0 * hb_int)
    }
    if params:
        p.update(params)

    ts = np.asarray(ts)
    dt = float(np.mean(np.diff(ts)))  # if ts is in ms, dt is in ms here

    Inp_A_rate = average_rates['Inp_A']
    Inp_B_rate = average_rates['Inp_B']
    hb_A_rate  = average_rates['hb_A']
    hb_B_rate  = average_rates['hb_B']

    alpha_gate = dt / p['tau_gate']
    gate_A = gated_trace_fast(Inp_B_rate, p['inp_thresh'], alpha_gate, p['inp_gain'], p['max_gate'])
    gate_B = gated_trace_fast(Inp_A_rate, p['inp_thresh'], alpha_gate, p['inp_gain'], p['max_gate'])

    hb_A_int = double_exponential_filter(hb_A_rate, dt, p['tau_rise'], p['tau_decay'])
    hb_B_int = double_exponential_filter(hb_B_rate, dt, p['tau_rise'], p['tau_decay'])

    signal_A = np.log1p(gate_A * (1.0 + p['hb_gain'] * hb_A_int))
    signal_B = np.log1p(gate_B * (1.0 + p['hb_gain'] * hb_B_int))

    noise_A = poisson_noise_rate_and_std_modulated(
        signal_A, dt,
        base_rate_hz=p['noise_base_hz'],
        rate_gain_hz=p['noise_rate_gain_hz'],
        rate_clip_hz=p['noise_rate_clip_hz'],
        base_std=p['noise_std_base'],
        std_gain=p['noise_std_gain'],
        std_clip=p['noise_std_clip'],
    )
    noise_B = poisson_noise_rate_and_std_modulated(
        signal_B, dt,
        base_rate_hz=p['noise_base_hz'],
        rate_gain_hz=p['noise_rate_gain_hz'],
        rate_clip_hz=p['noise_rate_clip_hz'],
        base_std=p['noise_std_base'],
        std_gain=p['noise_std_gain'],
        std_clip=p['noise_std_clip'],
    )

    if output_mode == "noise":
        motor_A, motor_B = noise_A, noise_B
    elif output_mode == "signal+noise":
        motor_A, motor_B = signal_A + noise_A, signal_B + noise_B
    elif output_mode == "signal":
        motor_A, motor_B = signal_A, signal_B
    else:
        raise ValueError(f"Unknown output_mode: {output_mode}")

    return ts, motor_A, motor_B


# ---------- plotting helpers ----------
def _collect_trial_traces(motor_diff, seqs, one_trial, n_trials_in_key=2):
    """
    Group by consecutive trials:
      n_trials_in_key=2: key = seqs[i] + seqs[i+1] ，“ trial (i+1)”
      n_trials_in_key=3: key = seqs[i-1]+seqs[i]+seqs[i+1]，“ trial (i+1)”
    Returns: trial_rates dict[key] -> (n_samples_in_group, one_trial)
    """
    seqs = np.asarray(seqs)

    trial_rates = {}
    nT = len(seqs)

    if n_trials_in_key == 2:
        # i from 1 .. nT-2  ( i+1 ，1)
        for i in range(1, nT - 1):
            key = str(seqs[i]) + str(seqs[i + 1])
            start = (i + 1) * one_trial
            end   = (i + 2) * one_trial
            if end > len(motor_diff):
                break
            trial_rates.setdefault(key, []).append(np.cumsum(motor_diff[start:end]))

    elif n_trials_in_key == 3:
        # i from 1 .. nT-2，key  (i-1,i,i+1)， trial (i+1)
        for i in range(1, nT - 1):
            key = str(seqs[i - 1]) + str(seqs[i]) + str(seqs[i + 1])
            start = (i + 1) * one_trial
            end   = (i + 2) * one_trial
            if end > len(motor_diff):
                break
            trial_rates.setdefault(key, []).append(np.cumsum(motor_diff[start:end]))

    else:
        raise ValueError("n_trials_in_key must be 2 or 3")

    # list -> np.array
    for k in list(trial_rates.keys()):
        trial_rates[k] = np.asarray(trial_rates[k])
        if trial_rates[k].size == 0:
            trial_rates.pop(k, None)

    return trial_rates


def merge_trial_rates(trial_rates_list):
    """
    trial_rates_list: list of dict[key] -> (n_trials, one_trial)
    return: merged dict[key] -> (sum_n_trials, one_trial)
    """
    merged = {}
    for d in trial_rates_list:
        for k, arr in d.items():
            if arr is None or len(arr) == 0:
                continue
            merged.setdefault(k, []).append(arr)
    for k in list(merged.keys()):
        merged[k] = np.concatenate(merged[k], axis=0)
    return merged


def collect_one_batch_trial_rates(data, runner, tool, n_trials_in_key=2, output_mode="noise", params=None):
    """
     batch  data/runner， trial_rates dict[key] -> (n_trials, one_trial)
    Reuses existing functions: process_motor_readout_final + _collect_trial_traces
    """
    average_rates = data['average_rates'].item()
    seqs = data['sequences']

    ts, mA, mB = process_motor_readout_final(
        average_rates,
        ts=runner.mon['ts'],
        params=params,
        output_mode=output_mode,
    )

    one_trial = int((tool.pre_stimulus_period + tool.stimulus_period + tool.delay_period) / bm.get_dt())
    motor_diff = mA - mB
    trial_rates = _collect_trial_traces(motor_diff, seqs, one_trial, n_trials_in_key=n_trials_in_key)
    return ts, one_trial, trial_rates

def _p_to_stars(p):
    if not np.isfinite(p):
        return "n/a"
    if p < 1e-4: return "****"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 0.05: return "*"
    return "ns"

def _holm_bonferroni(pvals):
    """
    Holm-Bonferroni correction (step-down).
    Returns adjusted p-values in original order.
    """
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]

    # raw Holm adjustment
    raw = (m - np.arange(m)) * ranked

    # enforce monotonic non-decreasing in sorted order
    adj_sorted = np.maximum.accumulate(raw)
    adj_sorted = np.clip(adj_sorted, 0, 1)

    adj = np.empty(m, dtype=float)
    adj[order] = adj_sorted
    return adj

def _pairwise_welch_holm(per_key_vals, keys):
    """
    per_key_vals: dict[key] -> 1D array of per-trial metric
    keys: list like ['RR','LR','RL','LL']
    returns: list of rows dict with p_holm + stars
    """
    rows = []
    for a, b in itertools.combinations(keys, 2):
        va = per_key_vals.get(a, None)
        vb = per_key_vals.get(b, None)
        if va is None or vb is None or len(va) < 2 or len(vb) < 2:
            t, p = np.nan, np.nan
            na, nb = (0 if va is None else len(va)), (0 if vb is None else len(vb))
        else:
            t, p = ttest_ind(va, vb, equal_var=False, nan_policy="omit")  # Welch
            na, nb = len(va), len(vb)
        rows.append({"A": a, "B": b, "t": t, "p": p, "nA": na, "nB": nb})

    pvals = np.array([r["p"] for r in rows], dtype=float)
    valid = np.isfinite(pvals)
    padj = np.full_like(pvals, np.nan, dtype=float)
    if valid.sum() > 0:
        padj[valid] = _holm_bonferroni(pvals[valid])

    for r, pa in zip(rows, padj):
        r["p_holm"] = pa
        r["sig_holm"] = _p_to_stars(pa)
    return rows

def _stats_text_block(rows, title="Pairwise (Holm)"):
    lines = [title]
    for r in rows:
        # e.g., RR vs LL: ** (p=0.0031)
        if np.isfinite(r["p_holm"]):
            lines.append(f"{r['A']} vs {r['B']}: {r['sig_holm']} (p={r['p_holm']:.2g})")
        else:
            lines.append(f"{r['A']} vs {r['B']}: n/a")
    return "\n".join(lines)

def plot_motor_results_from_trial_rates(ts, one_trial, trial_rates, tool,
                                       n_trials_in_key=2,
                                       max_plot_ms=12000,
                                       colors=None,
                                       prefix="motor_multi",
                                       last_window_ms=1000,
                                       downsample_step=50):
    """
    Fig1: time course mean±std of cumsum(motor_A - motor_B)
    Fig2: final location bar (mean±std) + pairwise Welch t-test (Holm) among 4 conditions
    Fig3: last_window_ms net cumulative bar (mean±std) + pairwise Welch t-test (Holm) among 4 conditions

    trial_rates[key] shape = (n_trials, one_trial), values are cumsum(motor_A - motor_B).
    """
    if colors is None and n_trials_in_key == 2:
        colors = {'RR':'#56ac55','LR':'#87ad7e','RL':'#dd9eb3','LL':'#c169a5'}

    # keys order
    if n_trials_in_key == 2:
        keys = ['RR', 'LR', 'RL', 'LL']
    else:
        keys = ["RRR", "LRR", "RLR", "LLR",
                "RRL", "LRL", "RLL", "LLL"]

    ts_trial = np.asarray(ts[:one_trial])
    dt_ms = float(np.mean(np.diff(ts_trial)))  # ms

    if max_plot_ms is not None:
        mask = ts_trial <= max_plot_ms
        ts_plot = ts_trial[mask]
    else:
        mask = slice(None)
        ts_plot = ts_trial

    # ---------------- Fig1: time course ----------------
    plt.figure(figsize=(10, 8))
    for key in keys:
        if key not in trial_rates:
            continue

        y = np.mean(trial_rates[key], axis=0)[mask]
        n = trial_rates[key].shape[0]
        std = np.std(trial_rates[key], axis=0, ddof=1)[mask]

        ts_ds  = ts_plot[::downsample_step]
        y_ds   = y[::downsample_step]
        std_ds = std[::downsample_step]

        if colors is not None and key in colors:
            c = colors[key]
            plt.plot(ts_ds, y_ds, label=f"{key} (n={n})", color=c)
            plt.fill_between(ts_ds, y_ds - std_ds, y_ds + std_ds,
                             color=c, alpha=0.2, edgecolor="none", linewidth=0)
        else:
            plt.plot(ts_ds, y_ds, label=f"{key} (n={n})")
            plt.fill_between(ts_ds, y_ds - std_ds, y_ds + std_ds,
                             alpha=0.2, edgecolor="none", linewidth=0)

    for tline in [tool.pre_stimulus_period, tool.pre_stimulus_period + tool.stimulus_period]:
        if max_plot_ms is None or tline <= max_plot_ms:
            plt.axvline(tline, linestyle='dashed')

    plt.legend(ncol=2, fontsize=9)
    plt.title(f'Motor Output (group by {n_trials_in_key}-trial key)')
    plt.xlabel('Time (ms)')
    plt.ylabel('Avoiding Distance (a.u.)')
    plt.tight_layout()
    plt.savefig(f'{prefix}_{n_trials_in_key}trial_cumsum_upto_{max_plot_ms}ms.pdf')

    # helper for bar colors
    x = np.arange(len(keys))
    bar_colors = [colors[k] for k in keys] if (colors is not None and all(k in colors for k in keys)) else None

    # ---------------- Fig2: final location bar + pairwise stats ----------------
    final_means = []
    final_stds  = []
    final_ns    = []
    per_key_final = {}

    for k in keys:
        if k in trial_rates and len(trial_rates[k]) > 0:
            v = trial_rates[k][:, -1]  # per-trial final cumsum
            per_key_final[k] = v
            final_means.append(np.mean(v))
            final_stds.append(np.std(v, ddof=1))
            final_ns.append(len(v))
        else:
            final_means.append(np.nan)
            final_stds.append(np.nan)
            final_ns.append(0)

    plt.figure(figsize=(max(8, 0.35 * len(keys)), 6))
    ax2 = plt.gca()
    ax2.bar(x, final_means, yerr=final_stds, color=bar_colors)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{k}\n(n={final_ns[i]})" for i, k in enumerate(keys)], rotation=45, ha='right')
    ax2.set_title(f'Final distance (group by {n_trials_in_key}-trial key)')
    ax2.set_ylabel('Avoiding distance (a.u.)')

    # pairwise only for 4-condition case
    if len(keys) == 4:
        rows2 = _pairwise_welch_holm(per_key_final, keys)

        # print
        print(f"\n[Pairwise stats | Fig2] {prefix} | final cumsum | Welch + Holm")
        for r in rows2:
            print(f"  {r['A']} vs {r['B']}: t={r['t']:.3g}, p={r['p']:.3g}, "
                  f"p_holm={r['p_holm']:.3g}, {r['sig_holm']} (n={r['nA']},{r['nB']})")

        # save CSV
        try:
            import pandas as pd
            df2 = pd.DataFrame(rows2, columns=["A","B","nA","nB","t","p","p_holm","sig_holm"])
            csv2 = f"{prefix}_{n_trials_in_key}trial_Fig2_final_pairwise_welch_holm.csv"
            df2.to_csv(csv2, index=False)
            print(f"[Saved] {csv2}")
        except Exception as e:
            print("[Warn] Fig2: failed to save CSV (pandas missing?):", e)

        # add textbox
        ax2.text(1.02, 0.98, _stats_text_block(rows2, title="Pairwise (Holm):"),
                 transform=ax2.transAxes, va="top", ha="left", fontsize=9)

    plt.tight_layout()
    plt.savefig(f'{prefix}_{n_trials_in_key}trial_final_location.pdf')

    # ---------------- Fig3: last 1s net cumulative bar + pairwise stats ----------------
    w = int(round(last_window_ms / dt_ms))
    w = max(1, min(w, one_trial - 1))

    last_means = []
    last_stds  = []
    last_ns    = []
    per_key_last = {}

    for k in keys:
        if k in trial_rates and len(trial_rates[k]) > 0:
            # last window net accumulation: Δcumsum over last window
            v = trial_rates[k][:, -1] - trial_rates[k][:, -1 - w]
            per_key_last[k] = v
            last_means.append(np.mean(v))
            last_stds.append(np.std(v, ddof=1))
            last_ns.append(len(v))
        else:
            last_means.append(np.nan)
            last_stds.append(np.nan)
            last_ns.append(0)

    plt.figure(figsize=(max(8, 0.35 * len(keys)), 6))
    ax3 = plt.gca()
    ax3.bar(x, last_means, yerr=last_stds, color=bar_colors)
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"{k}\n(n={last_ns[i]})" for i, k in enumerate(keys)], rotation=45, ha='right')
    ax3.set_title(f'Last {last_window_ms} ms net cumulative (MERGED, group by {n_trials_in_key}-trial key)')
    ax3.set_ylabel(f'Δcumsum over last {last_window_ms} ms')

    if len(keys) == 4:
        rows3 = _pairwise_welch_holm(per_key_last, keys)

        print(f"\n[Pairwise stats | Fig3] {prefix} | last{last_window_ms}ms Δcumsum | Welch + Holm")
        for r in rows3:
            print(f"  {r['A']} vs {r['B']}: t={r['t']:.3g}, p={r['p']:.3g}, "
                  f"p_holm={r['p_holm']:.3g}, {r['sig_holm']} (n={r['nA']},{r['nB']})")

        try:
            import pandas as pd
            df3 = pd.DataFrame(rows3, columns=["A","B","nA","nB","t","p","p_holm","sig_holm"])
            csv3 = f"{prefix}_{n_trials_in_key}trial_Fig3_last{last_window_ms}ms_pairwise_welch_holm.csv"
            df3.to_csv(csv3, index=False)
            print(f"[Saved] {csv3}")
        except Exception as e:
            print("[Warn] Fig3: failed to save CSV (pandas missing?):", e)

        ax3.text(1.02, 0.98, _stats_text_block(rows3, title="Pairwise (Holm):"),
                 transform=ax3.transAxes, va="top", ha="left", fontsize=9)

    plt.tight_layout()
    plt.savefig(f'{prefix}_{n_trials_in_key}trial_last{last_window_ms}ms_location.pdf')
