import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

import numpy as np

import brainpy as bp
import brainpy.math as bm

import tqdm

import datetime
date_str = datetime.datetime.now().strftime("%Y%m%d")


class Tool:
  def __init__(self, seq):
    self.pre_stimulus_period = 1000.
    self.stimulus_period = 1_000*10.
    self.delay_period = 1_000*24.
    self.freq_variance = 1.
    self.freq_interval = 20.
    self.seq = seq
    self.total_period = len(seq) * (self.pre_stimulus_period
                                    + self.stimulus_period + self.delay_period)

  def generate_freqs(self, R, L):
    # stimulus period
    n_stim = int(self.stimulus_period / self.freq_interval)
    n_interval = int(self.freq_interval / bm.get_dt())

    # pre stimulus period
    freqs_pre = np.zeros(int(self.pre_stimulus_period / bm.get_dt()))
    # post stimulus period
    freqs_delay = np.zeros(int(self.delay_period / bm.get_dt()))

    # ---- You can expose these as self.tau_rise / self.tau_decay / self.rise_frac for tuning ----
    rise_frac = getattr(self, "rise_frac", 0.4)  # matches your previous n_stim//8
    tau_rise  = getattr(self, "tau_rise",  self.stimulus_period * 0.3)  # default: faster rise
    tau_decay = getattr(self, "tau_decay", self.stimulus_period * 0.8)  # default: slower decay

    # rise length / decay length (at least 2 points to avoid divide-by-zero / too short)
    n_rise = max(2, int(round(n_stim * rise_frac)))
    n_rise = min(n_rise, n_stim - 2)  # leave some length for the decay segment
    n_decay = n_stim - n_rise

    # time interval per 'freq step' (not bm.get_dt(), but freq_interval)
    dt_step = self.freq_interval

    def exp_rise(mean, n, tau):
      T = (n - 1) * dt_step
      t = np.arange(n) * dt_step
      # normalize so the last point equals mean
      denom = (1.0 - np.exp(-T / tau)) if T > 0 else 1.0
      y = mean * (1.0 - np.exp(-t / tau)) / denom
      return y.reshape(-1, 1)

    def exp_decay_to_zero(mean, n, tau):
      T = (n - 1) * dt_step
      t = np.arange(n) * dt_step
      # normalize so the first point equals mean and the last point equals 0
      a = np.exp(-t / tau)
      aT = np.exp(-T / tau) if T > 0 else 0.0
      denom = (1.0 - aT) if T > 0 else 1.0
      y = mean * (a - aT) / denom
      return y.reshape(-1, 1)

    for i in range(len(self.seq)):
      mean = R if self.seq[i] == 'R' else L

      freqs_stim1 = exp_rise(mean, n_rise, tau_rise) \
                    + np.random.normal(0, self.freq_variance, (n_rise, 1))
      freqs_stim2 = exp_decay_to_zero(mean, n_decay, tau_decay) \
                    + np.random.normal(0, self.freq_variance, (n_decay, 1))

      freqs_stim = np.concatenate([freqs_stim1, freqs_stim2], axis=0)

      # original logic: interpolate each freq_interval to bm.get_dt() resolution
      freqs_stim = np.tile(freqs_stim, (1, n_interval)).flatten()
      one_freqs = np.concatenate([freqs_pre, freqs_stim, freqs_delay], axis=0)
      all_freqs = one_freqs if i == 0 else np.concatenate([all_freqs, one_freqs], axis=0)
    return bm.asarray(all_freqs)

  def visualize_results(self, mon, monlist, params, t_start=0., title=None, index=0, savefig=False, filename=None):
      fig, axes = plt.subplots(3, 1, figsize=(12, 10))
      mon_ts = mon['ts']
      rates = {key: bp.measure.firing_rate(mon[key], width=30.) for key in monlist}

      # subfigure1
      ax = axes[0]
      plot_keys = [monlist[0], monlist[1]]
      for key in plot_keys:
          ax.plot(mon_ts, rates[key], label=key)
      ax.set_ylabel('Population activity in layer 1 [Hz]')
      ax.set_xlim(t_start, self.total_period + 1)
      self._add_vertical_lines(ax)
      ax.legend(loc='upper right')

      # subfigure2
      ax = axes[1]
      plot_keys = [monlist[2], monlist[3], monlist[4], monlist[5], monlist[6], monlist[7]]
      for key in plot_keys:
          ax.plot(mon_ts, rates[key], label=key)
      ax.set_ylabel('Population activity  in layer 2 [Hz]')
      ax.set_xlim(t_start, self.total_period + 1)
      self._add_vertical_lines(ax)
      ax.legend(loc='upper right')

      # subfigure3
      ax = axes[2]
      plot_keys = [monlist[8], monlist[9], monlist[10], monlist[11]]
      for key in plot_keys:
          ax.plot(mon_ts, rates[key], label=key)
      ax.set_ylabel('Population activity  in layer 2 [Hz]')
      ax.set_xlim(t_start, self.total_period + 1)
      self._add_vertical_lines(ax)
      ax.legend(loc='upper right')
      ax.set_xlabel("Time [ms]")

      plt.tight_layout(rect=[0, 0.15, 1, 0.95])

      # print parameters
      parameters_text = (
          f"Parameters:\n"
          f"Pre-stimulus period: {self.pre_stimulus_period} ms, "
          f"Stimulus period: {self.stimulus_period} ms, "
          f"Delay period: {self.delay_period} ms, "
          f"Frequency variance: {self.freq_variance}, "
          f"Frequency interval: {self.freq_interval} ms, "
          f"Sequence: {self.seq}\n"
          "Model Parameters:\n" +
          ', '.join([f"{key}: {value}" for key, value in params.items()])
      )
      plt.figtext(0.5, 0.01, parameters_text, wrap=True, horizontalalignment='center', fontsize=8)

      # save in PDF format
      if savefig:
          if filename is None:
              filename = f'visualize_results_{index}.pdf'
          plt.savefig(filename, format='pdf')

  def _add_vertical_lines(self, ax):
      trial_period = self.pre_stimulus_period + self.stimulus_period + self.delay_period
      for i in range(len(self.seq)):
          start = self.pre_stimulus_period + i * trial_period
          ax.axvline(start, linestyle='dashed')
          ax.axvline(start + self.stimulus_period, linestyle='dashed')

  def visualize_current_bar_chart(self, currents, labels, mon_ts, plot_cond='R', title=None, savefig=False, filename=None):
        """
        Plot a bar chart of currents for the 2nd L trial under each condition, showing stimulus and delay
        mean currents with error bars (STD), arranged as 4 subplots with LL and RL on separate rows.

        Args:
        currents: List[np.ndarray]
            List of current arrays; each element is [num_time_steps, num_neurons] or [num_time_steps].
        labels: List[str]
            List of labels for each current trace (legend).
        mon_ts: np.ndarray
            Time vector matching the current traces.
        title: str, optional
            Overall figure title.
        savefig: bool, optional
            Whether to save as PDF (default: False).
        filename: str, optional
            Output filename (default: 'current_bar_chart.pdf').
        """
        # Validate input consistency
        if not currents or not labels or mon_ts is None:
            raise ValueError("currents, labels, and mon_ts must be provided and not empty.")

        #  trial
        trial_period = self.pre_stimulus_period + self.stimulus_period + self.delay_period

        #  'LL'  'RL'  L trial
        selected_trials = {}
        if plot_cond == 'R':
            conditions = ['RR', 'LR']
        else:
            conditions = ['LL', 'RL']
        for condition in conditions:
            count = 0
            for i in range(1, len(self.seq)):
                # ，， 'LL'  'RL'
                trial = self.seq[i-1] + self.seq[i]
                if trial == condition:
                    count += 1
                    if count == 2:
                        selected_trials[condition] = i
                        break
            if condition not in selected_trials:
                raise ValueError(f" {condition}  L trial。")

        # ， trial  stimulus  delay
        results = {}
        for condition, trial_idx in selected_trials.items():
            trial_start = self.pre_stimulus_period + trial_idx * trial_period
            stim_start = trial_start + self.stimulus_period//3
            stim_end = stim_start + self.stimulus_period//3 # 1/3 --> 2/3*stimulus_period # 200ms
            delay_start = stim_end + self.delay_period - self.stimulus_period//3 # 200ms
            delay_end = stim_end + self.delay_period

            #  mon_ts
            stim_mask = (mon_ts >= stim_start) & (mon_ts < stim_end)
            delay_mask = (mon_ts >= delay_start) & (mon_ts < delay_end)

            stim_avgs = []
            stim_stds = []
            delay_avgs = []
            delay_stds = []

            for current in currents:
                #  current ，
                if current.ndim == 2:
                    mean_current = current.mean(axis=1)
                elif current.ndim == 1:
                    mean_current = current

                # Stimulus
                stim_values = mean_current[stim_mask]
                stim_avg = np.mean(stim_values)
                stim_std = np.std(stim_values)

                # Delay
                delay_values = mean_current[delay_mask]
                delay_avg = np.mean(delay_values)
                delay_std = np.std(delay_values)

                stim_avgs.append(stim_avg)
                stim_stds.append(stim_std)
                delay_avgs.append(delay_avg)
                delay_stds.append(delay_std)

            results[condition] = {
                "stim_mean": stim_avgs,
                "stim_std": stim_stds,
                "delay_mean": delay_avgs,
                "delay_std": delay_stds
            }

        #  4 ：2  x 2
        fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharey=True)

        #
        bar_width = 0.8  #  Stimulus  Delay ，
        stim_color = 'skyblue'
        delay_color = 'orange'
        x = np.arange(len(labels))

        #  LL ( 1 )  RL ( 2 )
        for i, condition in enumerate(conditions):
            stim_avgs = results[condition]["stim_mean"]
            stim_stds = results[condition]["stim_std"]
            delay_avgs = results[condition]["delay_mean"]
            delay_stds = results[condition]["delay_std"]

            # ========== (i, 0) Stimulus subplot ==========
            ax_stim = axs[i, 0]
            bars_stim = ax_stim.bar(
                x,
                stim_avgs,
                bar_width,
                yerr=stim_stds,
                capsize=5,
                color=stim_color,
                edgecolor='black'
            )
            ax_stim.set_xticks(x)
            ax_stim.set_xticklabels(labels, rotation=15)
            ax_stim.axhline(0, color='black', linewidth=1)
            ax_stim.set_title(f'Condition {condition} - Stimulus')
            # ，
            ax_stim.bar_label(
                bars_stim,
                fmt='%.2f',
                label_type='edge',
                color='black',  #  stim_color，
                rotation=90,
                padding=3
            )

            # ========== (i, 1) Delay subplot ==========
            ax_delay = axs[i, 1]
            bars_delay = ax_delay.bar(
                x,
                delay_avgs,
                bar_width,
                yerr=delay_stds,
                capsize=5,
                color=delay_color,
                edgecolor='black'
            )
            ax_delay.set_xticks(x)
            ax_delay.set_xticklabels(labels, rotation=15)
            ax_delay.axhline(0, color='black', linewidth=1)
            ax_delay.set_title(f'Condition {condition} - Delay')
            #
            ax_delay.bar_label(
                bars_delay,
                fmt='%.2f',
                label_type='edge',
                color='black',  # ， delay_color
                rotation=90,
                padding=3
            )

        #  y
        axs[0, 0].set_ylabel('Average Current')
        axs[1, 0].set_ylabel('Average Current')

        if title:
            fig.suptitle(title, fontsize=14)

        plt.tight_layout(rect=[0, 0, 1, 0.95])

        #
        if savefig:
            if filename is None:
                filename = 'current_bar_chart.pdf'
            plt.savefig(filename, format='pdf')
        else:
            plt.show()

def calculate_firing_rate(names, runner):
    num_neurons = runner.mon[names[0] + '.spike'].shape[1]
    rates = []
    for i in tqdm.tqdm(range(num_neurons), desc=f"Calculating firing rates for {names[0]}"):
        rates.append(bp.measure.firing_rate(runner.mon[names[0] + '.spike'][:, i].reshape(-1, 1), width=30.))
    return np.asarray(rates)

def average_with_timebins(rate, time_bin):
    num_bins = rate.shape[1] // time_bin
    rate_bin = rate.reshape(rate.shape[0], num_bins, time_bin)
    rate_bin = np.mean(rate_bin, axis=2)
    return rate_bin


def add_trial_lines(ax, trial_len, pre_stimulus_period, stimulus_period, delay_period, dt):
    points = [
        pre_stimulus_period,
        pre_stimulus_period + stimulus_period,
        pre_stimulus_period + stimulus_period + delay_period
    ]
    for p in points:
        ax.axvline(p, linestyle='--', color='lightgray')
        ax.axvline(p + trial_len * dt, linestyle='--', color='lightgray')

def fill_stimulus(ax, trial_len, pre_stimulus_period, stimulus_period, dt):
    ymin, ymax = ax.get_ylim()
    for shift in [0, trial_len * dt]:
        ax.fill_betweenx([ymin, ymax],
                         shift + pre_stimulus_period,
                         shift + pre_stimulus_period + stimulus_period,
                         color='gray', alpha=0.2)

def compute_y_range(traj_dict, unique_conds):
    all_vals = np.concatenate(
        [np.array(traj_dict[c]) for c in unique_conds if len(traj_dict[c]) > 0],
        axis=0
    )
    ymin, ymax = np.min(all_vals), np.max(all_vals)
    margin = 0.01 * (ymax - ymin) if ymax > ymin else 1.0
    return ymin - margin, ymax + margin

def compute_y_range_nan_safe(traj_dict, unique_conds):
    vals = []
    for c in unique_conds:
        if len(traj_dict[c]) == 0:
            continue
        arr = np.asarray(traj_dict[c])
        if arr.size == 0:
            continue
        arr = arr[np.isfinite(arr)]
        if arr.size > 0:
            vals.append(arr)
    if len(vals) == 0:
        return (-1, 1)  # fallback
    all_vals = np.concatenate(vals)
    ymin, ymax = np.nanmin(all_vals), np.nanmax(all_vals)
    if not np.isfinite(ymin) or not np.isfinite(ymax) or ymax == ymin:
        return (ymin - 1, ymax + 1) if np.isfinite(ymin) else (-1, 1)
    margin = 0.01 * (ymax - ymin)
    return ymin - margin, ymax + margin


def plot_avg_trajectory_with_std(
    time,
    trial_len,
    data_dict,
    unique_conds,
    color_map,
    ylims,
    title,
    ylabel,
    filename_base,
    pre_stimulus_period,
    stimulus_period,
    delay_period,
    dt,
    downsample_step=5
):
    def _draw(ax, xlim=None):
        time_ds = time[::downsample_step]

        for cond in unique_conds:
            if len(data_dict[cond]) == 0:
                continue
            data = np.asarray(data_dict[cond])  # (n_trials, T)

            # nan-safe
            mean_curve = np.nanmean(data, axis=0)
            std_curve  = np.nanstd(data, axis=0, ddof=1) if data.shape[0] > 1 else np.nanstd(data, axis=0)

            mean_ds = mean_curve[::downsample_step]
            std_ds  = std_curve[::downsample_step]

            ax.plot(time_ds, mean_ds, label=f"{cond} (n={data.shape[0]})", color=color_map[cond])
            ax.fill_between(
                time_ds,
                mean_ds - std_ds,
                mean_ds + std_ds,
                color=color_map[cond],
                alpha=0.2,
                edgecolor="none",
                linewidth=0,
                zorder=2
            )

        ax.set_title(title)
        ax.set_xlabel('Time [ms]')
        ax.set_ylabel(ylabel)
        ax.set_ylim(ylims)  #  y ， tight /

        add_trial_lines(ax, trial_len, pre_stimulus_period, stimulus_period, delay_period, dt)

        #  ylim
        fill_stimulus(ax, trial_len, pre_stimulus_period, stimulus_period, dt)

        if xlim is not None:
            ax.set_xlim(*xlim)

        ax.legend()

    # --- full ---
    fig, ax = plt.subplots(figsize=(12, 5))
    _draw(ax, xlim=(time[0], time[-1]))
    full_path = f'{filename_base}_full.pdf'
    fig.savefig(full_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Saved: {full_path}")

    # --- half (second trial window) ---
    fig, ax = plt.subplots(figsize=(12, 5))
    _draw(ax, xlim=(time[trial_len], time[-1]))
    half_path = f'{filename_base}_half.pdf'
    fig.savefig(half_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Saved: {half_path}")


# ----------------- your loader, updated to collect trial-level traces -----------------
def load_plot_rates(group_name, tool, batch_num=40, results_dir="../results", prefix=""):
    dt = float(bm.get_dt())
    trial_len = int((tool.pre_stimulus_period + tool.stimulus_period + tool.delay_period) / dt)  # points in one trial
    T = trial_len * 2  # two consecutive trials concatenated

    unique_conds = ["RR", "LR", "RL", "LL"]
    # you can change colors if you want; keeping RdBu-like mapping
    color_map = {
        'RR': '#56ac55',
        'LR': '#87ad7e',
        'RL': '#dd9eb3',
        'LL': '#c169a5'
    }

    # collect trial-level trajectories for std shading
    traj = {c: [] for c in unique_conds}

    for b_i in range(batch_num):
        data = np.load(f'{results_dir}/rates_260126_batch{b_i}.npz', allow_pickle=True)
        average_rates = data['average_rates'].item()
        seqs = data['sequences']

        rate_A = average_rates[group_name[0]]
        rate_B = average_rates[group_name[1]]
        rate_diff = rate_A - rate_B

        # each i contributes a 2-trial window: [i, i+1]
        for i in range(1, len(seqs) - 1):
            key = str(seqs[i]) + str(seqs[i + 1])  # ensure string
            if key not in traj:
                continue
            seg = rate_diff[i * trial_len:(i + 2) * trial_len]
            if seg.shape[0] != T:
                continue
            traj[key].append(seg)

    # build time axis (ms)
    time = np.arange(T) * dt

    # y-lims based on all collected traces
    ylims = compute_y_range_nan_safe(traj, unique_conds)

    title = f"{group_name[0]} - {group_name[1]}"
    ylabel = "Firing rate diff (Hz)"  # rate_A - rate_B (keep consistent)
    filename_base = f"../figures/{date_str}_{prefix}{group_name[0]}_{group_name[1]}"

    plot_avg_trajectory_with_std(
        time=time,
        trial_len=trial_len,
        data_dict=traj,
        unique_conds=unique_conds,
        color_map=color_map,
        ylims=ylims,
        title=title,
        ylabel=ylabel,
        filename_base=filename_base,
        pre_stimulus_period=tool.pre_stimulus_period,
        stimulus_period=tool.stimulus_period,
        delay_period=tool.delay_period,
        dt=dt,
        downsample_step=50
    )
