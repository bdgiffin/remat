from contextlib import contextmanager
from math import pi, sin
import os
import sys
from matplotlib import colors as mcolors
from matplotlib.cm import ScalarMappable
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np

sys.path.append("../../install/package/")

import REMAT

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)


@contextmanager
def suppress_c_stdout(enabled=True):
    if not enabled:
        yield
        return

    sys.stdout.flush()
    saved_stdout_fd = os.dup(1)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 1)
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved_stdout_fd, 1)
        os.close(saved_stdout_fd)
        os.close(devnull_fd)


def set_material_parameters(relaxation_time, overflow_limit):
    parameter_values = {
        b"body_force_x": 0.0,
        b"body_force_y": 0.0,
        b"mass_damping_factor": 0.0,
        b"dt_scale_factor": 1.0,
        b"density": 1.0,
        b"youngs_modulus": 1.0,
        b"poissons_ratio": 0.25,
        b"truss_density": 1.0,
        b"truss_youngs_modulus": 1.0,
        b"area": 1.0,
        b"relaxation_time": relaxation_time,
        b"mat_overflow_limit": overflow_limit,
    }
    for name, value in parameter_values.items():
        REMAT.API.define_parameter(name, value)


def create_geometry():
    gauge_length = 1
    coordinates = np.array([[0.0, 0.0], [gauge_length, 0.0]], dtype=np.double)
    velocities = np.array([[0.0, 0.0], [0.0, 0.0]], dtype=np.double)
    fixity = np.zeros_like(coordinates, dtype=np.bool_)
    fixity[0, :] = True
    fixity[1, 1] = True

    connectivity = np.zeros((0, 4), dtype=np.int32)
    contacts = []
    truss_connectivity = np.array([[0, 1]], dtype=np.int32)

    REMAT.create_geometry(
        coordinates,
        velocities,
        fixity,
        connectivity,
        contacts,
        truss_connectivity,
    )
    return coordinates


def build_bc_function(left_x, epsilon0, name):

    def right_node_step(time, x, _y):
        return 0.0 if time == 0.0 else epsilon0 * (x - left_x)

    def right_node_constant_rate(time, x, _y):
        return (epsilon0 * time) * (x - left_x)

    def right_node_sinusoidal(time, x, _y):
        omega = 2.0 * pi  # 1 Hz by default
        eps_t = epsilon0 * sin(omega * time)
        return eps_t * (x - left_x)

    def right_node_clipped_sinusoid(time, x, _y):
        omega = 2.0 * pi
        amin = -0.5 * epsilon0
        amax = 0.5 * epsilon0
        clipped = max(min(epsilon0 * sin(omega * time), amax), amin)
        return clipped * (x - left_x)

    bc_map = {
        "right_node_step": right_node_step,
        "right_node_constant_rate": right_node_constant_rate,
        "right_node_sinusoidal": right_node_sinusoidal,
        "right_node_clipped_sinusoid": right_node_clipped_sinusoid,
    }
    if name not in bc_map:
        raise ValueError(f"Unknown boundary condition '{name}'")
    return bc_map[name]


def validate_required_params(**params):
    for name, value in params.items():
        if value is None:
            raise ValueError(f"'{name}' must be provided.")


def run_truss_relaxation(
    *,
    dt,
    Nsteps,
    Nsub_steps,
    epsilon0,
    bc_name,
    relaxation_time,
    overflow_limit,
    record_states=None,
    include_backward=True,
):

    validate_required_params(
        dt=dt,
        Nsteps=Nsteps,
        Nsub_steps=Nsub_steps,
        epsilon0=epsilon0,
        relaxation_time=relaxation_time,
        overflow_limit=overflow_limit,
    )

    REMAT.API.set_integrator_type(b"fixed_truss_visco")
    set_material_parameters(relaxation_time, overflow_limit)
    coordinates = create_geometry()

    left_x = coordinates[0, 0]
    bc_function = build_bc_function(left_x, epsilon0, bc_name)
    REMAT.define_displacement_bc(np.array([1], dtype=np.int32), 0, bc_function)

    REMAT.API.initialize()

    state_names = list(record_states) if record_states else []
    recorded_state_values = {}
    for name in state_names:
        recorded_state_values[name] = []
    forward_times = []

    def capture_state():
        for name in state_names:
            field = REMAT.get_field(b"truss", name)
            if field is None:
                raise RuntimeError(f"Truss field '{name}' is not available.")
            recorded_state_values[name].append(float(field[0]))

    for _ in range(1, Nsteps + 1):
        current_time = REMAT.API.update_state(dt,Nsub_steps,REMAT.PASS_FORWARD)
        forward_times.append(current_time)
        if state_names:
            capture_state()

    if include_backward:
        for _ in range(Nsteps, 0, -1):
            REMAT.API.update_state(dt,Nsub_steps,REMAT.PASS_BACKWARD)

    result = {"forward_time": np.asarray(forward_times)}
    if state_names:
        state_history = {}
        for name in state_names:
            state_history[name] = np.asarray(recorded_state_values[name])
        result["state_history"] = state_history
    return result


STATE_TO_PLOT = "dual_viscous_strain"
OUTPUT_FILE = "plot_fig_dual_vs_tau.pdf"
MIN_POSITIVE_MAGNITUDE = 1.0e-8
DISPLAY_YMIN = 1.0e-2
DISPLAY_YMAX = 1.0e3
BUFFER_LEVEL = 5e2
DEFAULT_RATE_MARKER_OFFSET = 0.55
DEFAULT_RATE_MARKER_TIME_FACTOR = 7.0
TAU_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "tau_gray",
    ["0.12", "0.68"],
)

SCENARIOS = [
    {
        "relaxation_time": 0.03,
        "dt": 1.0e-3,
        "Nsteps": 1000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_step",
        "overflow_limit": 1e6,
        "rate_marker_offset": 0.45,
        "rate_marker_time": 0.19,
    },
    {
        "relaxation_time": 0.1,
        "dt": 1.0e-3,
        "Nsteps": 1100,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_step",
        "overflow_limit": 1e6,
        "rate_marker_offset": 0.70,
        "rate_marker_time": 0.80,
    },
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 3000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_step",
        "overflow_limit": 1e6,
        "rate_marker_offset": 0.85,
        "rate_marker_time": 1.6,
    },
    {
        "relaxation_time": 1,
        "dt": 1.0e-3,
        "Nsteps": 5000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_step",
        "overflow_limit": 1e6,
        "rate_marker_offset": 0.75,
        "rate_marker_time": 2.60,
    },
]


TAU_NORM = mcolors.LogNorm(
    vmin=min(params["relaxation_time"] for params in SCENARIOS),
    vmax=max(params["relaxation_time"] for params in SCENARIOS),
)


def information_rate(relaxation_time):
    return 1.0 / (relaxation_time * np.log(2.0))


def format_tau(value):
    return rf"$\tau={value:g}$"


def local_semilog_rate(times, values, index, half_window):
    lower = max(0, index - half_window)
    upper = min(times.size, index + half_window + 1)
    if upper - lower < 3:
        return np.nan

    log_value = np.log(values[lower:upper])
    if not np.all(np.isfinite(log_value)):
        return np.nan

    rate, _intercept = np.polyfit(times[lower:upper], log_value, 1)
    return rate


def choose_rate_marker(
    times,
    values,
    relaxation_time,
    marker_offset,
    marker_time,
    x_limit,
    y_limit,
):
    half_window = max(5, min(60, times.size // 25))
    if times.size <= 2 * half_window + 1:
        return None

    if marker_time is None:
        marker_time = DEFAULT_RATE_MARKER_TIME_FACTOR * relaxation_time
    target_time = np.clip(
        marker_time,
        times[half_window],
        times[-half_window - 1],
    )
    candidate_indices = np.arange(half_window, times.size - half_window)
    marker_options = []
    duration = max(0.015, min(0.22 * relaxation_time, 0.18))

    for index in candidate_indices:
        rate = local_semilog_rate(times, values, int(index), half_window)
        if not np.isfinite(rate) or not (0.2 <= rate <= 250.0):
            continue
        if times[index] + duration >= x_limit:
            continue
        marker_base = values[index] * marker_offset
        if marker_base * np.exp(rate * duration) >= 0.80 * y_limit:
            continue
        if marker_base <= 1.2 * DISPLAY_YMIN:
            continue
        score = abs(np.log(times[index] / target_time))
        marker_options.append((score, int(index), float(rate), duration))

    if not marker_options:
        return None
    return min(marker_options, key=lambda item: item[0])[1:]


def format_rate(value):
    if value >= 10.0:
        return f"{value:.0f}"
    return f"{value:.1f}"


def add_rate_marker(ax, x, y, rate, color, marker_offset, x_limit, duration):
    y = y * marker_offset
    x_right = x + duration
    y_top = y * np.exp(rate * duration)
    if x_right * 1.05 > 0.92 * x_limit:
        text_x = x - 0.03
        ha = "right"
    else:
        text_x = x_right + 0.03
        ha = "left"

    ax.plot([x, x_right], [y, y], color=color, linewidth=0.8, alpha=0.95)
    ax.plot([x_right, x_right], [y, y_top], color=color, linewidth=0.8, alpha=0.95)
    ax.plot([x, x_right], [y, y_top], color=color, linewidth=0.8, alpha=0.95)
    ax.text(
        text_x,
        np.sqrt(y * y_top),
        rf"$r\approx{format_rate(rate)}\,\mathrm{{s}}^{{-1}}$",
        color=color,
        fontsize=8,
        ha=ha,
        va="center",
    )


def line_angle_degrees(ax, x, y, rate, duration):
    x2 = x + duration
    y2 = y * np.exp(rate * duration)
    p1 = ax.transData.transform((x, y))
    p2 = ax.transData.transform((x2, y2))
    return np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))


def add_tau_label(ax, x, y, rate, color, relaxation_time, duration):
    label_y = min(y * 1.35, DISPLAY_YMAX / 1.25)
    rotation = line_angle_degrees(ax, x, label_y, rate, duration)
    ax.text(
        x,
        label_y,
        format_tau(relaxation_time),
        color=color,
        fontsize=9,
        rotation=rotation,
        rotation_mode="anchor",
        ha="center",
        va="bottom",
    )


def main():
    fig, ax = plt.subplots(figsize=(6, 3.5))
    histories = []

    for params in SCENARIOS:
        with suppress_c_stdout():
            result = run_truss_relaxation(
                dt=params["dt"],
                Nsteps=params["Nsteps"],
                Nsub_steps=params["Nsub_steps"],
                epsilon0=params["epsilon0"],
                bc_name=params["bc_name"],
                relaxation_time=params["relaxation_time"],
                include_backward=False,
                record_states=(STATE_TO_PLOT,),
                overflow_limit=params["overflow_limit"],
            )

        history = result["state_history"][STATE_TO_PLOT]
        times = np.arange(1, history.size + 1) * params["dt"]
        magnitude = np.maximum(np.abs(history), MIN_POSITIVE_MAGNITUDE)
        color = TAU_CMAP(TAU_NORM(params["relaxation_time"]))

        ax.plot(
            times,
            magnitude,
            linewidth=1.6,
            marker=None,
            color=color,
            label="_nolegend_",
        )

        histories.append((params, times, magnitude, color))

    axis_dt = SCENARIOS[0]["dt"]
    x_max = 3400 * axis_dt
    ax.set_yscale("log")
    ax.set_xlim(-200 * axis_dt, x_max)
    ax.set_ylim(DISPLAY_YMIN, DISPLAY_YMAX)

    for _params, times, magnitude, color in histories:
        marker_offset = _params.get("rate_marker_offset", DEFAULT_RATE_MARKER_OFFSET)
        marker_time = _params.get("rate_marker_time")
        marker = choose_rate_marker(
            times,
            magnitude,
            _params["relaxation_time"],
            marker_offset,
            marker_time,
            x_max,
            DISPLAY_YMAX,
        )
        if marker is None:
            continue
        marker_index, rate, duration = marker
        add_rate_marker(
            ax,
            times[marker_index],
            magnitude[marker_index],
            rate,
            color,
            marker_offset,
            x_max,
            duration=duration,
        )
        add_tau_label(
            ax,
            times[marker_index],
            magnitude[marker_index],
            rate,
            color,
            _params["relaxation_time"],
            duration,
        )

    ax.axhline(
        BUFFER_LEVEL,
        color="0.35",
        linestyle=":",
        linewidth=1.0,
        label=rf"example buffer level $y_{{\max}}={BUFFER_LEVEL:g}$",
    )

    ax.set_xlabel("time (s)", fontsize="large")
    ax.set_ylabel(r"ancillary magnitude $|\varepsilon^{v*}|$", fontsize="large")
    ax.minorticks_on()
    ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0))
    ax.yaxis.set_minor_locator(
        mticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1)
    )
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.tick_params(axis="y", which="major", length=5, width=0.8)
    ax.tick_params(axis="y", which="minor", length=3, width=0.6)

    sm = ScalarMappable(norm=TAU_NORM, cmap=TAU_CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(r"relaxation time $\tau$ (s)", fontsize="medium")

    # ax.legend(loc="lower right", fontsize="small")

    fig.tight_layout()
    fig.savefig(OUTPUT_FILE, format="pdf", dpi=200)
    fig.clf()


if __name__ == "__main__":
    main()
