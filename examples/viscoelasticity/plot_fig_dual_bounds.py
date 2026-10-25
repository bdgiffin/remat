from math import pi, sin
import sys

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

SCENARIOS = [
    # {
    #     "relaxation_time": 0.3,
    #     "dt": 1.0e-1,
    #     "Nsteps": 20,
    #     "Nsub_steps": 1,
    #     "epsilon0": 0.1,
    #     "bc_name": "right_node_step",
    #     "overflow_limit": 1e6,
    # },
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 3000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_step",
        "overflow_limit": 1e6,
    },
    # {
    #     "relaxation_time": 0.3,
    #     "dt": 1.0e-4,
    #     "Nsteps": 25000,
    #     "Nsub_steps": 1,
    #     "epsilon0": 0.1,
    #     "bc_name": "right_node_step",
    #     "overflow_limit": 1e6,
    # },
    #     {
    #     "relaxation_time": 0.3,
    #     "dt": 1.0e-5,
    #     "Nsteps": 210000,
    #     "Nsub_steps": 1,
    #     "epsilon0": 0.1,
    #     "bc_name": "right_node_step",
    #     "overflow_limit": 1e6,
    # },

    # {
    #     "relaxation_time": 0.3,
    #     "dt": 1.0e-6,
    #     "Nsteps": 360000,
    #     "Nsub_steps": 1,
    #     "epsilon0": 0.1,
    #     "bc_name": "right_node_step",
    #     "overflow_limit": 1e6,
    # },
]


def format_tau(value,dt):
    # shows the relaxation time in math mode and dt size in the legend
    # show dt in 10^-n format if possible
    if dt == 1e-3:
        return rf"$\tau = {value:.2f}, \Delta t =  10^{{-3}}$ "
    elif dt == 1e-4:
        return rf"$\tau = {value:.2f}, \Delta t =  10^{{-4}}$ "
    elif dt == 1e-5:
        return rf"$\tau = {value:.2f}, \Delta t =  10^{{-5}}$ "
    else: 
        return rf"$\tau = {value:.2f}, \Delta t = {dt}$ "


H = 1.0e-6          # lattice spacing of the decimal fixed-point lattice (Fixed<10,-6>)
INT32_MAX = 2**31 - 1
M_MAX = 2**30       # largest admissible ancillary mantissa, half of the int32 range
T_ANCHOR = 1.0      # time at which the homogeneous-growth reference passes through the simulation


def rational_factor(dt, relaxation_time):
    """Rational relaxation factor p/q of eq:rational-approximation, as computed in REMAT."""
    q = int(np.sqrt(INT32_MAX))
    p = int(np.exp(-dt / relaxation_time) * q)
    return p, q


def compute_predicted_bounds(times, relaxation_time, dt, y_anchor):
    # Homogeneous growth at the rate 1/tau (eq:dual-lower), through the simulated value at T_ANCHOR
    lower = y_anchor * np.exp((times - T_ANCHOR) / relaxation_time)

    # Worst-case bound eq:dual-upper with m_0 = 0, converted from mantissas to values
    p, q = rational_factor(dt, relaxation_time)
    ratio = q / p
    n = np.rint(times / dt)
    upper = H * (q - 1) * (ratio**n - 1.0) / (ratio - 1.0)

    # Flush interval eq:flush-time for the sizing level M_MAX
    bound = np.log(1.0 + (ratio - 1.0) * M_MAX / (q - 1)) / np.log(ratio)
    t_max = (np.ceil(bound) - 1) * dt
    return lower, upper, t_max


def main():
    fig, ax = plt.subplots(figsize=(6, 3.5))  # same size as plot_fig_dual_vs_tau.py
    params = SCENARIOS[0]
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
    history = np.abs(result["state_history"][STATE_TO_PLOT])
    times = np.arange(1, history.size + 1) * params["dt"]
    y_anchor = history[np.argmin(np.abs(times - T_ANCHOR))]
    lower, upper, t_max = compute_predicted_bounds(times, params["relaxation_time"], params["dt"], y_anchor)

    ax.plot(times, history, linewidth=1.6, color="#2b738eff",
            label=format_tau(params["relaxation_time"], params["dt"]))
    ax.plot(times, lower, linestyle="--", linewidth=1.2, color="#f9826bff",
            label=r"homogeneous growth $\propto e^{t/\tau}$")
    ax.plot(times, upper, linestyle="--", linewidth=1.2, color="k",
            label="worst-case bound")
    ax.axhline(H * M_MAX, color="gray", linestyle=":", linewidth=1.0)
    ax.text(0.02 * times[-1], H * M_MAX * 1.25, r"sizing level $h\,m_{\max}$", fontsize="small", va="bottom")
    ax.axvline(t_max, color="gray", linestyle=":", linewidth=1.0)
    ax.text(t_max - 0.03, 2.0e-2, r"$N_{\max}\Delta t$", fontsize="medium", va="bottom", ha="right")

    ax.set_xlabel("time (s)", fontsize="large")
    ax.set_ylabel(r"ancillary magnitude $|\varepsilon^{v*}|$", fontsize="large")
    ax.legend(loc="lower right", fontsize="medium", frameon=False)
    ax.set_yscale("log")
    ax.set_xlim(0.0, params["Nsteps"] * params["dt"])
    ax.set_ylim(1.0e-2, 1.0e5)

    fig.tight_layout()
    fig.savefig("dual_bounds.pdf", dpi=200)
    fig.savefig("dual_bounds.png", dpi=200)


if __name__ == "__main__":
    main()
