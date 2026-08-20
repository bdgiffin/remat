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

# Assumption: enforce a tiny positive floor to keep log-scale bounds well defined.
MIN_POSITIVE_BOUND = 1.0e-6


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


def estimate_initial_dual(first_value, dt, relaxation_time):
    if first_value <= 0.0:
        # Assumption: a tiny positive value when the dual variable value starts at or less than zero
        return MIN_POSITIVE_BOUND
    # Assumption: approximate y(t=0) from the first stored state
    return first_value / np.exp(dt / relaxation_time)


def compute_predicted_bounds(times, relaxation_time, y0, dt, overflow_limit):
    # The Lyapunov growth rate lambda=1/tau
    lower = np.maximum(y0 * np.exp(times / relaxation_time), MIN_POSITIVE_BOUND)

    # Assumption: model worst-case round-off accumulation with R(dt)=exp(dt/tau)/dt.
    r_max = np.exp(dt / relaxation_time) / dt
    upper = (y0 + r_max * relaxation_time) * np.exp(times / relaxation_time) - r_max * relaxation_time

    if overflow_limit is not None:
        upper = np.minimum(upper, overflow_limit)

    base = y0 + r_max * relaxation_time
    limit_plus_round_off = None if overflow_limit is None else overflow_limit + r_max * relaxation_time
    t_max = None
    if limit_plus_round_off and base > 0:
        t_max = (np.log(limit_plus_round_off) - np.log(base)) * relaxation_time

    return lower, upper, t_max


def main():
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    axis_dt = SCENARIOS[0]["dt"]
    bounds_label_added = False
    y_min_values = []
    y_max_values = []

    for params in SCENARIOS:
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
        steps = np.arange(1, history.size + 1) * params["dt"]

        (line,) = ax.plot(
            steps,
            history,
            linewidth=1.1,
            marker=None,
            label=format_tau(params["relaxation_time"],params["dt"]),
            color='green'
        )

        color = line.get_color()
        
        initial_dual = estimate_initial_dual(history[0], params["dt"], params["relaxation_time"])
        lower, upper, _t_max = compute_predicted_bounds(
            steps,
            params["relaxation_time"],
            initial_dual,
            params["dt"],
            params["overflow_limit"],
        )

        lower_label = "predicted bounds (lower/upper)" if not bounds_label_added else None
        ax.plot(
            steps,
            lower,
            linestyle="--",
            linewidth=1.0,
            color=color,
            alpha=0.9,
            label=lower_label,
        )
        ax.plot(
            steps,
            upper,
            linestyle="--",
            linewidth=1.0,
            color=color,
            alpha=0.9,
            label=None,
        )
        bounds_label_added = True

        y_min_values.append(np.min([history.min(), lower.min()]))
        y_max_values.append(np.max([history.max(), upper.max()]))

    ax.set_xlabel("time (s)", fontsize="large")
    ax.set_ylabel("ancillary variable", fontsize="large")
    # ax.legend(loc="lower right", bbox_to_anchor=(0.0, 1.0), fontsize=10)
    ax.legend(loc="lower right", fontsize='medium')

    # plot y=200 line
    ax.axhline(y=200, color='gray', linestyle=':', linewidth=1.0)



    ax.set_xlim(-200 * axis_dt, 3400 * axis_dt)
    # ax.set_xlim(-.001, max(max(np.arange(1, params["Nsteps"] + 1) * params["dt"]) for params in SCENARIOS) + .001)

    ax.set_yscale("log")

    ymin = max(MIN_POSITIVE_BOUND, min(y_min_values) * 0.8)
    ymax = max(y_max_values) * 1.2
    # ax.set_ylim(ymin, ymax)
    ax.set_ylim(10**-2, 1000)
    ax.set_ylim(10**-6, 10000)
    ax.text(
        0.09,
        -0.065,
        r"$t_{\max}$",
        transform=ax.transAxes,
        fontsize='13',
    )


    # ax.set_ylim(10**-2, 1000)


    fig.tight_layout()
    fig.savefig("dual_bounds.pdf", dpi=200)


if __name__ == "__main__":
    main()
