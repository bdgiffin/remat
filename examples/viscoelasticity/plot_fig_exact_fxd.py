from math import pi, sin
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
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
        omega = 0.2 * pi
        eps_t = epsilon0 * sin(omega * time)
        return eps_t * (x - left_x)

    def right_node_clipped_sinusoid(time, x, _y):
        omega = 0.2 * pi
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


def compute_exact_stress(times, params):
    tau = params["relaxation_time"]
    eps0 = params["epsilon0"]
    if params["bc_name"] == "right_node_constant_rate":
        return eps0 * tau * (1.0 - np.exp(-times / tau))
    raise ValueError("Exact solution is only implemented for constant strain rate.")


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
        current_time = REMAT.API.update_state(+dt, Nsub_steps)
        forward_times.append(current_time)
        if state_names:
            capture_state()

    result = {"forward_time": np.asarray(forward_times)}
    if state_names:
        state_history = {}
        for name in state_names:
            state_history[name] = np.asarray(recorded_state_values[name])
        result["state_history"] = state_history
    return result


def compute_relative_metrics(reference, candidate):
    diff = candidate - reference
    denom = max(np.linalg.norm(reference), np.finfo(np.float64).eps)
    return {
        "max_abs": float(np.max(np.abs(diff))),
        "rel_l2": float(np.linalg.norm(diff) / denom),
    }


STATE_TO_PLOT = "axial_stress"
MODE_COLORS = {
    "exact": "#222222",
    "fixed": "#f9826bff",
}

SCENARIOS = [
    {
        "description": r"Constant strain rate, $\tau=1$, $\Delta t=10^{-2}$",
        "relaxation_time": 0.1,
        "dt": 1.0e-2,
        "Nsteps": 1000,
        "Nsub_steps": 1,
        "epsilon0": 0.001,
        "bc_name": "right_node_constant_rate",
        "overflow_limit": 1.0e0,
    },
]


def plot_exact_vs_fixed(params, fixed_result):
    fig, (ax_state, ax_error) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(5.0, 5.0),
        gridspec_kw={"height_ratios": [3, 1]},
    )

    times = fixed_result["forward_time"]
    fixed_history = fixed_result["state_history"][STATE_TO_PLOT]
    exact_history = compute_exact_stress(times, params)
    metrics = compute_relative_metrics(exact_history, fixed_history)

    ax_state.plot(
        times,
        exact_history,
        linewidth=1.6,
        label="exact",
        color=MODE_COLORS["exact"],
    )
    ax_state.plot(
        times,
        fixed_history,
        linewidth=1.6,
        label="fixed",
        color=MODE_COLORS["fixed"],
    )

    diff = fixed_history - exact_history
    ax_error.plot(
        times,
        diff,
        linewidth=1.2,
        color="#4a4a4aff",
        label=r"$\sigma_{xx}^\mathrm{fixed}-\sigma_{xx}^\mathrm{exact}$",
    )

    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    ax_error.yaxis.set_major_formatter(fmt)

    ax_error.axhline(0.0, color="#00000080", linewidth=1.0, linestyle=":")

    ax_state.set_ylabel(r"axial stress ($\sigma_{xx}$)", fontsize="large")
    ax_state.set_xlim(0, params["dt"] * params["Nsteps"])
    ax_state.legend(loc="best", fontsize="medium")

    ax_error.set_xlabel("time (s)", fontsize="large")
    ax_error.set_ylabel(r"$\sigma_{xx}^\mathrm{fixed}-\sigma_{xx}^\mathrm{exact}$", fontsize="large")
    ax_error.legend(loc="best", fontsize="medium")

    title = (
        f"{params['description']}\n"
        f"max|Δ|={metrics['max_abs']:.2e}, "
        f"rel‖Δ‖₂={metrics['rel_l2']:.2e}"
    )
    ax_state.set_title(title, fontsize="medium")

    fig.tight_layout()
    return fig


def summarize_diagnostics(params, fixed_result):
    times = fixed_result["forward_time"]
    fixed_history = fixed_result["state_history"][STATE_TO_PLOT]
    exact_history = compute_exact_stress(times, params)
    metrics = compute_relative_metrics(exact_history, fixed_history)

    print("-" * 80)
    print(f"Scenario: {params['description']}")
    print(f"  duration = {params['dt'] * params['Nsteps']:.3f}s, dt = {params['dt']}, tau = {params['relaxation_time']}")
    print(f"  forward max|fixed - exact| = {metrics['max_abs']:.6e}")
    print(f"  forward relative L2 error   = {metrics['rel_l2']:.6e}")
    print("-" * 80)


def main():
    for params in SCENARIOS:
        fixed_result = run_truss_relaxation(
            dt=params["dt"],
            Nsteps=params["Nsteps"],
            Nsub_steps=params["Nsub_steps"],
            epsilon0=params["epsilon0"],
            bc_name=params["bc_name"],
            relaxation_time=params["relaxation_time"],
            record_states=(STATE_TO_PLOT,),
            overflow_limit=params["overflow_limit"],
        )
        summarize_diagnostics(params, fixed_result)
        fig = plot_exact_vs_fixed(params, fixed_result)
        suffix = params["bc_name"].replace("right_node", "")
        fig.savefig(
            f"exact_vs_fxd_{suffix}.svg",
            metadata={"Title": str(params["description"])},
            dpi=200,
        )
        fig.clf()


if __name__ == "__main__":
    main()
