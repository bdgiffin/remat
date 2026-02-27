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
    set_integrator_type=None,
    include_backward=True,
):
    validate_required_params(
        dt=dt,
        Nsteps=Nsteps,
        Nsub_steps=Nsub_steps,
        epsilon0=epsilon0,
        relaxation_time=relaxation_time,
        overflow_limit=overflow_limit,
        set_integrator_type=set_integrator_type,
    )

    REMAT.API.set_integrator_type(set_integrator_type)
    set_material_parameters(relaxation_time, overflow_limit)
    coordinates = create_geometry()

    left_x = coordinates[0, 0]
    bc_function = build_bc_function(left_x, epsilon0, bc_name)
    REMAT.define_displacement_bc(np.array([1], dtype=np.int32), 0, bc_function)

    REMAT.API.initialize()

    state_names = list(record_states) if record_states else []
    forward_state_values = {}
    backward_state_values = {}

    for name in state_names:
        forward_state_values[name] = []
        backward_state_values[name] = []

    forward_times = []
    backward_times = []

    def capture_state(storage):
        for name in state_names:
            field = REMAT.get_field(b"truss", name)
            if field is None:
                raise RuntimeError(f"Truss field '{name}' is not available.")
            storage[name].append(float(field[0]))

    for _ in range(1, Nsteps + 1):
        current_time = REMAT.API.update_state(+dt, Nsub_steps)
        forward_times.append(current_time)
        if state_names:
            capture_state(forward_state_values)

    if include_backward:
        for _ in range(Nsteps, 0, -1):
            current_time = REMAT.API.update_state(-dt, Nsub_steps)
            backward_times.append(current_time)
            if state_names:
                capture_state(backward_state_values)

    result = {"forward_time": np.asarray(forward_times)}
    if include_backward:
        result["backward_time"] = np.asarray(backward_times)
    if state_names:
        state_history = {}
        for name in state_names:
            state_history[name] = {"forward": np.asarray(forward_state_values[name])}
            if include_backward:
                state_history[name]["backward"] = np.asarray(
                    backward_state_values[name]
                )
        result["state_history"] = state_history
    return result


STATE_TO_PLOT = "lambda_adjoint"
PRECISION_MODES = [
    ("adjoint float / truss float", b"float_truss_visco"),
    ("adjoint float / truss fixed", b"fixed_truss_visco_adj_float"),
    ("adjoint fixed / truss fixed", b"fixed_truss_visco"),
]
MODE_COLORS = {
    "adjoint float / truss float": "#2b738e",
    "adjoint float / truss fixed": "#f9826b",
    "adjoint fixed / truss fixed": "#6f6f6f",
}
MODE_OVERFLOW_LIMITS = {
    "adjoint float / truss float": 1.0e12,
    "adjoint float / truss fixed": 1.0e12,
    "adjoint fixed / truss fixed": 1.0e12,
}

SCENARIOS = [
    {
        "description": r"Backward adjoint, step strain, $\tau=0.1$, $\Delta t=10^{-3}$",
        "relaxation_time": 0.1,
        "dt": 1.0e-3,
        "Nsteps": 450,
        "Nsub_steps": 1,
        "epsilon0": 2.0,
        "bc_name": "right_node_step",
    },
]


def plot_precision_histories(params, precision_results):
    fig, ax_state = plt.subplots(1, 1, figsize=(6.2, 4.0))

    for mode_name, _ in PRECISION_MODES:
        histories = precision_results[mode_name]["state_history"][STATE_TO_PLOT]
        backward_time = precision_results[mode_name]["backward_time"]
        order = np.argsort(backward_time)
        ax_state.plot(
            backward_time[order],
            histories["backward"][order],
            linewidth=1.6,
            label=f"{mode_name} (backward)",
            color=MODE_COLORS[mode_name],
        )

    ax_state.set_ylabel(r"adjoint state ($\lambda$)", fontsize="large")
    ax_state.set_xlabel("time (s)", fontsize="large")
    ax_state.set_xlim(0.0, params["dt"] * params["Nsteps"])
    ax_state.legend(loc="best", fontsize="small")
    fig.tight_layout()
    return fig


def run_precision_suite(params):
    results = {}
    for mode_name, integrator in PRECISION_MODES:
        results[mode_name] = run_truss_relaxation(
            dt=params["dt"],
            Nsteps=params["Nsteps"],
            Nsub_steps=params["Nsub_steps"],
            epsilon0=params["epsilon0"],
            bc_name=params["bc_name"],
            relaxation_time=params["relaxation_time"],
            include_backward=True,
            record_states=(STATE_TO_PLOT,),
            set_integrator_type=integrator,
            overflow_limit=MODE_OVERFLOW_LIMITS[mode_name],
        )
    return results


def summarize_diagnostics(params, precision_results):
    print("-" * 80)
    print(f"Scenario: {params['description']}")
    print(
        f"  duration = {params['dt'] * params['Nsteps']:.3f}s, "
        f"dt = {params['dt']}, tau = {params['relaxation_time']}"
    )
    for mode_name, _ in PRECISION_MODES:
        history = precision_results[mode_name]["state_history"][STATE_TO_PLOT]["backward"]
        print(
            f"  {mode_name}: "
            f"min={np.min(history): .6e}, max={np.max(history): .6e}, final={history[-1]: .6e}"
        )
    print("-" * 80)


def main():
    for params in SCENARIOS:
        precision_results = run_precision_suite(params)
        summarize_diagnostics(params, precision_results)
        fig = plot_precision_histories(params, precision_results)
        fig.savefig(
            "flt_vs_fxd_adjoint_lambda_backward.svg",
            metadata={"Title": str(params["description"])},
            dpi=200,
        )
        fig.clf()


if __name__ == "__main__":
    main()
