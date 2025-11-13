from math import pi, sin
import sys
import matplotlib.pyplot as plt
import numpy as np

sys.path.append("../../install/package/")

import REMAT


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
    include_backward=False,
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


STATE_X_TO_PLOT = "axial_strain"
STATE_Y_TO_PLOT = "axial_stress"
set_integrator_type = b"float_truss_visco"
overflow_limit_value = 1e6

SCENARIOS = [
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 15000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_step",
        "overflow_limit": overflow_limit_value,
    },
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 15000,
        "Nsub_steps": 1,
        "epsilon0": 0.02,
        "bc_name": "right_node_constant_rate",
        "overflow_limit": overflow_limit_value,
    },
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 15000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_sinusoidal",
        "overflow_limit": overflow_limit_value,
    },
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 15000,
        "Nsub_steps": 1,
        "epsilon0": 0.2,
        "bc_name": "right_node_clipped_sinusoid",
        "overflow_limit": overflow_limit_value,
    },
]


def format_tau(value):
    return rf"$\tau = {value:.2f}$"


def strain_limits(params):
    if params["bc_name"] == "right_node_constant_rate":
        return (
            0.0,
            params["epsilon0"] * params["Nsteps"] * params["dt"] * 1.1,
        )
    if params["bc_name"] == "right_node_step":
        return (0.0, 0.12)
    return (-0.12, 0.12)


def stress_limits(params):
    if params["bc_name"] == "right_node_step":
        return (-0.02, 0.12)
    if params["bc_name"] == "right_node_constant_rate":
        return (0.0, 0.01)
    if params["bc_name"] == "right_node_sinusoidal":
        return (-0.09, 0.09)
    return (-0.06, 0.06)


def main():

    for params in SCENARIOS:
        fig, ax = plt.subplots(figsize=(5.0, 5.0))
        result = run_truss_relaxation(
            dt=params["dt"],
            Nsteps=params["Nsteps"],
            Nsub_steps=params["Nsub_steps"],
            epsilon0=params["epsilon0"],
            bc_name=params["bc_name"],
            relaxation_time=params["relaxation_time"],
            include_backward=True,
            record_states=(STATE_X_TO_PLOT, STATE_Y_TO_PLOT),
            set_integrator_type=set_integrator_type,
            overflow_limit=params["overflow_limit"],
        )

        histories = result["state_history"]
        strain_history = histories[STATE_X_TO_PLOT]["forward"]
        stress_history = histories[STATE_Y_TO_PLOT]["forward"]
        backward_strain = histories[STATE_X_TO_PLOT].get("backward")
        backward_stress = histories[STATE_Y_TO_PLOT].get("backward")

        ax.plot(
            strain_history,
            stress_history,
            linewidth=1.6,
            marker=None,
            label="forward",
            color="#2b738eff",
        )
        ax.plot(strain_history[0], stress_history[0], '*', color="#2b738eff")

        if backward_strain is not None and backward_stress is not None:
            ax.plot(
                backward_strain[::-1][1:],
                backward_stress[::-1][1:],
                linewidth=1.6,
                marker=None,
                dashes=(6, 8),
                label="backward",
                color="#f9826bff",
            )
            ax.plot(backward_strain[::-1][-1], backward_stress[::-1][-1], '*', color="#f9826bff")

        ax.set_xlabel("axial strain", fontsize="large")
        ax.set_ylabel(r"axial stress ($\sigma_{xx}$)", fontsize="large")
        ax.set_title(format_tau(params["relaxation_time"]))
        legend = ax.legend(loc="upper left", fontsize="medium")
        legend.get_texts()[0].set_color("#2b738eff")
        if backward_strain is not None:
            legend.get_texts()[1].set_color("#f9826bff")
            legend.get_lines()[1].set_linestyle("--")
        ax.set_xlim(*strain_limits(params))
        ax.set_ylim(*stress_limits(params))
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        fig.tight_layout()
        suffix = params["bc_name"].replace("right_node", "")
        fig.savefig(f"stress_vs_strain_{suffix}.svg", dpi=200)
        fig.clf()


if __name__ == "__main__":
    main()
