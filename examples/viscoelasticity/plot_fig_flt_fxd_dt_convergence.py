from math import pi, sin
from contextlib import contextmanager
import os
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
        current_time = REMAT.API.update_state(dt,Nsub_steps,REMAT.PASS_FORWARD)
        forward_times.append(current_time)
        if state_names:
            capture_state(forward_state_values)

    if include_backward:
        for _ in range(Nsteps, 0, -1):
            current_time = REMAT.API.update_state(dt,Nsub_steps,REMAT.PASS_BACKWARD)
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


def compute_error_components(reference, candidate):
    diff = candidate - reference
    denom = max(np.linalg.norm(reference), np.finfo(np.float64).eps)
    max_abs = float(np.max(np.abs(diff)))
    return max_abs, float(denom), float(max_abs / denom)


def compute_max_relative_error(reference, candidate):
    return compute_error_components(reference, candidate)[2]


STATE_TO_PLOT = "axial_stress"
PRECISION_MODES = [
    ("float", b"float_truss_visco"),
    ("fixed", b"fixed_truss_visco"),
]
MODE_COLORS = {
    "float": "#2b738eff",
    "fixed": "#f9826bff",
}
MODE_OVERFLOW_LIMITS = {
    "float": 1.0e56,
    "fixed": 1000.0,
}

# Evaluate Delta t from 1e-5 to 1e-1 s with a fixed duration.
DT_VALUES = np.logspace(-5, -1, num=9)
DURATION = 15.0
SCENARIO = {
    "description": r"Constant strain rate, $\tau=0.3$, $\Delta t/\tau$ convergence",
    "relaxation_time": 0.3,
    "Nsub_steps": 1,
    "epsilon0": 0.1,
    "bc_name": "right_node_constant_rate",
}


def compute_analytical_stress(times, params):
    tau = params["relaxation_time"]
    eps0 = params["epsilon0"]
    if params["bc_name"] == "right_node_step":
        # step strain: sigma = E * eps0 * exp(-t/tau), E=1 here
        return eps0 * np.exp(-times / tau)
    elif params["bc_name"] == "right_node_constant_rate":
        # constant strain rate: epsilon = eps0 * t, sigma = E * eps0 * tau * (1 - exp(-t/tau))
        return eps0 * tau * (1.0 - np.exp(-times / tau))
    else:
        raise ValueError("No analytical solution implemented for this bc.")


def run_precision_suite(params, quiet=True):
    results = {}
    for mode_name, integrator in PRECISION_MODES:
        with suppress_c_stdout(quiet):
            results[mode_name] = run_truss_relaxation(
                dt=params["dt"],
                Nsteps=params["Nsteps"],
                Nsub_steps=params["Nsub_steps"],
                epsilon0=params["epsilon0"],
                bc_name=params["bc_name"],
                relaxation_time=params["relaxation_time"],
                include_backward=False,
                record_states=(STATE_TO_PLOT,),
                set_integrator_type=integrator,
                overflow_limit=MODE_OVERFLOW_LIMITS[mode_name],
            )
    return results


def evaluate_time_step_convergence(dt_values, scenario):
    results = []
    for dt in dt_values:
        nsteps = max(1, int(np.round(DURATION / dt)))
        params = {**scenario, "dt": float(dt), "Nsteps": nsteps}
        precision_results = run_precision_suite(params)
        float_hist = precision_results["float"]["state_history"][STATE_TO_PLOT]["forward"]
        fixed_hist = precision_results["fixed"]["state_history"][STATE_TO_PLOT]["forward"]
        analytical_hist = compute_analytical_stress(precision_results["float"]["forward_time"], params)

        float_max_abs, reference_l2_norm, float_err = compute_error_components(
            analytical_hist, float_hist
        )
        fixed_max_abs, _, fixed_err = compute_error_components(
            analytical_hist, fixed_hist
        )

        results.append(
            {
                "dt": float(dt),
                "Nsteps": nsteps,
                "reference_l2_norm": reference_l2_norm,
                "float_max_abs_error": float_max_abs,
                "float_err": float_err,
                "fixed_max_abs_error": fixed_max_abs,
                "fixed_err": fixed_err,
            }
        )

        print("-" * 80)
        print(f"Δt = {dt:.2e} s, Nsteps = {nsteps}")
        print(f"  reference L2 norm = {reference_l2_norm:.6e}")
        print(
            f"  float max abs error = {float_max_abs:.6e}, "
            f"normalized = {float_err:.6e}"
        )
        print(
            f"  fixed max abs error = {fixed_max_abs:.6e}, "
            f"normalized = {fixed_err:.6e}"
        )
    return results


def print_convergence_orders(dt_results):
    if len(dt_results) < 2:
        print("Need at least two dt values to compute convergence orders.")
        return

    data = sorted(dt_results, key=lambda x: x["dt"], reverse=True)  # coarse to fine
    print("Convergence order (max relative error vs analytical):")
    print("  dt_coarse | dt_fine | float_order | fixed_order")
    print("  -----------------------------------------------")
    for i in range(len(data) - 1):
        coarse = data[i]
        fine = data[i + 1]
        f_order = np.nan
        g_order = np.nan
        if coarse["float_err"] > 0.0 and fine["float_err"] > 0.0:
            f_order = np.log(coarse["float_err"] / fine["float_err"]) / np.log(
                coarse["dt"] / fine["dt"]
            )
        if coarse["fixed_err"] > 0.0 and fine["fixed_err"] > 0.0:
            g_order = np.log(coarse["fixed_err"] / fine["fixed_err"]) / np.log(
                coarse["dt"] / fine["dt"]
            )
        print(
            f"  {coarse['dt']:.3e} | {fine['dt']:.3e} | {f_order:11.4f} | {g_order:11.4f}"
        )
    print("  -----------------------------------------------")


def plot_dt_errors(dt_results, scenario):
    fig, ax = plt.subplots(figsize=(4.5, 3.0))

    dts = np.asarray([entry["dt"] for entry in dt_results])
    normalized_steps = dts / scenario["relaxation_time"]
    max_rel_error_float = np.asarray([entry["float_err"] for entry in dt_results])
    max_rel_error_fixed = np.asarray([entry["fixed_err"] for entry in dt_results])

    ax.plot(
        normalized_steps,
        max_rel_error_float,
        marker="o",
        linewidth=1.6,
        label=r"float",
        color=MODE_COLORS["float"],
    )
    ax.plot(
        normalized_steps,
        max_rel_error_fixed,
        marker="s",
        linewidth=1.6,
        label=r"fixed",
        color=MODE_COLORS["fixed"],
    )

    fixed_min_error = np.min(max_rel_error_fixed)
    ax.axhline(
        fixed_min_error,
        color=MODE_COLORS["fixed"],
        linestyle="--",
        linewidth=1.0,
        alpha=0.85,
        label=r"fixed min.",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"normalized time step $\Delta t/\tau$", fontsize="large")
    ax.set_ylabel("normalized max stress error", fontsize="large")
    # ax.set_title(
    #     f"{scenario['description']}, duration={DURATION}s",
    #     fontsize="medium",
    # )
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.2)

    ax.legend(loc="best", fontsize="medium")
    fig.tight_layout()
    return fig


def main():
    dt_results = evaluate_time_step_convergence(DT_VALUES, SCENARIO)
    #print_convergence_orders(dt_results)
    fig = plot_dt_errors(dt_results, SCENARIO)
    fig.savefig(
        "flt_vs_fxd_dt_convergence.pdf",
        metadata={"Title": str(SCENARIO["description"])},
        dpi=200,
    )
    fig.clf()


if __name__ == "__main__":
    main()
