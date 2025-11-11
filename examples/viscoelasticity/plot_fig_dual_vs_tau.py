"""Compare dual-variable histories for multiple relaxation times without
touching the original truss_step_relaxation example."""

from __future__ import annotations

from math import pi, sin
from pathlib import Path
import sys
from typing import Iterable, Optional

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# REMAT setup (mirrors the example script but lives entirely in this module)

repo_root = Path(__file__).resolve().parents[2]
package_path = repo_root / "install" / "package"
sys.path.append(str(package_path))

exodus_available = False
if not sys.platform == "emscripten":
    try:
        from ExodusIO import ExodusIO  # type: ignore
    except ModuleNotFoundError:
        exodus_available = False
    else:
        exodus_available = True

import REMAT  # noqa: E402

gauge_length = 1.0
DEFAULT_DT = 1.0e-3
DEFAULT_NSTEPS = 10
DEFAULT_NSUB_STEPS = 100
DEFAULT_EPSILON0 = 1e-1
DEFAULT_BC = "right_node_step"
DEFAULT_RELAXATION_TIME = 1e-1
DEFAULT_OVERFLOW_LIMIT = 30.0


def _set_material_parameters(relaxation_time: float, overflow_limit: float) -> None:
    """Define the parameters used by the REMAT integrator."""
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


def _create_geometry():
    """Create the two-node truss geometry and register it with REMAT."""
    coordinates = np.array([[0.0, 0.0], [gauge_length, 0.0]], dtype=np.double)
    velocities = np.array([[0.0, 0.0], [0.0, 0.0]], dtype=np.double)
    fixity = np.zeros_like(coordinates, dtype=np.bool_)
    fixity[0, :] = True  # left node fully fixed
    fixity[1, 1] = True  # right node vertical DOF fixed

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


def _build_bc_function(left_x: float, epsilon0: float, name: str):
    """Return the displacement boundary condition selected by name."""

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
    try:
        return bc_map[name]
    except KeyError as exc:
        raise ValueError(f"Unknown boundary condition '{name}'") from exc


def run_truss_relaxation(
    dt: float = DEFAULT_DT,
    Nsteps: int = DEFAULT_NSTEPS,
    Nsub_steps: int = DEFAULT_NSUB_STEPS,
    epsilon0: float = DEFAULT_EPSILON0,
    BC_selected: str = DEFAULT_BC,
    relaxation_time: float = DEFAULT_RELAXATION_TIME,
    record_states: Optional[Iterable[str]] = None,
    include_backward: bool = True,
    exodus_path: Optional[str] = None,
    overflow_limit: float = DEFAULT_OVERFLOW_LIMIT,
    
):
    """
    Run the viscoelastic truss problem and optionally return truss state history.
    """

    REMAT.API.set_integrator_type(b"fixed_truss_visco")
    _set_material_parameters(relaxation_time, overflow_limit)
    coordinates = _create_geometry()

    left_x = coordinates[0, 0]
    bc_function = _build_bc_function(left_x, epsilon0, BC_selected)
    REMAT.define_displacement_bc(np.array([1], dtype=np.int32), 0, bc_function)

    REMAT.API.initialize()

    exo = None
    if exodus_path and exodus_available:
        exo = ExodusIO()
        exo.create(exodus_path)
        exo.output_state()

    state_names = tuple(record_states or ())
    recorded_state_values = {name: [] for name in state_names}
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
        if exo:
            exo.output_state()

    if include_backward:
        for _ in range(Nsteps, 0, -1):
            REMAT.API.update_state(-dt, Nsub_steps)
            if exo:
                exo.output_state()

    if exo:
        exo.output_state()
        exo.finalize()

    result = {"forward_time": np.asarray(forward_times)}
    if state_names:
        result["state_history"] = {
            name: np.asarray(values) for name, values in recorded_state_values.items()
        }
    return result


# ---------------------------------------------------------------------------
# Scenario configuration + plotting

STATE_TO_PLOT = "dual_viscous_strain"

SCENARIOS = [
    {
        "relaxation_time": 0.01,
        "dt": 1.0e-3,
        "Nsteps": 130,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "BC_selected": "right_node_step",
        "overflow_limit": 1e6,
    },
    {
        "relaxation_time": 0.1,
        "dt": 1.0e-3,
        "Nsteps": 950,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "BC_selected": "right_node_step",
        "overflow_limit": 1e6,
    },
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 2600,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "BC_selected": "right_node_step",
        "overflow_limit": 1e6,
    },
    {
        "relaxation_time": 0.5,
        "dt": 1.0e-3,
        "Nsteps": 4000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "BC_selected": "right_node_step",
        "overflow_limit": 1e6,
    },
    {
        "relaxation_time": 0.7,
        "dt": 1.0e-3,
        "Nsteps": 5560,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "BC_selected": "right_node_constant_rate",
        "overflow_limit": 1e6,
    },
    {
        "relaxation_time": 1,
        "dt": 1.0e-3,
        "Nsteps": 7500,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "BC_selected": "right_node_constant_rate",
        "overflow_limit": 1e6,
    },

    # {
    #     "relaxation_time": 0.5,
    #     "dt": 1.0e-3,
    #     "Nsteps": 3500,
    #     "Nsub_steps": 1,
    #     "epsilon0": 5.0e-2,
    #     "BC_selected": "right_node_constant_rate",
    # },
    # {
    #     "relaxation_time": 0.8,
    #     "dt": 1.0e-3,
    #     "Nsteps": 6000,
    #     "Nsub_steps": 1,
    #     "epsilon0": 5.0e-2,
    #     "BC_selected": "right_node_constant_rate",
    # },
    # {
    #     "relaxation_time": 0.9,
    #     "dt": 1.0e-3,
    #     "Nsteps": 6000,
    #     "Nsub_steps": 1,
    #     "epsilon0": 5.0e-2,
    #     "BC_selected": "right_node_sinusoidal",
    # },
    # {
    #     "relaxation_time": 1.0,
    #     "dt": 1.0e-3,
    #     "Nsteps": 6000,
    #     "Nsub_steps": 1,
    #     "epsilon0": 5.0e-2,
    #     "BC_selected": "right_node_sinusoidal",
    # },
]


def _format_tau(value: float) -> str:
    return rf"$\tau = {value:.2f}$"


def main():
    fig, ax = plt.subplots(figsize=(7.5, 4.0))

    for params in SCENARIOS:
        result = run_truss_relaxation(
            dt=params["dt"],
            Nsteps=params["Nsteps"],
            Nsub_steps=params["Nsub_steps"],
            epsilon0=params["epsilon0"],
            BC_selected=params["BC_selected"],
            relaxation_time=params["relaxation_time"],
            include_backward=False,
            record_states=(STATE_TO_PLOT,),
            # integrator_type="fixed_truss_visco",
            overflow_limit=params["overflow_limit"],
        )

        history = result["state_history"][STATE_TO_PLOT]
        steps = np.arange(1, history.size + 1)

        ax.plot(
            steps,
            history,
            linewidth=1.6,
            marker=None,
            label=_format_tau(params["relaxation_time"]),
        )

    ax.set_xlabel("time step")
    # ax.set_ylabel(STATE_TO_PLOT.replace("_", " ").title())
    ax.set_ylabel("dual variable")
    # ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend(loc="upper right", fontsize=10)
    ax.set_xlim(-100, 7000)
    ax.set_ylim(-10, 750)

    fig.tight_layout()
    fig.savefig("plot_fig_dual_vs_tau.svg", dpi=200)

    backend = plt.get_backend().lower()
    if "agg" not in backend:
        plt.show()


if __name__ == "__main__":
    main()

