"""
Minimal driver to compute objective f and its sensitivity df/dtau for the
uniaxial viscoelastic truss used in plot_fig_dual_vs_tau.py.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

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
    gauge_length = 1.0
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


def build_bc_function(left_x, epsilon0):
    def right_node_step(time, x, _y):
        return 0.0 if time == 0.0 else epsilon0 * (x - left_x)

    return right_node_step


def run_forward(dt, Nsteps, epsilon0, relaxation_time, overflow_limit, run_backward=True):
    REMAT.API.set_integrator_type(b"fixed_truss_visco")
    set_material_parameters(relaxation_time, overflow_limit)
    create_geometry()
    left_x = 0.0
    bc_function = build_bc_function(left_x, epsilon0)
    REMAT.define_displacement_bc(np.array([1], dtype=np.int32), 0, bc_function)
    REMAT.API.initialize()

    for _ in range(Nsteps):
        REMAT.API.update_state(dt, 1)
    if run_backward:
        for _ in range(Nsteps):
            REMAT.API.update_state(-dt, 1)

    df_dtau_val = REMAT.get_field(b"truss", "df_dtau")[0]
    return df_dtau_val


def main():
    # Nominal point (matches plot_fig_dual_vs_tau.py)
    tau = 0.01
    dt = 1.0e-3
    Nsteps = 100
    epsilon0 = 0.1
    overflow_limit = 90
    df_dtau = run_forward(dt, Nsteps, epsilon0, tau, overflow_limit, run_backward=True)
    print(f"df/dtau  adjoint method   = {df_dtau:.10f}")

    # Sweep tau logarithmically from 1e-3 to 10
    sweep_taus = np.logspace(-2, 1, num=30)
    df_vals = []
    for tau_i in sweep_taus:
        df_i = run_forward(dt,Nsteps, epsilon0, float(tau_i), overflow_limit, run_backward=True)
        df_vals.append(df_i)

    sweep_taus = np.asarray(sweep_taus)
    df_vals = np.asarray(df_vals)

    fig, ax0 = plt.subplots(1, 1, figsize=(7, 4.0))

    ax0.semilogx(sweep_taus, df_vals, marker="o", linewidth=1.2)
    ax0.set_xlabel(r"$\tau$")
    ax0.set_ylabel(r"$df/d\tau$")
    # ax0.grid(True, which="both", linestyle="--", alpha=0.5)

    fig.tight_layout()
    fig.savefig("sensitivity_vs_tau.svg", dpi=200)
    print("Saved sweep plot to sensitivity_vs_tau.svg")


if __name__ == "__main__":
    main()
