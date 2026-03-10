from math import *
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../install/package/")
import REMAT

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)


def define_parameters(tau, youngs_modulus, mat_overflow_limit):
    REMAT.API.define_parameter(b"body_force_x", 0.0)
    REMAT.API.define_parameter(b"body_force_y", 0.0)
    REMAT.API.define_parameter(b"mass_damping_factor", 0.0)
    REMAT.API.define_parameter(b"dt_scale_factor", 1.0)

    REMAT.API.define_parameter(b"density", 1.0)
    REMAT.API.define_parameter(b"youngs_modulus", youngs_modulus)
    REMAT.API.define_parameter(b"poissons_ratio", 0.25)
    REMAT.API.define_parameter(b"truss_density", 1.0)
    REMAT.API.define_parameter(b"truss_youngs_modulus", youngs_modulus)
    REMAT.API.define_parameter(b"area", 1.0)
    REMAT.API.define_parameter(b"relaxation_time", tau)
    REMAT.API.define_parameter(b"mat_overflow_limit", mat_overflow_limit)


def create_chain_problem(num_elements, impact_velocity, lumped_point_mass):
    num_nodes = num_elements + 1

    coordinates = np.zeros((num_nodes, 2), dtype=np.double)
    velocities = np.zeros((num_nodes, 2), dtype=np.double)
    fixity = np.zeros((num_nodes, 2), dtype=np.bool_)

    for i in range(num_nodes):
        coordinates[i, 0] = float(i)
        coordinates[i, 1] = 0.0

    # Left node fixed, all y-DoFs fixed for a 1D axial chain.
    fixity[0, :] = True
    fixity[:, 1] = True

    # Impact-like loading using initial velocity on the rightmost node.
    velocities[-1, 0] = impact_velocity

    connectivity = np.zeros((0, 4), dtype=np.int32)
    contacts = []
    truss_connectivity = np.zeros((num_elements, 2), dtype=np.int32)
    for e in range(num_elements):
        truss_connectivity[e, 0] = e
        truss_connectivity[e, 1] = e + 1

    REMAT.create_geometry(
        coordinates,
        velocities,
        fixity,
        connectivity,
        contacts,
        truss_connectivity,
    )

    if lumped_point_mass > 0.0:
        point_ids = np.array([num_nodes - 1], dtype=np.int32)
        point_mass = np.array([lumped_point_mass], dtype=np.double)
        REMAT.API.define_point_mass(point_ids, point_mass, 1)


def run_forward_only(
    *,
    num_elements,
    dt,
    nsteps,
    nsub_steps,
    tau,
    youngs_modulus,
    impact_velocity,
    lumped_point_mass,
    mat_overflow_limit,
    integrator_type,
    store_histories=False,
):
    REMAT.API.set_integrator_type(integrator_type)
    define_parameters(tau, youngs_modulus, mat_overflow_limit)
    create_chain_problem(num_elements, impact_velocity, lumped_point_mass)
    REMAT.API.initialize()

    loss = 0.0
    time_history = []
    mean_stress_history = []

    for _ in range(nsteps):
        # Built-in adjoint recurrence corresponds to summing stress at the current state
        # before the next constitutive update.
        stress = REMAT.get_field(b"truss", "axial_stress")
        loss += 0.5 * np.sum((stress * stress) / youngs_modulus)

        if store_histories:
            time_history.append(REMAT.API.get_time())
            mean_stress_history.append(float(np.mean(stress)))

        REMAT.API.update_state(+dt, nsub_steps)

    result = {"loss": float(loss)}
    if store_histories:
        result["time_history"] = np.asarray(time_history)
        result["mean_stress_history"] = np.asarray(mean_stress_history)
    return result


def run_forward_backward(
    *,
    num_elements,
    dt,
    nsteps,
    nsub_steps,
    tau,
    youngs_modulus,
    impact_velocity,
    lumped_point_mass,
    mat_overflow_limit,
    integrator_type,
    store_histories=False,
):
    REMAT.API.set_integrator_type(integrator_type)
    define_parameters(tau, youngs_modulus, mat_overflow_limit)
    create_chain_problem(num_elements, impact_velocity, lumped_point_mass)
    REMAT.API.initialize()

    loss = 0.0
    time_history = []
    mean_stress_history = []

    for _ in range(nsteps):
        stress = REMAT.get_field(b"truss", "axial_stress")
        loss += 0.5 * np.sum((stress * stress) / youngs_modulus)

        if store_histories:
            time_history.append(REMAT.API.get_time())
            mean_stress_history.append(float(np.mean(stress)))

        REMAT.API.update_state(+dt, nsub_steps)

    for _ in range(nsteps):
        REMAT.API.update_state(-dt, nsub_steps)

    grad_tau = float(np.sum(REMAT.get_field(b"truss", "df_dtau")))

    result = {
        "loss": float(loss),
        "grad_tau": grad_tau,
    }
    if store_histories:
        result["time_history"] = np.asarray(time_history)
        result["mean_stress_history"] = np.asarray(mean_stress_history)
    return result


def finite_difference_gradients(
    *,
    num_elements,
    dt,
    nsteps,
    nsub_steps,
    tau,
    youngs_modulus,
    impact_velocity,
    lumped_point_mass,
    mat_overflow_limit,
    fd_step_tau,
    integrator_type,
):
    plus_tau = run_forward_only(
        num_elements=num_elements,
        dt=dt,
        nsteps=nsteps,
        nsub_steps=nsub_steps,
        tau=tau + fd_step_tau,
        youngs_modulus=youngs_modulus,
        impact_velocity=impact_velocity,
        lumped_point_mass=lumped_point_mass,
        mat_overflow_limit=mat_overflow_limit,
        integrator_type=integrator_type,
    )["loss"]
    minus_tau = run_forward_only(
        num_elements=num_elements,
        dt=dt,
        nsteps=nsteps,
        nsub_steps=nsub_steps,
        tau=tau - fd_step_tau,
        youngs_modulus=youngs_modulus,
        impact_velocity=impact_velocity,
        lumped_point_mass=lumped_point_mass,
        mat_overflow_limit=mat_overflow_limit,
        integrator_type=integrator_type,
    )["loss"]
    fd_tau = (plus_tau - minus_tau) / (2.0 * fd_step_tau)
    return float(fd_tau)


def optimize_case(
    *,
    case_name,
    num_elements,
    dt,
    nsteps,
    nsub_steps,
    tau0,
    youngs_modulus,
    max_iters,
    lr_tau,
    min_tau,
    optimize_tau,
    impact_velocity,
    lumped_point_mass,
    mat_overflow_limit,
    integrator_type,
    fd_check,
    fd_step_tau,
    store_histories=False,
):
    tau = tau0

    hist_iter = []
    hist_loss = []
    hist_tau = []
    hist_grad_tau = []

    first_history = None
    last_history = None

    print("-" * 90)
    print(f"Case: {case_name} | elements={num_elements}")
    print("-" * 90)

    for it in range(max_iters):
        run = run_forward_backward(
            num_elements=num_elements,
            dt=dt,
            nsteps=nsteps,
            nsub_steps=nsub_steps,
            tau=tau,
            youngs_modulus=youngs_modulus,
            impact_velocity=impact_velocity,
            lumped_point_mass=lumped_point_mass,
            mat_overflow_limit=mat_overflow_limit,
            integrator_type=integrator_type,
            store_histories=store_histories and (it == 0 or it == max_iters - 1),
        )

        if store_histories and it == 0:
            first_history = {
                "time_history": run["time_history"],
                "mean_stress_history": run["mean_stress_history"],
            }
        if store_histories and it == max_iters - 1:
            last_history = {
                "time_history": run["time_history"],
                "mean_stress_history": run["mean_stress_history"],
            }

        loss = run["loss"]
        grad_tau = run["grad_tau"]

        hist_iter.append(it)
        hist_loss.append(loss)
        hist_tau.append(tau)
        hist_grad_tau.append(grad_tau)

        print(
            f"iter={it:03d}  loss={loss: .8e}  tau={tau: .8e}  grad_tau={grad_tau: .8e}"
        )

        if fd_check and it == 0:
            fd_tau = finite_difference_gradients(
                num_elements=num_elements,
                dt=dt,
                nsteps=nsteps,
                nsub_steps=nsub_steps,
                tau=tau,
                youngs_modulus=youngs_modulus,
                impact_velocity=impact_velocity,
                lumped_point_mass=lumped_point_mass,
                mat_overflow_limit=mat_overflow_limit,
                fd_step_tau=fd_step_tau,
                integrator_type=integrator_type,
            )

            rel_tau = abs(grad_tau - fd_tau) / max(abs(fd_tau), 1.0e-14)
            print(f"FD check: fd_tau={fd_tau: .8e}, rel_err_tau={rel_tau: .3e}")

        if optimize_tau:
            tau = max(min_tau, tau - lr_tau * grad_tau)

    result = {
        "case_name": case_name,
        "num_elements": num_elements,
        "iters": np.asarray(hist_iter),
        "loss": np.asarray(hist_loss),
        "tau": np.asarray(hist_tau),
        "grad_tau": np.asarray(hist_grad_tau),
        "first_history": first_history,
        "last_history": last_history,
    }
    return result


def plot_optimization_result(result, out_prefix):
    case_name = result["case_name"].replace(" ", "_").lower()

    fig1, ax1 = plt.subplots(figsize=(6.0, 4.0))
    ax1.plot(result["iters"], result["loss"], linewidth=1.6, color="#2b738e")
    ax1.set_xlabel("optimization iteration", fontsize="large")
    ax1.set_ylabel("loss", fontsize="large")
    ax1.set_title(f"{result['case_name']}: loss history", fontsize="medium")
    fig1.tight_layout()
    fig1.savefig(f"{out_prefix}_{case_name}_loss.svg", dpi=200)
    plt.close(fig1)

    fig2, ax_tau = plt.subplots(figsize=(6.0, 3.6))
    ax_tau.plot(result["iters"], result["tau"], linewidth=1.6, color="#f9826b")
    ax_tau.set_xlabel("optimization iteration", fontsize="large")
    ax_tau.set_ylabel(r"$\\tau$", fontsize="large")
    ax_tau.set_title(f"{result['case_name']}: parameter history", fontsize="medium")
    fig2.tight_layout()
    fig2.savefig(f"{out_prefix}_{case_name}_params.svg", dpi=200)
    plt.close(fig2)

    if result["first_history"] is not None and result["last_history"] is not None:
        fig3, ax3 = plt.subplots(figsize=(6.0, 4.0))
        ax3.plot(
            result["first_history"]["time_history"],
            result["first_history"]["mean_stress_history"],
            linewidth=1.4,
            label="iteration 0",
            color="#2b738e",
        )
        ax3.plot(
            result["last_history"]["time_history"],
            result["last_history"]["mean_stress_history"],
            linewidth=1.4,
            label=f"iteration {len(result['iters']) - 1}",
            color="#f9826b",
        )
        ax3.set_xlabel("time (s)", fontsize="large")
        ax3.set_ylabel("mean axial stress", fontsize="large")
        ax3.set_title(f"{result['case_name']}: forward stress trace", fontsize="medium")
        ax3.legend(loc="best", fontsize="medium")
        fig3.tight_layout()
        fig3.savefig(f"{out_prefix}_{case_name}_stress.svg", dpi=200)
        plt.close(fig3)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Adjoint-based viscoelastic optimization for tau (1-element and multi-element truss chains)."
    )
    parser.add_argument("--dt", type=float, default=1.0e-3)
    parser.add_argument("--nsteps", type=int, default=200)
    parser.add_argument("--nsub-steps", type=int, default=1)

    parser.add_argument("--tau0", type=float, default=0.30)
    parser.add_argument("--youngs-modulus", type=float, default=10.0)
    parser.add_argument("--max-iters", type=int, default=50)
    parser.add_argument("--lr-tau", type=float, default=1.0e-3)
    parser.add_argument("--min-tau", type=float, default=1.0e-4)
    parser.add_argument("--disable-optimize-tau", action="store_true")

    parser.add_argument("--impact-velocity", type=float, default=2.0)
    parser.add_argument("--point-mass", type=float, default=1.0)
    parser.add_argument("--mat-overflow-limit", type=float, default=1.0e6)

    parser.add_argument("--num-elements-one", type=int, default=1)
    parser.add_argument("--num-elements-many", type=int, default=8)
    parser.add_argument("--skip-one", action="store_true")
    parser.add_argument("--skip-many", action="store_true")

    parser.add_argument("--integrator-type", type=str, default="fixed_truss_visco_adj_float")
    parser.add_argument("--fd-check", action="store_true")
    # Fixed-point truss mode benefits from larger FD perturbations.
    parser.add_argument("--fd-step-tau", type=float, default=1.0e-2)

    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--plot-prefix", type=str, default="adjoint_visco_optimize")
    return parser.parse_args()


def main():
    args = parse_args()

    optimize_tau = not args.disable_optimize_tau
    if not optimize_tau:
        raise ValueError("Tau optimization must be enabled.")

    case_results = []

    if not args.skip_one:
        one_case = optimize_case(
            case_name="One Element + Lumped Mass",
            num_elements=args.num_elements_one,
            dt=args.dt,
            nsteps=args.nsteps,
            nsub_steps=args.nsub_steps,
            tau0=args.tau0,
            youngs_modulus=args.youngs_modulus,
            max_iters=args.max_iters,
            lr_tau=args.lr_tau,
            min_tau=args.min_tau,
            optimize_tau=optimize_tau,
            impact_velocity=args.impact_velocity,
            lumped_point_mass=args.point_mass,
            mat_overflow_limit=args.mat_overflow_limit,
            integrator_type=args.integrator_type.encode("utf-8"),
            fd_check=args.fd_check,
            fd_step_tau=args.fd_step_tau,
            store_histories=args.plot,
        )
        case_results.append(one_case)

    if not args.skip_many:
        many_case = optimize_case(
            case_name="Many Elements + Lumped Mass",
            num_elements=args.num_elements_many,
            dt=args.dt,
            nsteps=args.nsteps,
            nsub_steps=args.nsub_steps,
            tau0=args.tau0,
            youngs_modulus=args.youngs_modulus,
            max_iters=args.max_iters,
            lr_tau=args.lr_tau,
            min_tau=args.min_tau,
            optimize_tau=optimize_tau,
            impact_velocity=args.impact_velocity,
            lumped_point_mass=args.point_mass,
            mat_overflow_limit=args.mat_overflow_limit,
            integrator_type=args.integrator_type.encode("utf-8"),
            fd_check=args.fd_check,
            fd_step_tau=args.fd_step_tau,
            store_histories=args.plot,
        )
        case_results.append(many_case)

    print("=" * 90)
    print("Final summary")
    print("=" * 90)
    for result in case_results:
        print(
            f"{result['case_name']}: "
            f"loss={result['loss'][-1]: .8e}, "
            f"tau={result['tau'][-1]: .8e}"
        )

    if args.plot:
        for result in case_results:
            plot_optimization_result(result, args.plot_prefix)


if __name__ == "__main__":
    main()
