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


STATE_TO_PLOT = "axial_strain"

SCENARIOS = [

    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 15000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_step",
        "overflow_limit": 1e6,
    },
        {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 15000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_constant_rate",
        "overflow_limit": 1e6,
    },
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 15000,
        "Nsub_steps": 1,
        "epsilon0": 0.1,
        "bc_name": "right_node_sinusoidal",
        "overflow_limit": 1e6,
    },
    {
        "relaxation_time": 0.3,
        "dt": 1.0e-3,
        "Nsteps": 15000,
        "Nsub_steps": 1,
        "epsilon0": 0.2,
        "bc_name": "right_node_clipped_sinusoid",
        "overflow_limit": 1e6,
    },
]


def format_tau(value):
    return rf"$\tau = {value:.2f}$"


def main():

    for params in SCENARIOS:
        fig, ax = plt.subplots(figsize=(3,3))
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

        ax.plot(
            steps,
            history,
            linewidth=2.0,
            marker=None,
            label=format_tau(params["relaxation_time"]), color='#2b738eff'
        )
        # align= right for the x axis label :
        ax.set_xlabel("time(s)", fontsize="large", horizontalalignment='left', x=0.001)
        ax.set_ylabel(r"axial strain ($\varepsilon_{xx}$)", fontsize="large")
        # ax.legend(loc="upper right", fontsize=10)
        axis_dt = SCENARIOS[0]["dt"]
        ax.set_xlim(0,15)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_position('zero')
        ax.spines['bottom'].set_position(('data', 0)) #'zero')
        # ax.set_yscale("log")
        if params["bc_name"] == "right_node_step":
            ax.set_ylim(0,.12)
        elif params["bc_name"] == "right_node_constant_rate":
            ax.set_ylim(0,params["epsilon0"]*params["Nsteps"]*params["dt"]*1.1)
        else:
            ax.set_ylim(-0.12,0.12)
            # clear the x axis tick  just the zero tick at first
            ticks = ax.xaxis.get_major_ticks()
            for tick in ax.get_xticklabels():
                if tick.get_text() == '0':
                    tick.set_visible(False)
        

        fig.tight_layout()
        fig.savefig(f"strain_{params['bc_name'].replace('right_node', '')}.pdf", dpi=200)
        fig.clf()


if __name__ == "__main__":
    main()
