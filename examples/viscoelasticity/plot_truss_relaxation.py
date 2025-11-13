from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pyexodus

def load_all_state_histories(state_names, exo_path):
    # Load all requested state variables in a single pass from the Exodus file.
    with pyexodus.exodus(str(exo_path), mode="r", array_type="numpy") as exo:
        ds = exo._f
        time_steps = np.asarray(ds.variables["time_whole"][:-1], dtype=float)
        if time_steps.size == 0:
            raise RuntimeError("No time steps found in the Exodus file.")
        
        block_id = 1
        num_steps = time_steps.size
        
        # Initialize dictionary to store all state data
        state_data = {name: np.empty(num_steps, dtype=float) for name in state_names}
        
        # Load all state variables for all time steps
        for step_index in range(1, num_steps + 1):
            for state_name in state_names:
                state_data[state_name][step_index - 1] = exo.get_element_variable_values(
                    block_id, state_name, step_index
                )[0]
        
    return time_steps, state_data

def plot_state(time_steps, stress, state_name):
    time_steps = time_steps - 1
    center_x = 0.5 * (time_steps.min() + time_steps.max())
    right = time_steps >= center_x
    left = time_steps <= center_x

    y0, y1 = np.nanmin(stress[left]), np.nanmax(stress[left])
    pad = 0.05 * (y1 - y0 if y1 > y0 else 1.0)
    y0, y1 = y0 - pad, y1 + pad

    fig, (ax1, ax2) = plt.subplots(
        1, 2, sharey=True, figsize=(9, 4),
        gridspec_kw={"wspace": 0.00}
    )

    ax1.plot(time_steps[left], stress[left], marker="o", linewidth=1.2, markersize=2)
    
    # Set ylabel based on state name
    if state_name == "axial_stress":
        ax1.set_ylabel(r"stress ($\sigma_{xx}$)")
    elif state_name == "axial_strain":
        ax1.set_ylabel(r"strain ($\varepsilon_{xx}$)")
    elif state_name == "viscous_strain":
        ax1.set_ylabel(r"viscous strain ($\varepsilon_{xx}^{v}$)")
    elif state_name == "dual_viscous_strain":
        ax1.set_ylabel(r"dual variable")

    ax1.set_xlabel("forward time step")
    ax1.set_ylim(y0, y1)

    ax2.plot(time_steps[left][::-1], stress[right], marker="o", linewidth=1.2, markersize=2)
    ax2.set_xlabel("backward time step")
    ax2.xaxis.set_inverted(True) 

    fig.savefig(f"truss_step_relaxation_{state_name}.pdf", dpi=200)

def plot_state_covering(time_steps, stress, state_name):
    time_steps = time_steps - 1
    center_x = 0.5 * (time_steps.min() + time_steps.max())
    right = time_steps > center_x
    left = time_steps < center_x

    y0, y1 = np.nanmin(stress[left]), np.nanmax(stress[left])
    pad = 0.15 * (y1 - y0 if y1 > y0 else 1.0)
    y0, y1 = y0 - pad, y1 + pad
    


    fig, ax1 = plt.subplots()
    ax1.set_ylim(y0, y1)
    
    color1 = 'k'  # Blue for forward
    color2 = 'tab:red'  # Orange for backward
    line_fwd, = ax1.plot(time_steps[left], stress[left], color=color1, label='forward \n↦⟶⟶⟶')
    line_bwd, = ax1.plot(time_steps[left][::-1], stress[right], color=color2, linestyle='--', label='backward \n⟵⟵⟵↤')


    ax1.set_xlabel('time step', fontsize = 'large', color=color1)
    # ax1.set_ylabel('Absorbance (a.u.)', fontsize = 'large')
    if state_name == "axial_stress":
        ax1.set_ylabel(r"stress ($\sigma_{xx}$)", fontsize = 'large')
    elif state_name == "axial_strain":
        ax1.set_ylabel(r"strain ($\varepsilon_{xx}$)", fontsize = 'large')
    elif state_name == "viscous_strain":
        ax1.set_ylabel(r"viscous strain ($\varepsilon_{xx}^{v}$)", fontsize = 'large')
    elif state_name == "dual_viscous_strain":
        ax1.set_ylabel(r"dual variable", fontsize = 'large')
    ax1.tick_params(axis='x', colors=color1)

   
    fwd_legend = ax1.legend(handles=[line_fwd], loc='upper left')
    fwd_legend.get_texts()[0].set_color(color1)   # "forward →"
    # Add the legend manually to the Axes.
    ax1.add_artist(fwd_legend)
    leg = ax1.legend(handles=[line_bwd], loc='lower right')
    leg.get_texts()[0].set_color(color2)   # "backward ←"


    plt.tight_layout()
    fig.savefig(f"truss_step_relaxation_{state_name}.pdf", dpi=200)

exo_path = "truss_step_relaxation.exo"

# Request list of all states to plot
state_names = ["axial_stress", "axial_strain", "viscous_strain", "dual_viscous_strain"]

# Load all state histories at once
time_steps, all_state_data = load_all_state_histories(state_names, exo_path)

# Loop through each state and create a plot
for state_name in state_names:
    plot_state_covering(time_steps, all_state_data[state_name], state_name)


