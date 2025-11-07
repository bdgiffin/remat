# %%
# from pathlib import Path
# import matplotlib.pyplot as plt
# import numpy as np
# import pyexodus

# def load_stress_history(exo_path):
#     with pyexodus.exodus(str(exo_path), mode="r", array_type="numpy") as exo:
#         ds = exo._f
#         time_steps = np.asarray(ds.variables["time_whole"][:-1], dtype=float)
#         if time_steps.size == 0:
#             raise RuntimeError("No time steps found in the Exodus file.")
#         block_id = 1
#         stress = np.empty_like(time_steps, dtype=float)
#         for step_index in range(1, time_steps.size):
#             stress[step_index - 1] = exo.get_element_variable_values(block_id, "axial_stress", step_index)[0]
#     return time_steps, stress

# def plot_stress(time_steps, stress):
#     fig, ax = plt.subplots()
#     time_steps = time_steps - 1

#     ax.plot(time_steps, stress, marker="o", linewidth=1.2, markersize=2)
#     ax.set_xlabel("time step")
#     ax.set_ylabel(r"stress ($\sigma_{xx}$)")

#     # --- center + gradient for x > center WITHOUT changing data limits ---
#     center_x = 0.5 * (time_steps.min() + time_steps.max())
#     y0, y1 = ax.get_ylim()

#     # Draw gradient in axes coords (y from 0 to 1), so limits stay untouched
#     grad = np.linspace(1.0, 0.0, 256)[None, :]  # 1xN row, horizontal fade
#     ax.imshow(
#         grad,
#         extent=[center_x, time_steps.max(), 0,1],  # x in data, y in axes coords
#         transform=ax.get_xaxis_transform(),         # <- key: axes-coord in y
#         origin="lower",
#         cmap="Greys",
#         alpha=0.12,
#         aspect="auto",
#         interpolation="bilinear",
#         zorder=0,
#     )

#     ax.axvline(center_x, color="gray", linestyle="--", linewidth=0.8)

#     # Make sure original data limits are used
#     ax.set_ylim(y0, y1)
#     ax.set_xlim(time_steps.min()-5, time_steps.max()+5)
#     t = ax.text(41, 0.002, "Backward",
#             ha="left", va="center", rotation=0, size=12,
#             bbox=dict(boxstyle="rarrow,pad=0.3",
#                       fc="lightgray", ec="black", lw=.9))

#     fig.tight_layout()
#     fig.savefig("truss_step_relaxation.pdf", dpi=200)

# exo_path = "truss_step_relaxation.exo"
# time_steps, stress = load_stress_history(exo_path)
# plot_stress(time_steps, stress)
#%%
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pyexodus

def load_state_history(state_name,exo_path):
    with pyexodus.exodus(str(exo_path), mode="r", array_type="numpy") as exo:
        ds = exo._f
        time_steps = np.asarray(ds.variables["time_whole"][:-1], dtype=float)
        if time_steps.size == 0:
            raise RuntimeError("No time steps found in the Exodus file.")
        block_id = 1
        state_data = np.empty_like(time_steps, dtype=float)
        for step_index in range(1, time_steps.size+1):
            state_data[step_index - 1] = exo.get_element_variable_values(block_id, state_name, step_index)[0]
        # print(f"Loaded state '{state_data}' from Exodus file.")
    return time_steps, state_data

def plot_stress(time_steps, stress):
    time_steps = time_steps - 1
    center_x = 0.5 * (time_steps.min() + time_steps.max())
    right = time_steps >= center_x
    left = time_steps <= center_x

    y0, y1 = np.nanmin(stress), np.nanmax(stress[left])
    pad = 0.05 * (y1 - y0 if y1 > y0 else 1.0)
    y0, y1 = y0 - pad, y1 + pad

    fig, (ax1, ax2) = plt.subplots(
        1, 2, sharey=True, figsize=(9, 4),
        gridspec_kw={"wspace": 0.00}
    )

    ax1.plot(time_steps[left], stress[left], marker="o", linewidth=1.2, markersize=2)
    ax1.set_ylabel(r"stress ($\sigma_{xx}$)")
    # ax1.set_ylabel(r"dual variable")

    ax1.set_xlabel("forward time step")
    ax1.set_ylim(y0, y1)
    # ax1.set_xlim(time_steps.min(), center_x)
    # ax1.axvline(center_x, color="gray", linestyle="--", linewidth=0.8)


    ax2.plot(time_steps[left][::-1], stress[right], marker="o", linewidth=1.2, markersize=2)
    ax2.set_xlabel("backward time step")
    # ax2.set_ylim(y0, y1)
    ax2.xaxis.set_inverted(True) 
    # ax2.set_xlim(center_x, time_steps.max())
    # ax2.axvline(center_x, color="gray", linestyle="--", linewidth=0.8)

    # fig.tight_layout()
    fig.savefig("truss_step_relaxation.pdf", dpi=200)

exo_path = "truss_step_relaxation.exo"
ts, s = load_state_history("axial_stress", exo_path)


# ts, s = load_state_history("dual_viscous_strain", exo_path)

plot_stress(ts, s)


