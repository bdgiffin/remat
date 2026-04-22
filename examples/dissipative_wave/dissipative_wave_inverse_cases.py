import argparse
import shlex
import subprocess
import sys
from pathlib import Path


THIS_DIR = Path(__file__).resolve().parent
INVERSE_SCRIPT = THIS_DIR / "dissipative_wave_inverse.py"
OUTPUT_DIR = THIS_DIR / "inverse_case_outputs"


# Five preset inverse problems with different difficulty/conditioning tradeoffs.
CASES = {
    "case1_quick": {
        "description": "Fast coarse run to sanity-check inversion behavior.",
        "args": [
            "--nx", "22", "--ny", "8",
            "--nsteps", "70", "--nsub-steps", "1", "--dt", "2.2e-3",
            "--n-layers", "4", "--n-sensors", "8",
            "--impact-velocity", "0.75", "--source-window-fraction", "0.10",
            "--true-layers", "0.80,1.40,0.72,1.28",
            "--init-layers", "1.22,0.86,1.12,0.88",
            "--true-tau", "0.10", "--init-tau", "0.28",
            "--max-iters", "30",
            "--lr-layers", "18.0", "--lr-tau", "0.08",
            "--min-layer", "0.2", "--max-layer", "2.2",
            "--min-tau", "0.01", "--max-tau", "0.6",
        ],
    },
    "case2_balanced": {
        "description": "Balanced mid-size mesh and longer horizon.",
        "args": [
            "--nx", "36", "--ny", "12",
            "--nsteps", "130", "--nsub-steps", "1", "--dt", "1.7e-3",
            "--n-layers", "4", "--n-sensors", "12",
            "--impact-velocity", "0.85", "--source-window-fraction", "0.10",
            "--true-layers", "0.78,1.48,0.70,1.30",
            "--init-layers", "1.18,0.90,1.08,0.92",
            "--true-tau", "0.09", "--init-tau", "0.26",
            "--max-iters", "45",
            "--lr-layers", "0.10", "--lr-tau", "0.04",
            "--min-layer", "0.2", "--max-layer", "2.2",
            "--min-tau", "0.01", "--max-tau", "0.6",
        ],
    },
    "case3_high_contrast": {
        "description": "Higher-contrast layered profile for stronger reflections (with L-BFGS-B boundary restarts).",
        "args": [
            "--nx", "150", "--ny", "30",
            "--nsteps", "1000", "--nsub-steps", "1", "--dt", "4e-3",
            "--n-layers", "2", "--n-sensors", "5",
            "--impact-velocity", "1.0", "--source-window-fraction", "0.09",
            "--true-layers", "15,7.7",
            "--init-layers", "10.0,4.0",
            "--true-tau", "0.08", "--init-tau", "0.1",
            "--max-iters", "25",
            "--min-layer", "1.0", "--max-layer", "20.0",
            "--min-tau", "0.03", "--max-tau", "0.18",
            "--lbfgsb-restarts", "4",
            "--lbfgsb-boundary-perturb", "0.2",
            "--overflow-limit", "10",
            "--mat-overflow-limit", "10",
        ],
    },
    "case3_new": {
        "description": "Higher-contrast layered profile for stronger reflections (with L-BFGS-B boundary restarts).",
        "args": [
            "--nx", "150", "--ny", "30",
            "--nsteps", "1000", "--nsub-steps", "1", "--dt", "4e-3",
            "--n-layers", "4", "--n-sensors", "5",
            "--impact-velocity", "1.0", "--source-window-fraction", "0.09",
            "--true-layers", "15,7.7,4.4,1.8",
            "--init-layers", "10.0,4.0,5.0,3.0",
            "--true-tau", "0.08", "--init-tau", "0.1",
            "--max-iters", "25",
            "--min-layer", "1.0", "--max-layer", "20.0",
            "--min-tau", "0.03", "--max-tau", "0.18",
            "--lbfgsb-restarts", "4",
            "--lbfgsb-boundary-perturb", "0.2",
            "--overflow-limit", "10",
            "--mat-overflow-limit", "10",
        ],
    },
    "case4_tau_focus": {
        "description": "Easier layer field, harder tau mismatch (tau-identification focus).",
        "args": [
            "--nx", "34", "--ny", "12",
            "--nsteps", "140", "--nsub-steps", "1", "--dt", "1.8e-3",
            "--n-layers", "4", "--n-sensors", "10",
            "--impact-velocity", "0.90", "--source-window-fraction", "0.10",
            "--true-layers", "1.00,1.08,0.94,1.00",
            "--init-layers", "1.00,1.00,1.00,1.00",
            "--true-tau", "0.07", "--init-tau", "0.34",
            "--max-iters", "45",
            "--lr-layers", "3.0", "--lr-tau", "1",
            "--min-layer", "0.2", "--max-layer", "2.2",
            "--min-tau", "0.01", "--max-tau", "0.6",
        ],
    },
    "case5_layer_focus": {
        "description": "Tau starts near true value; emphasis on recovering layer coefficients.",
        "args": [
            "--nx", "38", "--ny", "14",
            "--nsteps", "150", "--nsub-steps", "1", "--dt", "1.6e-3",
            "--n-layers", "4", "--n-sensors", "12",
            "--impact-velocity", "0.95", "--source-window-fraction", "0.10",
            "--true-layers", "0.66,1.56,0.74,1.38",
            "--init-layers", "1.24,0.86,1.16,0.90",
            "--true-tau", "0.12", "--init-tau", "0.14",
            "--max-iters", "50",
            "--lr-layers", "12.0", "--lr-tau", "0.01",
            "--min-layer", "0.2", "--max-layer", "2.2",
            "--min-tau", "0.01", "--max-tau", "0.6",
        ],
    },
}


def build_command(case_name, extra_args):
    case = CASES[case_name]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_file = OUTPUT_DIR / f"{case_name}_summary.png"
    setup_plot_file = OUTPUT_DIR / f"{case_name}_setup.svg"

    cmd = [
        sys.executable,
        str(INVERSE_SCRIPT),
        *case["args"],
        "--plot-file",
        str(plot_file),
        "--setup-plot-file",
        str(setup_plot_file),
        *extra_args,
    ]
    return cmd


def print_case_table():
    print("Available inverse presets:\n")
    for name, spec in CASES.items():
        print(f"- {name}: {spec['description']}")
        cmd = build_command(name, [])
        print("  ", shlex.join(cmd))
    print("\nTip: append overrides after '--', e.g.")
    print("  python dissipative_wave_inverse_cases.py --case case2_balanced -- --max-iters 20")


def run_case(case_name, extra_args, dry_run=False):
    cmd = build_command(case_name, extra_args)
    print(f"\n[{case_name}] {CASES[case_name]['description']}")
    print(shlex.join(cmd))
    if not dry_run:
        subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Launcher for multiple dissipative_wave_inverse presets. "
            "Use --list to inspect all commands, then run one with --case."
        )
    )
    parser.add_argument("--list", action="store_true", help="List all preset problems and commands.")
    parser.add_argument("--case", type=str, default="case2_balanced", choices=sorted(CASES.keys()),
                        help="Run one preset case.")
    parser.add_argument("--run-all", action="store_true", help="Run all 5 preset cases sequentially.")
    parser.add_argument("--dry-run", action="store_true", help="Print command(s) only, do not execute.")

    args, extra = parser.parse_known_args()
    if extra and extra[0] == "--":
        extra = extra[1:]

    if args.list:
        print_case_table()
        if not args.run_all:
            return

    if args.run_all:
        for case_name in CASES.keys():
            run_case(case_name, extra, dry_run=args.dry_run)
        return

    run_case(args.case, extra, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
