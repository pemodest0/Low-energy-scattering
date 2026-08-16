# infra — reproducible machines for the lab

One script sets up every machine that runs this laboratory. It is idempotent:
run it again any time.

| file | what it does |
|---|---|
| `setup_lab.sh` | Ubuntu/Debian installer. Auto-detects the target (`mac`, `wsl`, `linux`). Installs the toolchain, clones the repo, builds the Python environment, **runs the 25 tests as a burn-in check**, then sets up SSH, Tailscale, Claude Code and boot-time services. |
| `CHECKLIST_MAC.md` | Step-by-step for turning an Intel MacBook running Ubuntu 24.04 into the 24/7 bench, including the purple-screen recovery path. |
| `setup_desktop_windows.ps1` | Mirrors the lab on a Windows desktop via WSL2 Ubuntu, plus Tailscale to reach the bench. |

## Quick start

**Linux / Mac-with-Ubuntu**

```bash
curl -fsSL https://raw.githubusercontent.com/pemodest0/Low-energy-scattering/main/infra/setup_lab.sh -o /tmp/setup_lab.sh
bash /tmp/setup_lab.sh            # add --role=mac to force the Mac tweaks
```

**Windows desktop** — Terminal as Administrator:

```powershell
Set-ExecutionPolicy -Scope Process Bypass -Force
.\setup_desktop_windows.ps1
```

## Flags

| flag | effect |
|---|---|
| `--role=mac\|wsl\|linux` | skip auto-detection |
| `--no-tailscale` | do not install the private network |
| `--repo-only` | only clone the lab, build the venv and run the tests |

## Shell shortcuts installed

| alias | action |
|---|---|
| `lab` | cd into the lab with the venv activated |
| `lab-test` | run the 25 tests |
| `lab-run` | `python main.py --sem-ajustes` — regenerate figures and tables |
| `lab-app` | Streamlit lab on :8501 |
| `lab-jup` | JupyterLab on :8888 |
| `lab-up` | pull, reinstall deps, re-test |

## The three machines

| machine | role |
|---|---|
| Intel MacBook + Ubuntu 24.04 | the 24/7 bench, reachable over Tailscale |
| Windows desktop + WSL2 | local mirror / workstation |
| Heaviside (IFSC-USP) | heavy jobs |

## Why Ubuntu 24.04 LTS and not Mint

The purple screen after a Mint update is a new kernel fighting the Mac's
graphics. A *definitive* machine wants the opposite of surprise: 5 years of
conservative updates, and the distribution every scientific-computing tutorial
assumes. The installer keeps old kernels around and leaves the GRUB menu
visible for 5 seconds — that menu is the escape hatch if a future kernel
misbehaves.
