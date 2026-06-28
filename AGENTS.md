# AGENTS.md — peptide-tag

Canonical context for AI coding agents. (Claude Code reads `CLAUDE.md`, which imports this.)

## Project
A Python (Flask) web app for peptide tagging. Server-rendered (templates + static).

## Stack
- Python (pinned in `runtime.txt`). Flask app in `app.py`. Deps in `requirements.txt`. `templates/` (Jinja) + `static/`.
- Containerized: `Dockerfile` (Debian bookworm base). `Procfile` for the web process.

## Commands
- Install: `pip install -r requirements.txt`
- Run (local): `python app.py` (or gunicorn per `Procfile`).
- Container: build per `Dockerfile`.

## Deploy
- Railway, built from the `Dockerfile`. Push to `main`; if Railway is not tracking latest, force-deploy the main HEAD via GraphQL `serviceInstanceDeployV2` with the commit SHA.

## Conventions (IMPORTANT)
- **Do NOT reintroduce the RDKit apt block** (`libgl1-mesa-glx`). That package was removed from Debian bookworm and it broke the Docker build; it was deliberately deleted. Keep the Dockerfile lean.
- **NEVER use em dashes.**
- Secrets in `.env` only (gitignored), never committed.

## Do-not
- No `libgl1-mesa-glx` / vestigial RDKit apt lines. No committed secrets. No em dashes.
