# AGENTS.md - peptide-tag

Canonical context for AI coding agents. (Claude Code reads `CLAUDE.md`, which imports this.)

## Project
Type your name, get the peptide it spells, and find where that sequence turns up
inside a real human protein. Flask, server rendered.

The product promise is that everything on the page is **true**. An earlier
version generated fake "therapeutic applications" attributed to a source called
"AI Analysis", alongside a molecular formula that contradicted the molecular
weight printed beside it. Do not reintroduce anything in that spirit. If a fact
is not computed or looked up, it does not go on the page.

## Layout
- `app.py` - routes only.
- `chemistry.py` - residue table, name to sequence, SMILES, formula, pI, hydropathy.
- `proteome.py` - loads reviewed human Swiss-Prot and finds the longest run of a
  name inside it.
- `artwork.py` - the residue-chain drawing: inline SVG, specimen card, social PNG.
- `structure.py` - ideal-geometry 3D model (PeptideBuilder).
- `printing.py` - watertight STL for 3D printing (trimesh).
- `data/human_swissprot.fasta.gz` - 7.2 MB, committed on purpose. Refresh with
  `python -m proteome --refresh`.
- `static/fonts/` - two OFL fonts, committed because the slim container has none
  and Pillow needs real files to draw the social card.

## Routes
`/`, `/p/<name>`, `/p/<name>/tag.svg`, `/p/<name>/og.png`, `/p/<name>/model.pdb`,
`/p/<name>/model.stl`, `/healthz`.

Results are a **GET**. They were a POST once, which meant no shareable link and
no social preview, and the tag lived in a process-local dict that 404ed on every
restart. Keep results addressable and stateless.

## Commands
- Install: `pip install -r requirements-dev.txt`
- Test: `.venv/bin/python -m pytest tests/ -q`
- Run: `python app.py` (port 5055; macOS AirPlay occupies 5000)
- Serve: gunicorn per `Procfile` / `Dockerfile`

## Deploy
Railway, built from the `Dockerfile`. Push to `main`; if Railway is not tracking
latest, force-deploy the main HEAD via GraphQL `serviceInstanceDeployV2` with the
commit SHA. `SHOP_URL` is optional: set it and an "order a print" link appears.

## Conventions (IMPORTANT)
- **NEVER use em dashes.**
- **No apt packages in the Dockerfile.** Every dependency ships a manylinux
  wheel and the fonts are committed. In particular do not reintroduce the old
  RDKit apt block (`libgl1-mesa-glx`); that package is gone from Debian bookworm
  and it broke the build.
- **RDKit is a test-only dependency** (`requirements-dev.txt`). It independently
  rebuilds every peptide and checks our SMILES and formulas against its own. Do
  not import it from application code: the app must not carry 128 MB to compute
  a formula.
- **Do not use `<path:name>` in routes.** It swallows `/model.pdb` and renders it
  as a name. The default string converter is deliberate.
- Chemistry changes need a passing `tests/test_chemistry.py`. It cross-checks
  against RDKit, so a wrong side chain or a wrong water-loss correction fails.
- Secrets in `.env` only (gitignored), never committed.

## Analytics
PostHog, cookieless, EU, shared "Robin Apps" project, split by the `app`
super-property. The snippet lives in `templates/base.html`. Keep it on both pages.

## Do-not
- No `libgl1-mesa-glx` or any other apt line. No committed secrets. No em dashes.
- No invented biology, no therapeutic claims, no fabricated citations.
