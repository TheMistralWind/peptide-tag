"""Peptide Tag: your name, as a molecule.

Every result page is a plain GET at /p/<name>, so a peptide can be linked,
bookmarked and previewed. The previous version answered a POST and cached the
artwork in a process-local dict, which meant no shareable URL at all and a 404
on every download link the moment the container restarted.
"""

from __future__ import annotations

import os
import re
from urllib.parse import quote

from flask import (Flask, Response, abort, make_response, redirect,
                   render_template, request, url_for)

import artwork
import chemistry as C
import printing
import proteome
import structure

app = Flask(__name__)

# Parse the proteome at import rather than on the first request. Under gunicorn
# --preload this happens once in the parent, so workers share the 11 MB blob
# copy-on-write and nobody's first search pays the 0.1s.
proteome.get()

MAX_INPUT = 80
CACHE = "public, max-age=31536000, immutable"


def canonical(name: str) -> str:
    """Collapse whitespace and cap length, so URLs stay tidy and stable.

    Slashes are dropped rather than escaped. The routes use the default string
    converter precisely so that /p/Robin/model.pdb resolves to the model and not
    to a peptide named "Robin/model.pdb", which is what a path converter did.
    """
    cleaned = re.sub(r"[/\\]", " ", name or "")
    return re.sub(r"\s+", " ", cleaned.strip())[:MAX_INPUT]


def sentence(seq: str, props: dict) -> str:
    """One true sentence about this peptide, computed, never invented."""
    n = props["length"]
    pi = props["isoelectric_point"]
    g = props["gravy"]

    if pi > 8.5:
        charge = ("It carries a net positive charge in blood, which is the sort "
                  "of thing that makes a peptide stick to DNA and to cell membranes")
    elif pi < 6.0:
        charge = ("It carries a net negative charge in blood, which is the sort "
                  "of thing that makes a peptide bind metal ions")
    else:
        charge = ("It sits close to neutral at the pH of blood, which is also "
                  "the pH at which it would be least soluble")

    if g > 0.5:
        water = "It is greasy enough that it would rather sit in a membrane than in water"
    elif g < -1.0:
        water = "It is water-loving, and would dissolve readily"
    else:
        water = "It is neither strongly water-loving nor greasy"

    if n <= 5:
        size = f"At {n} residues it is far too short to fold into anything"
    elif n <= 15:
        size = (f"At {n} residues it is still too short to fold on its own, "
                "though peptides this size are exactly the size of most hormones")
    else:
        size = (f"At {n} residues it is long enough that it might adopt "
                "a shape, given something to fold against")

    return f"{size}. {charge}. {water}."


def build(name: str) -> dict | None:
    """Everything the result page needs. None when there is nothing to show."""
    peptide = C.text_to_sequence(name)
    if not peptide.sequence:
        return None

    seq = peptide.sequence
    props = C.analyse(seq)
    standard = C.to_standard(seq)
    fragment = proteome.longest_fragment(standard)

    residues = C.composition(seq)
    for row in residues:
        row["charge"] = C.CHARGE.get(row["letter"], 0)
        row["partial"] = row["letter"] in C.PARTIAL_CHARGE

    return {
        "name": name,
        "peptide": peptide,
        "seq": seq,
        "props": props,
        "residues": residues,
        "sentence": sentence(seq, props),
        "smiles": C.build_smiles(seq),
        "fragment": fragment,
        "projected": standard != seq,
        "chain_svg": artwork.chain_svg(seq, max_width=760, max_cols=10),
        "chain_svg_narrow": artwork.chain_svg(seq, max_width=330, max_cols=6),
    }


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        name = canonical(request.form.get("name") or request.form.get("username", ""))
        if not C.text_to_sequence(name).sequence:
            return render_template(
                "index.html", error="That has no letters in it, so there is "
                                    "nothing to spell a peptide with.",
                value=name), 400
        return redirect(url_for("peptide", name=name))
    return render_template("index.html")


@app.route("/p/<name>")
def peptide(name: str):
    clean = canonical(name)
    if not clean:
        return redirect(url_for("index"))
    if clean != name:
        return redirect(url_for("peptide", name=clean))
    data = build(clean)
    if not data:
        return render_template(
            "index.html", error="That has no letters in it, so there is "
                                "nothing to spell a peptide with.",
            value=clean), 404
    return render_template("result.html", **data)


@app.route("/p/<name>/tag.svg")
def tag_svg(name: str):
    data = build(canonical(name))
    if not data:
        abort(404)
    svg = artwork.specimen_svg(data["name"], data["seq"], data["props"])
    r = make_response(svg)
    r.headers["Content-Type"] = "image/svg+xml; charset=utf-8"
    r.headers["Cache-Control"] = CACHE
    if request.args.get("download"):
        safe = re.sub(r"[^A-Za-z0-9_-]", "_", data["name"]) or "peptide"
        r.headers["Content-Disposition"] = f'attachment; filename="{safe}.svg"'
    return r


@app.route("/p/<name>/og.png")
def og_png(name: str):
    data = build(canonical(name))
    if not data:
        abort(404)
    png = artwork.social_png(data["name"], data["seq"], data["props"])
    return Response(png, mimetype="image/png",
                    headers={"Cache-Control": CACHE})


@app.route("/p/<name>/model.pdb")
def model_pdb(name: str):
    shape = request.args.get("shape", "helix")
    if shape not in structure.SHAPES:
        shape = "helix"
    data = build(canonical(name))
    if not data:
        abort(404)
    pdb = structure.model_pdb(data["seq"], shape)
    if not pdb:
        abort(404)
    return Response(pdb, mimetype="chemical/x-pdb",
                    headers={"Cache-Control": CACHE})


@app.route("/p/<name>/model.stl")
def model_stl(name: str):
    shape = request.args.get("shape", "helix")
    if shape not in structure.SHAPES:
        shape = "helix"
    data = build(canonical(name))
    if not data:
        abort(404)
    stl = printing.model_stl(data["seq"], shape)
    if not stl:
        abort(404)
    safe = re.sub(r"[^A-Za-z0-9_-]", "_", data["name"]) or "peptide"
    return Response(stl, mimetype="model/stl", headers={
        "Cache-Control": CACHE,
        "Content-Disposition": f'attachment; filename="{safe}-{shape}.stl"',
    })


@app.route("/healthz")
def healthz():
    p = proteome.get()
    return {"ok": True, "proteins": len(p.proteins), "residues": p.residues}


@app.errorhandler(404)
def not_found(_):
    return render_template("index.html",
                           error="There is nothing at that address."), 404


@app.errorhandler(500)
def server_error(_):
    return render_template("index.html",
                           error="Something broke on our side. Try again."), 500


@app.template_filter("formula")
def formula_filter(value: str) -> str:
    from markupsafe import Markup
    return Markup(C.format_formula_html(value))


@app.context_processor
def helpers():
    # SHOP_URL is unset until there is a shop. The "order a print" link simply
    # does not render until Railway has the variable, so nothing here promises
    # something that cannot be bought yet.
    return {"quote": quote, "shop_url": os.environ.get("SHOP_URL", "")}


if __name__ == "__main__":
    proteome.get()  # warm the index before serving
    # 5000 is taken by AirPlay on macOS, so local runs default elsewhere.
    port = int(os.environ.get("PORT", 5055))
    app.run(host="0.0.0.0", port=port, debug=False)
