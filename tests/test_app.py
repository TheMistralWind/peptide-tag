"""Route-level checks, mostly guarding bugs the old version actually shipped."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app import app as flask_app  # noqa: E402


@pytest.fixture()
def client():
    flask_app.config["TESTING"] = True
    with flask_app.test_client() as c:
        yield c


def test_landing_page(client):
    r = client.get("/")
    assert r.status_code == 200
    assert b"Your name is hiding inside a real protein" in r.data


@pytest.mark.parametrize("junk", ["123", "!!!", "   ", "...", "42"])
def test_letterless_input_does_not_500(client, junk):
    """The old app returned HTTP 500 for any input without letters."""
    r = client.post("/", data={"name": junk})
    assert r.status_code == 400
    assert b"no letters" in r.data


def test_post_redirects_to_a_shareable_url(client):
    """The old app answered the POST directly, so results had no address."""
    r = client.post("/", data={"name": "Robin"})
    assert r.status_code == 302
    assert r.headers["Location"] == "/p/Robin"


def test_result_page_is_a_plain_get(client):
    r = client.get("/p/Robin")
    assert r.status_code == 200
    assert b"RONIN" in r.data
    assert b"752.9 Da" in r.data


def test_result_page_carries_social_tags(client):
    r = client.get("/p/Robin")
    body = r.data.decode()
    assert 'property="og:image"' in body
    assert "/p/Robin/og.png" in body
    assert 'name="twitter:card"' in body


def test_asset_routes_are_not_swallowed_by_the_name_route(client):
    """A path converter matched "Robin/model.pdb" as a name. A string one does not."""
    assert client.get("/p/Robin/model.pdb").mimetype == "chemical/x-pdb"
    assert client.get("/p/Robin/og.png").mimetype == "image/png"
    assert client.get("/p/Robin/tag.svg").mimetype == "image/svg+xml"


def test_tag_svg_has_a_viewbox(client):
    """Without a viewBox the old tag was clipped instead of scaled."""
    body = client.get("/p/Robin/tag.svg").data.decode()
    assert "viewBox=" in body


def test_tag_svg_is_inline_unless_download_is_asked_for(client):
    """It was served as an attachment while also being used as an img src."""
    assert "Content-Disposition" not in client.get("/p/Robin/tag.svg").headers
    r = client.get("/p/Robin/tag.svg?download=1")
    assert r.headers["Content-Disposition"].startswith("attachment")


def test_og_image_is_a_png_not_an_svg(client):
    """No major platform renders SVG in og:image."""
    r = client.get("/p/Robin/og.png")
    assert r.status_code == 200
    assert r.data[:8] == b"\x89PNG\r\n\x1a\n"


def test_helix_and_strand_are_different_structures(client):
    helix = client.get("/p/Robin/model.pdb?shape=helix").data
    strand = client.get("/p/Robin/model.pdb?shape=strand").data
    assert helix and strand and helix != strand


def test_unknown_shape_falls_back_rather_than_erroring(client):
    r = client.get("/p/Robin/model.pdb?shape=../../etc/passwd")
    assert r.status_code == 200


def test_stl_is_binary_and_attached(client):
    r = client.get("/p/Robin/model.stl")
    assert r.status_code == 200
    assert len(r.data) > 10_000
    assert r.headers["Content-Disposition"].startswith("attachment")


def test_names_with_slashes_do_not_break_routing(client):
    r = client.get("/p/Robin%2Fmodel.pdb")
    assert r.status_code in (200, 301, 302, 404)


def test_accented_name_round_trips(client):
    r = client.get("/p/H%C3%A5kan")
    assert r.status_code == 200
    assert b"HAKAN" in r.data


def test_page_does_not_dump_the_smiles_string(client):
    """A raw SMILES in the body is what overflowed the card by 232px."""
    body = client.get("/p/Robin").data.decode()
    assert "Copy SMILES" not in body
    assert "<code>" not in body


def test_no_fabricated_science_remains(client):
    body = client.get("/p/Robin").data.decode().lower()
    for phrase in ["therapeutic potential", "ai analysis", "therapeutic applications",
                   "biological pathways", "predicted function"]:
        assert phrase not in body, phrase


def test_social_urls_respect_the_forwarded_protocol(client):
    """Railway forwards plain HTTP, so without ProxyFix og:image goes out
    as http:// and crawlers are picky about that."""
    r = client.get("/p/Robin", headers={
        "X-Forwarded-Proto": "https",
        "X-Forwarded-Host": "peptide-tag-production.up.railway.app",
    })
    body = r.data.decode()
    assert "https://peptide-tag-production.up.railway.app/p/Robin/og.png" in body
    assert "http://peptide-tag-production" not in body


def test_interest_form_appears_only_while_there_is_no_shop(client, monkeypatch):
    assert b"count me in" in client.get("/p/Robin").data
    monkeypatch.setenv("SHOP_URL", "https://example.com/shop")
    body = client.get("/p/Robin").data
    assert b"count me in" not in body
    assert b"order a print" in body


def test_interest_accepts_an_address_once(client, tmp_path, monkeypatch):
    import interest
    monkeypatch.setattr(interest, "DATA_DIR", tmp_path)
    monkeypatch.setattr(interest, "_recent", {})
    r = client.post("/interest", data={"email": "a@b.com", "name": "Robin"})
    assert r.status_code == 200 and r.get_json()["status"] == "added"
    r = client.post("/interest", data={"email": "A@B.com", "name": "Robin"})
    assert r.get_json()["status"] == "already"
    assert interest.count() == 1


@pytest.mark.parametrize("bad", ["", "nope", "a@b", "a b@c.com", "@b.com", "a@@b.com"])
def test_interest_rejects_rubbish(client, tmp_path, monkeypatch, bad):
    import interest
    monkeypatch.setattr(interest, "DATA_DIR", tmp_path)
    monkeypatch.setattr(interest, "_recent", {})
    r = client.post("/interest", data={"email": bad})
    assert r.status_code == 400
    assert r.get_json()["status"] == "invalid"


def test_interest_is_rate_limited(client, tmp_path, monkeypatch):
    import interest
    monkeypatch.setattr(interest, "DATA_DIR", tmp_path)
    monkeypatch.setattr(interest, "_recent", {})
    codes = [client.post("/interest", data={"email": f"p{i}@q.com"}).status_code
             for i in range(interest.RATE_LIMIT + 2)]
    assert 429 in codes


def test_interest_export_needs_the_token(client, monkeypatch):
    monkeypatch.delenv("EXPORT_TOKEN", raising=False)
    assert client.get("/interest.csv").status_code == 404
    monkeypatch.setenv("EXPORT_TOKEN", "s3cret")
    assert client.get("/interest.csv").status_code == 404
    assert client.get("/interest.csv?token=wrong").status_code == 404
    r = client.get("/interest.csv?token=s3cret")
    assert r.status_code == 200
    assert r.data.startswith(b"created_at,email")


def test_an_address_can_be_removed(tmp_path, monkeypatch):
    import interest
    monkeypatch.setattr(interest, "DATA_DIR", tmp_path)
    interest.add("gone@example.com")
    assert interest.delete("GONE@example.com") is True
    assert interest.delete("gone@example.com") is False


def test_removal_endpoint_needs_the_token(client, tmp_path, monkeypatch):
    import interest
    monkeypatch.setattr(interest, "DATA_DIR", tmp_path)
    monkeypatch.setattr(interest, "_recent", {})
    interest.add("bye@example.com")
    monkeypatch.setenv("EXPORT_TOKEN", "s3cret")
    assert client.post("/interest/remove",
                       data={"email": "bye@example.com"}).status_code == 404
    r = client.post("/interest/remove?token=s3cret", data={"email": "bye@example.com"})
    assert r.status_code == 200 and r.get_json()["removed"] is True


def test_result_page_shows_the_helix_without_javascript(client):
    """The drawn helix must be in the HTML, not only after 3Dmol loads."""
    body = client.get("/p/Robin").data.decode()
    assert "viewer-still" in body
    assert "Idealised alpha helix" in body


def test_healthz_reports_a_loaded_proteome(client):
    data = client.get("/healthz").get_json()
    assert data["ok"] is True
    assert data["proteins"] > 20_000


def test_long_input_is_capped(client):
    r = client.post("/", data={"name": "A" * 500})
    assert r.status_code == 302
    assert len(r.headers["Location"].split("/p/")[1]) <= 80
