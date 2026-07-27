"""A list of people who would buy a print.

There is no shop. The point of this is to find out whether building one is
worth it before building it, so all it records is an email address, the peptide
the person was looking at when they asked, and when.

SQLite on a mounted volume rather than a Postgres service: this will hold tens
of rows, a second Railway service would cost more per month than the answer is
worth, and two gunicorn workers against one WAL-mode SQLite file is fine. If
DATA_DIR is not a real volume the data will not survive a redeploy, so the
health endpoint reports which path is in use.
"""

from __future__ import annotations

import os
import re
import sqlite3
import threading
import time
from pathlib import Path

DATA_DIR = Path(os.environ.get("DATA_DIR", "/data"))
FALLBACK_DIR = Path(__file__).resolve().parent / "var"

# Deliberately permissive. The job is to reject obvious rubbish, not to
# adjudicate RFC 5321, and a real typo is better caught by nothing arriving.
EMAIL = re.compile(r"^[^@\s,;]{1,64}@[^@\s,;.]+(\.[^@\s,;.]+)+$")
MAX_EMAIL = 254

_lock = threading.Lock()
_recent: dict[str, list[float]] = {}
RATE_LIMIT = 5           # submissions
RATE_WINDOW = 300.0      # per five minutes, per address


def _dir() -> Path:
    if DATA_DIR.is_dir() and os.access(DATA_DIR, os.W_OK):
        return DATA_DIR
    FALLBACK_DIR.mkdir(parents=True, exist_ok=True)
    return FALLBACK_DIR


def db_path() -> Path:
    return _dir() / "interest.sqlite3"


def durable() -> bool:
    """False when we fell back to container-local storage, which is ephemeral."""
    return _dir() == DATA_DIR


def _connect() -> sqlite3.Connection:
    conn = sqlite3.connect(db_path(), timeout=5.0)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    return conn


def init() -> None:
    with _connect() as conn:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS signups (
                id INTEGER PRIMARY KEY,
                email TEXT NOT NULL UNIQUE,
                name TEXT,
                sequence TEXT,
                wants TEXT,
                created_at TEXT NOT NULL DEFAULT (datetime('now'))
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS signups_created ON signups(created_at)")


def valid(email: str) -> bool:
    return bool(email) and len(email) <= MAX_EMAIL and bool(EMAIL.match(email))


def throttled(key: str) -> bool:
    now = time.time()
    with _lock:
        hits = [t for t in _recent.get(key, []) if now - t < RATE_WINDOW]
        if len(hits) >= RATE_LIMIT:
            _recent[key] = hits
            return True
        hits.append(now)
        _recent[key] = hits
        # Keep the table from growing without bound on a long-lived process.
        if len(_recent) > 5000:
            for k in [k for k, v in _recent.items()
                      if not v or now - v[-1] > RATE_WINDOW]:
                _recent.pop(k, None)
    return False


def add(email: str, name: str = "", sequence: str = "", wants: str = "print") -> str:
    """Record interest. Returns 'added', 'already' or 'invalid'."""
    email = (email or "").strip().lower()
    if not valid(email):
        return "invalid"
    init()
    with _connect() as conn:
        existing = conn.execute(
            "SELECT 1 FROM signups WHERE email = ?", (email,)).fetchone()
        if existing:
            return "already"
        conn.execute(
            "INSERT INTO signups (email, name, sequence, wants) VALUES (?, ?, ?, ?)",
            (email, (name or "")[:80], (sequence or "")[:80], (wants or "")[:40]),
        )
    return "added"


def count() -> int:
    try:
        init()
        with _connect() as conn:
            return conn.execute("SELECT COUNT(*) FROM signups").fetchone()[0]
    except sqlite3.Error:
        return 0


def rows() -> list[tuple]:
    init()
    with _connect() as conn:
        return conn.execute(
            "SELECT created_at, email, name, sequence, wants "
            "FROM signups ORDER BY id DESC").fetchall()


def delete(email: str) -> bool:
    """So somebody can be taken off the list when they ask."""
    init()
    with _connect() as conn:
        cur = conn.execute("DELETE FROM signups WHERE email = ?",
                           ((email or "").strip().lower(),))
        return cur.rowcount > 0
