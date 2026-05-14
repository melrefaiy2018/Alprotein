"""
Recent-files list persisted at ``~/.alprotein/recent.json``.

Tracks both PDB files and ``.alproj`` projects opened via the workbench.
The list is most-recent-first and capped at :data:`MAX_RECENT`.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

MAX_RECENT = 10
_RECENT_DIR = Path.home() / ".alprotein"
_RECENT_FILE = _RECENT_DIR / "recent.json"


@dataclass(frozen=True)
class RecentEntry:
    """One recent-files entry."""

    path: str
    kind: str  # "pdb" | "project"
    opened_at: str  # ISO-8601 UTC


class RecentFiles:
    """Persistent, capped MRU list of opened files."""

    def __init__(self, storage_path: Optional[Path] = None) -> None:
        self._path = Path(storage_path) if storage_path else _RECENT_FILE
        self._entries: List[RecentEntry] = []
        self._load()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def add(self, path: str | Path, kind: str = "pdb") -> None:
        """Record an opened file; deduplicates by path."""
        p = str(Path(path).expanduser().resolve())
        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        # Remove any existing entry with the same path (we'll re-add at the top).
        self._entries = [e for e in self._entries if e.path != p]
        self._entries.insert(0, RecentEntry(path=p, kind=kind, opened_at=now))
        del self._entries[MAX_RECENT:]
        self._save()

    def entries(self) -> List[RecentEntry]:
        """Return existing entries in MRU order; missing files are pruned lazily."""
        live = [e for e in self._entries if Path(e.path).exists()]
        if len(live) != len(self._entries):
            self._entries = live
            self._save()
        return list(self._entries)

    def clear(self) -> None:
        self._entries = []
        self._save()

    # ------------------------------------------------------------------
    # Storage
    # ------------------------------------------------------------------

    def _load(self) -> None:
        if not self._path.exists():
            return
        try:
            raw = json.loads(self._path.read_text("utf-8"))
        except (OSError, json.JSONDecodeError) as e:
            logger.warning("Could not read %s: %s — starting empty", self._path, e)
            return
        items = raw.get("entries", []) if isinstance(raw, dict) else raw
        loaded: List[RecentEntry] = []
        for item in items:
            if not isinstance(item, dict):
                continue
            path = item.get("path")
            if not isinstance(path, str):
                continue
            loaded.append(
                RecentEntry(
                    path=path,
                    kind=str(item.get("kind", "pdb")),
                    opened_at=str(item.get("opened_at", "")),
                )
            )
        self._entries = loaded[:MAX_RECENT]

    def _save(self) -> None:
        try:
            self._path.parent.mkdir(parents=True, exist_ok=True)
            payload = {"entries": [asdict(e) for e in self._entries]}
            self._path.write_text(json.dumps(payload, indent=2), "utf-8")
        except OSError as e:
            logger.warning("Could not write %s: %s", self._path, e)


__all__ = ["RecentFiles", "RecentEntry", "MAX_RECENT"]
