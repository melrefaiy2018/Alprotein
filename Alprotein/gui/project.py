"""
.alproj — Alprotein project file.

A project bundles everything needed to reopen a session:

* the input PDB (embedded verbatim, so the project is portable);
* the calculation parameters (Tools panel + axis range);
* user-applied site-energy overrides;
* any computed numerical arrays (Hamiltonian, eigenvalues, eigenvectors,
  absorption / fluorescence spectra, exciton distributions);
* a free-text notebook;
* basic metadata (version, timestamps).

On disk it is a zip archive containing:

    project.json        # metadata + parameters + overrides + notebook
    structure.pdb       # raw PDB text
    arrays.h5           # numerical arrays (created lazily — only if any
                          calculation has been run)

The format is deliberately simple so it can be inspected with standard
tooling (``unzip -l``, ``h5dump``) and so we can evolve the schema without
breaking forward compatibility — readers ignore unknown keys.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

logger = logging.getLogger(__name__)


PROJECT_VERSION = 1
PROJECT_SUFFIX = ".alproj"
_METADATA_NAME = "project.json"
_PDB_NAME = "structure.pdb"
_ARRAYS_NAME = "arrays.h5"


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class Project:
    """In-memory representation of a project."""

    # Metadata --------------------------------------------------------------
    name: str = ""
    version: int = PROJECT_VERSION
    created: str = field(default_factory=_utcnow)
    last_modified: str = field(default_factory=_utcnow)

    # Source ----------------------------------------------------------------
    source_pdb_path: str = ""  # absolute path the project was first loaded from
    pdb_text: str = ""         # raw PDB contents, embedded for portability

    # Parameters & user edits ----------------------------------------------
    parameters: Dict[str, Any] = field(default_factory=dict)
    site_energy_overrides: Dict[str, float] = field(default_factory=dict)
    vacuum_energies: Dict[str, float] = field(default_factory=dict)

    # Computed arrays — stored as numpy or None ----------------------------
    hamiltonian: Optional[np.ndarray] = None
    eigenvalues: Optional[np.ndarray] = None
    eigenvectors: Optional[np.ndarray] = None
    site_energies: Dict[str, float] = field(default_factory=dict)
    calculated_site_energies: Dict[str, float] = field(default_factory=dict)
    pigment_labels: Optional[list] = None
    spectrum_data: Optional[Dict[str, np.ndarray]] = None
    exciton_distributions: Optional[Dict[str, np.ndarray]] = None
    exciton_labels: Optional[list] = None

    # Misc ------------------------------------------------------------------
    notebook: str = ""
    view_state: Dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Disk I/O
    # ------------------------------------------------------------------

    def save(self, path: str | Path) -> Path:
        """Write the project to ``path``. Returns the absolute path written."""
        target = Path(path).expanduser().resolve()
        if target.suffix != PROJECT_SUFFIX:
            target = target.with_suffix(PROJECT_SUFFIX)
        target.parent.mkdir(parents=True, exist_ok=True)

        self.last_modified = _utcnow()

        meta = self._metadata_dict()
        with ZipFile(target, "w", compression=ZIP_DEFLATED) as zf:
            zf.writestr(_METADATA_NAME, json.dumps(meta, indent=2))
            if self.pdb_text:
                zf.writestr(_PDB_NAME, self.pdb_text)
            if self._has_arrays():
                arrays_bytes = self._dump_arrays_h5()
                zf.writestr(_ARRAYS_NAME, arrays_bytes)

        logger.info("Project saved to %s", target)
        return target

    @classmethod
    def load(cls, path: str | Path) -> "Project":
        """Load a project from ``path``."""
        source = Path(path).expanduser().resolve()
        if not source.exists():
            raise FileNotFoundError(source)

        with ZipFile(source, "r") as zf:
            names = set(zf.namelist())
            if _METADATA_NAME not in names:
                raise ValueError(
                    f"{source} is not a valid {PROJECT_SUFFIX} file: missing {_METADATA_NAME}"
                )
            meta = json.loads(zf.read(_METADATA_NAME).decode("utf-8"))

            pdb_text = ""
            if _PDB_NAME in names:
                pdb_text = zf.read(_PDB_NAME).decode("utf-8")

            arrays: Dict[str, Any] = {}
            if _ARRAYS_NAME in names:
                arrays = cls._load_arrays_h5(zf.read(_ARRAYS_NAME))

        # Filter unknown metadata keys forward-compatibly.
        proj = cls._from_metadata(meta, pdb_text=pdb_text, arrays=arrays)
        logger.info("Project loaded from %s", source)
        return proj

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _metadata_dict(self) -> Dict[str, Any]:
        return {
            "schema": "alprotein.project",
            "schema_version": PROJECT_VERSION,
            "name": self.name,
            "created": self.created,
            "last_modified": self.last_modified,
            "source_pdb_path": self.source_pdb_path,
            "parameters": self.parameters,
            "site_energy_overrides": self.site_energy_overrides,
            "vacuum_energies": self.vacuum_energies,
            "site_energies": self.site_energies,
            "calculated_site_energies": self.calculated_site_energies,
            "pigment_labels": list(self.pigment_labels) if self.pigment_labels else None,
            "exciton_labels": list(self.exciton_labels) if self.exciton_labels else None,
            "notebook": self.notebook,
            "view_state": self.view_state,
        }

    @classmethod
    def _from_metadata(
        cls,
        meta: Mapping[str, Any],
        *,
        pdb_text: str,
        arrays: Mapping[str, Any],
    ) -> "Project":
        proj = cls(
            name=str(meta.get("name", "")),
            version=int(meta.get("schema_version", PROJECT_VERSION)),
            created=str(meta.get("created", _utcnow())),
            last_modified=str(meta.get("last_modified", _utcnow())),
            source_pdb_path=str(meta.get("source_pdb_path", "")),
            pdb_text=pdb_text,
            parameters=dict(meta.get("parameters") or {}),
            site_energy_overrides=dict(meta.get("site_energy_overrides") or {}),
            vacuum_energies=dict(meta.get("vacuum_energies") or {}),
            site_energies=dict(meta.get("site_energies") or {}),
            calculated_site_energies=dict(meta.get("calculated_site_energies") or {}),
            pigment_labels=list(meta.get("pigment_labels") or []) or None,
            exciton_labels=list(meta.get("exciton_labels") or []) or None,
            notebook=str(meta.get("notebook", "")),
            view_state=dict(meta.get("view_state") or {}),
        )
        proj.hamiltonian = arrays.get("hamiltonian")
        proj.eigenvalues = arrays.get("eigenvalues")
        proj.eigenvectors = arrays.get("eigenvectors")
        proj.spectrum_data = arrays.get("spectrum_data")
        proj.exciton_distributions = arrays.get("exciton_distributions")
        return proj

    def _has_arrays(self) -> bool:
        return any(
            x is not None
            for x in (
                self.hamiltonian,
                self.eigenvalues,
                self.eigenvectors,
                self.spectrum_data,
                self.exciton_distributions,
            )
        )

    def _dump_arrays_h5(self) -> bytes:
        """Serialize numerical arrays to an in-memory HDF5 blob."""
        import io

        import h5py

        buf = io.BytesIO()
        with h5py.File(buf, "w") as f:
            if self.hamiltonian is not None:
                f.create_dataset("hamiltonian", data=np.asarray(self.hamiltonian))
            if self.eigenvalues is not None:
                f.create_dataset("eigenvalues", data=np.asarray(self.eigenvalues))
            if self.eigenvectors is not None:
                f.create_dataset("eigenvectors", data=np.asarray(self.eigenvectors))
            if self.spectrum_data:
                grp = f.create_group("spectrum_data")
                for key, val in self.spectrum_data.items():
                    if val is None:
                        continue
                    grp.create_dataset(key, data=np.asarray(val))
            if self.exciton_distributions:
                grp = f.create_group("exciton_distributions")
                for key, arr in self.exciton_distributions.items():
                    if arr is None:
                        continue
                    grp.create_dataset(str(key), data=np.asarray(arr))
        return buf.getvalue()

    @staticmethod
    def _load_arrays_h5(blob: bytes) -> Dict[str, Any]:
        import io

        import h5py

        out: Dict[str, Any] = {}
        with h5py.File(io.BytesIO(blob), "r") as f:
            for key in ("hamiltonian", "eigenvalues", "eigenvectors"):
                if key in f:
                    out[key] = np.array(f[key][...])
            if "spectrum_data" in f:
                grp = f["spectrum_data"]
                out["spectrum_data"] = {k: np.array(grp[k][...]) for k in grp.keys()}
            if "exciton_distributions" in f:
                grp = f["exciton_distributions"]
                out["exciton_distributions"] = {
                    k: np.array(grp[k][...]) for k in grp.keys()
                }
        return out


__all__ = ["Project", "PROJECT_VERSION", "PROJECT_SUFFIX"]
