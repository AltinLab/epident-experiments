from MDAnalysis.exceptions import NoDataError
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis
import freesasa
import numpy as np


class PatchedSASAAnalysis(SASAAnalysis):

    def _prepare(self):
        self.results.total_area = np.zeros(
            self.n_frames,
            dtype=float,
        )
        self.results.residue_area = np.zeros(
            (self.n_frames, len(self.universe.residues.resids)),
            dtype=float,
        )
        self.results.relative_residue_area = np.zeros(
            (self.n_frames, len(self.universe.residues.resids)),
            dtype=float,
        )

    def _single_frame(self):
        """Calculate data from a single frame of trajectory"""

        structure = freesasa.Structure()
        # FreeSasa structure accepts PDBS if not available requires to reconstruct the structure using `addAtom`
        for a in self.atomgroup:
            x, y, z = a.position
            try:
                resname = a.resname
            except NoDataError:
                resname = "ANY"  # Default classifier value

            structure.addAtom(
                a.type.rjust(2), resname, a.resnum.item(), a.segid, x, y, z
            )

        # Define 1 cpu for windows avoid freesasa code to calculate it.
        parametes = freesasa.Parameters()
        if self._is_windows():
            parametes.setNThreads(1)

        result = freesasa.calc(structure, parametes)

        residue_areas = [
            result.residueAreas()[s][r]
            for s in list(result.residueAreas().keys())
            for r in list(result.residueAreas()[s].keys())
        ]
        self.results.total_area[self._frame_index] = result.totalArea()

        # Defend agains residue counts mismatch
        if len(self.universe.residues.resids) != len(residue_areas):
            raise ValueError(
                f"Residude count do not match the expectation, residue SASA not in results { len(self.universe.residues.resids)} != {len(residue_areas)}"
            )
        else:
            self.results.residue_area[self._frame_index] = [
                r.total for r in residue_areas
            ]
            self.results.relative_residue_area[self._frame_index] = [
                r.relativeTotal for r in residue_areas
            ]
