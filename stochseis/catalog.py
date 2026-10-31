"""Catalog container and file format.

Events are stored in a numpy structured array (one row per event); slip
distributions of large events are stored in compressed-row form
(``slip_ptr``, ``slip_idx``, ``slip_val``).  Everything, including the
fault geometry and the parameter set, goes into one ``.npz`` file.
"""

from __future__ import annotations

import json

import numpy as np

from .geometry import Fault
from .params import Params

EVENT_DTYPE = np.dtype([
    ("t", float),          # origin time [yr]
    ("mag", float),        # moment magnitude
    ("m0", float),         # seismic moment [N m]
    ("kind", np.int8),     # 1 = fault event (grown rupture), 0 = sub-patch GR event
    ("nuc", int),          # nucleation patch index
    ("cs", float),         # centroid along strike [m]
    ("czd", float),        # centroid down dip [m]
    ("cx", float), ("cy", float), ("cz", float),   # Cartesian centroid [m]
    ("length", float),     # along-strike extent [m]
    ("width", float),      # down-dip extent [m]
    ("area", float),       # ruptured area [m^2]
    ("npatch", int),
    ("dmean", float),      # mean slip [m]
    ("dmax", float),       # maximum slip [m]
    ("dsig", float),       # target (= realised moment-weighted) stress drop [Pa]
    ("duration", float),   # rupture duration [s] (max onset time over the rupture)
    ("vr", float),         # rupture speed [m/s]
    ("surface", np.int8),  # 1 if surface breaking
    ("stop", np.int8),     # 0 hazard, 1 deficit, 2 boundary, -1 n/a
    ("rho", float),        # deficit ratio at arrest
    ("parent", int),       # index of triggering event, -1 for background
    ("gen", np.int16),     # generation in the ETAS cascade
    ("has_slip", np.int8), # 1 if the slip distribution is stored
])


def moment_magnitude(m0):
    return 2.0 / 3.0 * (np.log10(m0) - 9.05)


def moment_from_magnitude(mw):
    return 10.0 ** (1.5 * mw + 9.05)


class Catalog:
    def __init__(self, events, fault: Fault, params: Params, meta: dict,
                 slip_ptr, slip_idx, slip_val, slip_event, released, deficit0):
        self.events = events
        self.fault = fault
        self.params = params
        self.meta = meta
        self.slip_ptr = np.asarray(slip_ptr, int)
        self.slip_idx = np.asarray(slip_idx, int)
        self.slip_val = np.asarray(slip_val, float)
        self.slip_event = np.asarray(slip_event, int)   # event index of each stored slip map
        self.released = np.asarray(released, float)     # cumulative slip per patch at end [m]
        self.deficit0 = np.asarray(deficit0, float)     # initial deficit per patch [m]
        self._slip_lookup = {int(e): k for k, e in enumerate(self.slip_event)}

    # ------------------------------------------------------------------
    def __len__(self):
        return len(self.events)

    @property
    def duration(self) -> float:
        return float(self.params.duration_yr)

    def slip_of(self, event_index):
        """(patch indices, slip [m]) of a stored event, or None."""
        k = self._slip_lookup.get(int(event_index))
        if k is None:
            return None
        a, b = self.slip_ptr[k], self.slip_ptr[k + 1]
        return self.slip_idx[a:b], self.slip_val[a:b]

    def fault_events(self):
        return self.events[self.events["kind"] == 1]

    def final_deficit(self):
        """Slip deficit per patch at the end of the run [m]."""
        return self.fault.slip_rate * self.duration + self.deficit0 - self.released

    # ------------------------------------------------------------------
    def save(self, path: str) -> None:
        d = dict(events=self.events, slip_ptr=self.slip_ptr, slip_idx=self.slip_idx,
                 slip_val=self.slip_val, slip_event=self.slip_event,
                 released=self.released, deficit0=self.deficit0,
                 params_json=np.array(json.dumps(self.params.as_dict())),
                 meta_json=np.array(json.dumps(self.meta, default=float)))
        for k, v in self.fault.to_dict().items():
            d["fault_" + k] = v
        np.savez_compressed(path, **d)

    @classmethod
    def load(cls, path: str) -> "Catalog":
        z = np.load(path, allow_pickle=False)
        params = Params.from_dict(json.loads(str(z["params_json"])))
        meta = json.loads(str(z["meta_json"]))
        fd = {k[6:]: z[k] for k in z.files if k.startswith("fault_")}
        fault = Fault.from_dict(fd)
        return cls(z["events"], fault, params, meta, z["slip_ptr"], z["slip_idx"],
                   z["slip_val"], z["slip_event"], z["released"], z["deficit0"])
