# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""stochseis: a stochastic seismicity generator on fault patch graphs.

See README.md and docs/ for the model description.  Typical use::

    from stochseis import Params, run
    cat = run(Params(duration_yr=1000, seed=3))
    cat.save("cat.npz")

followed by ``python -m stochseis.cli analyze cat.npz --plots out/``.
"""

from .catalog import Catalog
from .params import Params
from .simulator import Simulation, run

__all__ = ["Params", "Catalog", "Simulation", "run"]
__version__ = "0.1.0"
