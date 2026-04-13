import matplotlib
import numpy as np

matplotlib.use("Agg")

import femlabpy as fp
from femlabpy.examples import run_bar01_nlbar


def test_readme_nonlinear_truss_example_runs():
    result = run_bar01_nlbar(plot=True)

    assert result["U_path"].shape[0] > 1
    assert result["F_path"].shape == result["U_path"].shape
    assert len(result["figures"]) == 2


def test_readme_low_level_q4_assembly_example_runs():
    nn = 4
    dof = 2
    K, p, q = fp.init(nn, dof)

    T = np.array([[1, 2, 3, 4, 1]], dtype=int)
    X = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=float)
    G = np.array([[210000.0, 0.3, 1.0, 1.0]], dtype=float)

    K = fp.kq4e(K, T, X, G)

    C = np.array([[1, 1, 0.0], [1, 2, 0.0], [4, 1, 0.0]], dtype=float)
    P = np.array([[2, 1, 1000.0], [3, 1, 1000.0]], dtype=float)

    p = fp.setload(p, P)
    K, p, bcwt = fp.setbc(K, p, C, dof)
    u = np.linalg.solve(K, p)

    assert q.shape == p.shape
    assert bcwt > 0.0
    assert np.isfinite(u).all()
