import numpy as np

import femlabpy as fp
from femlabpy.examples import (
    run_bar01_nlbar,
    run_cantilever,
    run_square_plastpe,
    run_square_plastps,
)


def test_matlab_data_aliases_expose_original_inputs():
    assert fp.canti()["T"].shape == (16, 5)
    assert fp.flow()["T1"].shape[1] == 5
    assert "T" in fp.bar01()
    assert "T" in fp.bar02()
    assert "T" in fp.bar03()
    assert "T" in fp.square()
    assert "T" in fp.hole(plane_strain=True)


def test_elastic_alias_matches_cantilever_runner():
    data = fp.canti()
    result = fp.elastic(
        data["T"], data["X"], data["G"], data["C"], data["P"], dof=int(data["dof"])
    )
    reference = run_cantilever(plot=False)
    np.testing.assert_allclose(result["u"], reference["u"])
    np.testing.assert_allclose(result["S"], reference["S"])


def test_elastic_auto_detects_t3_and_reduces_planar_gmsh_coordinates():
    X = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    T = np.array([[1, 2, 3, 1]], dtype=int)
    G = np.array([[210000.0, 0.3, 1.0]], dtype=float)
    C = np.array([[1, 1, 0.0], [1, 2, 0.0], [2, 2, 0.0]], dtype=float)
    P = np.array([[2, 1, 100.0]], dtype=float)

    result = fp.elastic(T, X, G, C, P, dof=2)

    assert result["u"].shape == (6, 1)
    assert result["S"].shape == (1, 3)
    assert result["E"].shape == (1, 3)
    assert result["data"]["X"].shape == (3, 2)
    assert np.isfinite(result["u"]).all()


def test_elastic_accepts_explicit_t3_etype():
    X = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    T = np.array([[1, 2, 3, 1]], dtype=int)
    G = np.array([[210000.0, 0.3, 1.0]], dtype=float)
    C = np.array([[1, 1, 0.0], [1, 2, 0.0], [2, 2, 0.0]], dtype=float)
    P = np.array([[2, 1, 100.0]], dtype=float)

    result = fp.elastic(T, X, G, C, P, dof=2, etype="t3")

    assert result["u"].shape == (6, 1)


def test_flow_aliases_execute():
    q4 = fp.flowq4(plot=False)
    t3 = fp.flowt3(plot=False)
    assert q4["u"].shape[0] == q4["data"]["X"].shape[0]
    assert t3["u"].shape[0] == t3["data"]["X"].shape[0]
    flow_data = fp.flow()
    q4_generic = fp.flowq4(
        flow_data["T1"],
        flow_data["X"],
        flow_data["G"],
        flow_data["C"],
        dof=int(flow_data["dof"]),
        plot=False,
    )
    t3_generic = fp.flowt3(
        flow_data["T2"],
        flow_data["X"],
        flow_data["G"],
        flow_data["C"],
        dof=int(flow_data["dof"]),
        plot=False,
    )
    np.testing.assert_allclose(q4_generic["u"], q4["u"])
    np.testing.assert_allclose(t3_generic["u"], t3["u"])


def test_nlbar_alias_matches_legacy_runner():
    data = fp.bar01()
    result = fp.nlbar(
        data["T"],
        data["X"],
        data["G"],
        data["C"],
        data["P"],
        no_loadsteps=int(data["no_loadsteps"][0, 0]),
        i_max=int(data["i_max"][0, 0]),
        i_d=int(data["i_d"][0, 0]),
        plotdof=int(data["plotdof"][0, 0]),
        tol=float(data["TOL"][0, 0]),
    )
    reference = run_bar01_nlbar(plot=False)
    np.testing.assert_allclose(result["u"], reference["u"])
    np.testing.assert_allclose(result["F_path"], reference["F_path"])


def test_plastic_aliases_match_legacy_runners():
    square_ps = fp.square(plane_strain=False)
    result_ps = fp.plastps(
        square_ps["T"],
        square_ps["X"],
        square_ps["G"],
        square_ps["C"],
        square_ps["P"],
        no_loadsteps=int(square_ps["no_loadsteps"][0, 0]),
        i_max=int(square_ps["i_max"][0, 0]),
        i_d=int(square_ps["i_d"][0, 0]),
        plotdof=int(square_ps["plotdof"][0, 0]),
        tol=float(square_ps["TOL"][0, 0]),
    )
    reference_ps = run_square_plastps(plot=False)
    np.testing.assert_allclose(result_ps["u"], reference_ps["u"])
    np.testing.assert_allclose(
        np.asarray(result_ps["E"]).reshape(-1),
        np.asarray(reference_ps["E"]).reshape(-1),
    )

    square_pe = fp.square(plane_strain=True)
    result_pe = fp.plastpe(
        square_pe["T"],
        square_pe["X"],
        square_pe["G"],
        square_pe["C"],
        square_pe["P"],
        no_loadsteps=int(square_pe["no_loadsteps"][0, 0]),
        i_max=int(square_pe["i_max"][0, 0]),
        i_d=int(square_pe["i_d"][0, 0]),
        plotdof=int(square_pe["plotdof"][0, 0]),
        tol=float(square_pe["TOL"][0, 0]),
    )
    reference_pe = run_square_plastpe(plot=False)
    np.testing.assert_allclose(result_pe["u"], reference_pe["u"])
    np.testing.assert_allclose(
        np.asarray(result_pe["E"]).reshape(-1),
        np.asarray(reference_pe["E"]).reshape(-1),
    )


def test_matlab_api_docstrings_present():
    for name in (
        "setpath",
        "load_gmsh",
        "load_gmsh2",
        "canti",
        "flow",
        "bar01",
        "bar02",
        "bar03",
        "square",
        "hole",
        "elastic",
        "flowq4",
        "flowt3",
        "nlbar",
        "plastps",
        "plastpe",
        "solve_nlbar",
        "solve_plastic",
    ):
        doc = getattr(fp, name).__doc__
        assert doc is not None
        assert len(doc.split()) >= 8
