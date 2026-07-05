"""Smoke test: build a tiny flat SYNTHETIC stack model through the wheel and
read positive volumes off it. Synthetic data only — no real dataset content."""

import json

import petekstatic


def test_flat_model_yields_positive_volumes():
    # A hand-authored flat block: 11x11 top at 2000 m, OWC at 2100 m, 1 km^2
    # footprint, 50 m gross, 5 layers, textbook priors.
    m = petekstatic.build_flat_model(
        n=11,
        depth_m=2000.0,
        owc_m=2100.0,
        area_m2=1_000_000.0,
        gross_height_m=50.0,
        nk=5,
        porosity=0.25,
        net_to_gross=0.8,
        water_saturation=0.3,
    )

    assert m.bulk_volume() > 0.0

    ip = m.in_place(boi=1.25)
    assert ip["grv_m3"] > 0.0
    assert ip["hcpv_m3"] > 0.0
    assert ip["ooip_sm3"] > 0.0
    # OOIP is HCPV divided by Boi (>1), so it must be the smaller number.
    assert ip["ooip_sm3"] < ip["hcpv_m3"]


def test_bundles_and_zone_rollup_serialize():
    m = petekstatic.build_flat_model()

    # Per-zone rollup is valid JSON with a positive total GRV.
    zoned = json.loads(m.in_place_by_zone())
    assert zoned["total"]["grv_m3"] > 0.0
    assert isinstance(zoned["zones"], list)

    # The map + volume bundles serialize to non-empty JSON strings.
    map_json = json.loads(m.map_bundle())
    assert isinstance(map_json, dict)

    section_json = json.loads(
        m.intersection_bundle([[0.1, 0.5], [0.9, 0.5]])
    )
    assert isinstance(section_json, dict)
