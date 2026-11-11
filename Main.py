
"""The code is written to produce the results of First-Exit-Aware Physical Trace Reduction
and Certified Mode Minimality
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CORE = ROOT / 'core'
sys.path.insert(0, str(CORE))

from source import canonical_verify_rq1_reachable as fixture
from source.verify_viable_full_disturbance import report as viable_exact
from source.verify_viable_family_interval import report as viable_family
from unified.exactbox import successor, residual
from unified.fixtures import fourteen, bridge, strict
from unified.efficient_pipeline import certify
from unified.global_minimum import exact_audit as affine_minimum_audit
from unified.pipeline import verify_pair
from unified.support import geometry as support_geometry, terminal_decisions, absorption
from unified.compare import benchmark as stage_benchmark
from unified.viable_unified import pipeline as polyhedral_viability
from unified.facetlp import solve as first_exit_facet_lp
from global_mimo.scientific_run import run as nonlinear_and_global_audit
from independent.viable_predicate_audit import verify as independent_predicate_verify


def write_csv(path: Path, rows: list[dict]):
    assert rows, f'Cannot write empty table: {path}'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def rational_facets():
    
    modes = fixture.build(6)
    reps = {m['grp']: m for m in modes if m['is_rep']}
    rows = []
    for m in modes:
        if m['is_rep']:
            continue
        j = reps[m['grp']]
        for t in range(4):
            for k in range(2):
                delta_a = [m['A'][k][c] - j['A'][k][c] for c in range(2)]
                delta_b = [m['B'][k][c] - j['B'][k][c] for c in range(2)]
                upper = (m['w'][k] - j['w'][k]
                         + sum(abs(delta_a[c]) * fixture.C[c] +
                               abs(delta_b[c]) * fixture.U[c] for c in range(2)))
                assert upper <= 0
                for sign in (-1, 1):
                    rows.append({'source': m['name'], 'representative': j['name'],
                                 'time': t, 'facet': f'{sign:+d}x{k+1}',
                                 'exact_guard_upper_bound_NOT_LP_optimum': str(upper)})
    assert len(rows) == 192
    return rows


def numerical_facets():
    
    models, pairs, X0, X, U, _ = fourteen()
    rows = []
    for src, rep in pairs:
        for result in first_exit_facet_lp(models[src], models[rep], X0, X, U, H=4):
            assert result['status'] in ('NUMERICAL_PASS', 'VACUOUS'), (src, rep, result)
            rows.append({'source': src, 'representative': rep,
                         'time': result['t'], 'facet': result['facet'],
                         'numerical_facet_margin': result.get('alpha', ''),
                         'solver_status': result['status']})
    assert len(rows) == 192
    return rows


def exact_and_numerical_evidence(out: Path, quick: bool = False):
    
    values: list[dict] = []

    def record(name, actual, expected=None, nature='exact rational'):
        if expected is not None:
            assert actual == expected, f'{name}: expected {expected!r}; got {actual!r}'
        values.append({'claim': name, 'computed_value': str(actual),
                       'verification': 'PASS', 'evidence_type': nature})
        print(f'  PASS  {name}: {actual}')

    print('\n[1/5] Rational fourteen-mode benchmark')
    models, pairs, X0, X, U, guard = fourteen()
    record('Original modes', len(models), 14)
    record('Candidate fixed representative pairs', len(pairs), 12)
    pair_rows = []
    for src, rep in pairs:
        verified = certify(models[src], models[rep], X0, X, U, 4, [guard], direct_guard=guard)
        assert verified['status'] == 'CERTIFIED_DIRECT'
        assert verified['lp_facet_calls'] == 0
        pair_rows.append({'source': src, 'representative': rep,
                          'certification': verified['status'],
                          'unnecessary_LPs_after_valid_direct_guard': 0,
                          'global_all_X_dominance': 'FAIL' if not residual(models[src], models[rep], X, U).within(models[rep].W) else 'PASS'})
    record('Pairwise fixed-representative certificates', len(pair_rows), 12)
    record('All-X global domination failures', sum(x['global_all_X_dominance'] == 'FAIL' for x in pair_rows), 12)
    record('Redundant LP checks with valid supplied direct guards', sum(r['unnecessary_LPs_after_valid_direct_guard'] for r in pair_rows), 0)
    rational_rows = rational_facets()
    record('Exact rational facet upper bounds', len(rational_rows), 192)

    print('\n[2/5] Exact viable strictness R1-T6')
    poly = polyhedral_viability(H=2, full_disturbance=True)
    assert poly['status'] == 'CERTIFIED' and poly['terminal_safe_nonempty']
    assert all(x['representative_facet_pass'] for x in poly['finite_stages'])
    assert viable_exact['outcome'] == 'PASS'
    record('R1-T6 exact finite facet tests', len(viable_exact['finite_first_exit_LP_facets']), 8)
    record('R1-T6 nonempty invariant continuation', poly['tail']['all_action_invariant'], True)
    record('Convex direct-initial-guard contradiction at 1/10000', viable_exact['strict_separation']['strict_gap'], '959/400000')
    record('Uniform parameter bound at 1/1000', viable_family['worst_positive_convex_guard_residual_gap'], '59/40000')
    record('Shared robust reach-and-continue', viable_exact['task']['all_future_states_safe_both_modes'], True)
    independently = independent_predicate_verify(F(1, 10000))
    record('Independent viable-predicate contradiction', independently['no_convex_direct_guard_gap'], '959/400000')
    write_csv(out / 'strictness_rational_values.csv', [
        {'epsilon': '1/10000', 'convex_guard_gap': viable_exact['strict_separation']['strict_gap'], 'method': 'R1-T6 Fraction witness'},
        {'epsilon': '1/1000', 'convex_guard_gap': viable_family['worst_positive_convex_guard_residual_gap'], 'method': 'R1-T6 monotone endpoint bound'}])

    print('\n[3/5] Global-MIMO component and mode-minimum certificates')
    affine = affine_minimum_audit()
    record('Nonempty original-mode subsets screened', affine['candidate_original_subsets_exhausted'], 16383)
    record('Original successor components', affine['exact_original_successor_components'], 2)
    record('Within-component intersection pairs', affine['within_group_intersecting_pairs'], 42)
    record('Cross-component disjoint pairs', affine['cross_group_disjoint_pairs'], 49)
    v = nonlinear_and_global_audit(ROOT)
    fourteen_global = v['original_14_mode_global_mimo_certificate']
    record('Global-MIMO necessary topology lower bound', fourteen_global['topological_lower_bound'], 2)
    record('Certified retained-mode upper bound', fourteen_global['independent_U_T1_upper_bound'], 2)
    record('Class-global exact physical-trace minimum', fourteen_global['minimum'], 2)
    record('Connected-outer one-step cost gap', v['sharp_connected_outer_cost_gap']['sharp_minimum_worst_case_cost_of_ANY_CONNECTED_OUTER_APPROXIMATION'], '7/100')
    record('Later-history initial components', v['history_conditioning']['initial_selected_probe_components'], 1)
    record('Later-history conditional components', v['history_conditioning']['history_conditioned_components'], 2)
    record('Nonlinear 3D/2-input minimum', v['coupled_nonlinear_3D_2input_global_minimum']['global_minimum_same_physical_state_connected_images'], 3)
    record('Single-task minimum (not universal trace minimum)', fourteen_global['one_prescribed_X_safety_task_minimum'], 1)

    print('  Additional finite-to-tail and restricted-scope regression cases')
    bi, bj, b0, bX, bU, bC = bridge()
    br = verify_pair(bi, bj, b0, bX, bU, 3, [bC], run_lp=False)
    record('Nonempty-terminal bridge case', br['status'], 'CERTIFIED')
    si, sj, s0, sX, sU = strict()
    vac = verify_pair(si, sj, s0, sX, sU, 3, [None], run_lp=False)
    record('Empty-safe-terminal negative-control status', vac['tail']['status'], 'CERTIFIED_EMPTY_TAIL')
    gg = support_geometry()
    record('Independent within-cluster graph edges', gg['within_cluster_intersections'], 42)
    td = terminal_decisions()
    record('Terminal decision intersection instances', td['checked_finite_set_intersection_instances'], 64)
    ab = absorption()
    record('Restricted absorbing-class horizon', ab['certified_horizon'], 2)
    return {'exact_viable': viable_exact, 'v': v, 'values': values}


def finite_state_comparators(out):
    from task_specific.methods import (
        unreachable_component_control, v_interval_gridded_fixture,
        grid_mimo_fixture, krook_rsbs_partition, reachable_simulation_2025,
        semi_complete_2026)
    rows = []
    for bench in ('v_gridded_source', 'v_gridded_representative',
                  'coupled_robust', 'coupled_sparse', 'coupled_recurrent',
                  'unreachable_component_control'):
        if bench.startswith('v_gridded_'):
            ts = v_interval_gridded_fixture(bench.removeprefix('v_gridded_'))
        elif bench == 'unreachable_component_control':
            ts = unreachable_component_control()
        else:
            ts = grid_mimo_fixture(bench)
        for method, func in (
            ('Krook-2024', lambda: krook_rsbs_partition(ts)),
            ('Ganty-2025', lambda: reachable_simulation_2025(ts)),
            ('Ganty-2026-reachable-region ', lambda: semi_complete_2026(ts, ts.reachable())),
            ('Ganty-2026-full-region control', lambda: semi_complete_2026(ts, None))):
            r = func()
            rows.append({'benchmark': bench, 'adapted_method': method,
                         'finite_states': len(ts.states), 'reachable_states': len(ts.reachable()),
                         'output_blocks': r.get('blocks', ''),
                         'output_simulation_pairs': r.get('simulating_pairs', ''),
                         'iterations': r.get('iterations', '')})
    assert len(rows) == 24
    write_csv(out / 'finite_state_adaptations.csv', rows)
    print('  PASS  24 finite-state adaptation outcomes across 6 benchmarks')


def fresh_timings(out: Path):
    r = stage_benchmark(repetitions=3)
    measures = ['global_check_median_ms', 'direct_guard_median_ms',
                'unified_fraction_outer_plus_tail_median_ms',
                'unified_finite_LP_plus_tail_median_ms']
    write_csv(out / 'fresh_offline_timings.csv',
              [{'stage': name, 'new_machine_median_ms': r[name],
                'repetitions': r['N_repetitions']}
               for name in measures])
    for k in measures:
        print(f'        {k}: {r[k]:.3f} ms')


def run_tests():
    suite = unittest.defaultTestLoader.discover(start_dir=str(ROOT / 'tests'), pattern='test_*.py')
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    if not result.wasSuccessful():
        raise RuntimeError(f'Unit tests failed: {len(result.failures)} failures / {len(result.errors)} errors')
    print(f'  PASS  {result.testsRun} source-based independent tests')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=ROOT / 'results', help='Generated outputs directory')
    args = p.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    begin = time.perf_counter()
    evidence = exact_and_numerical_evidence(out, quick=args.quick)
    print('\n[4/5] Algorithm-contract tests and contextual finite-state comparisons')
    run_tests()
    if not args.skip_comparators and not args.quick:
        finite_state_comparators(out)
    if not args.skip_timing and not args.quick:
        fresh_timings(out)
    if not args.skip_figures:
        from plot_paper import plot_all
        plot_all(evidence, out)
    print(f'Completed in {time.perf_counter()-begin:.1f} seconds. Outputs: {out}')


if __name__ == '__main__':
    main()
