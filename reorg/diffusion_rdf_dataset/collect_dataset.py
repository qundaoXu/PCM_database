#!/usr/bin/env python3
"""Collect, rebuild and verify the approved local RDF first-peak dataset.

Python >=3.10; dependencies: numpy, scipy. All candidate membership and source
hashes are frozen in collection_manifest.json before scoring. Original extxyz
frame bytes are copied verbatim; scientific coordinates are never rewritten.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import io
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import time
import zlib

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.signal import find_peaks
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
BIN = 0.05
EDGES = np.arange(71, dtype=float) * BIN
RADII = (EDGES[:-1] + EDGES[1:]) / 2
SHELLS = 4 * np.pi / 3 * np.diff(EDGES ** 3)
WINDOW = (RADII >= 2.0) & (RADII <= 3.3)
SYMBOLS = {'32': 'Ge', '34': 'Se', '51': 'Sb', '52': 'Te'}


def log(s):
    print(time.strftime('[%H:%M:%S] ') + str(s), flush=True)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def comp_key(counts):
    counts = {k: int(v) for k, v in counts.items() if v}
    gcd = math.gcd(*counts.values())
    return ';'.join(f'{k}:{v // gcd}' for k, v in sorted(counts.items()))


def comp_name(key):
    return ''.join(s + (v if int(v) != 1 else '') for s, v in
                   (p.split(':') for p in key.split(';')))


def geometry_id(symbols, pos, cell):
    xyz = np.rint(pos * 1e6).astype('<i8')
    labels = np.array(symbols, dtype='U3')
    order = np.lexsort((xyz[:, 2], xyz[:, 1], xyz[:, 0], labels))
    h = hashlib.sha256(np.rint(cell.reshape(-1) * 1e6).astype('<i8').tobytes())
    h.update(labels[order].astype('S3').tobytes())
    h.update(xyz[order].tobytes())
    return h.hexdigest()


def xyz_frames(f):
    """Yield raw frame bytes and parsed physical arrays, with no coordinate edit."""
    index = 0
    while line := f.readline():
        if not line.strip():
            continue
        n = int(line)
        header = f.readline()
        atoms = [f.readline() for _ in range(n)]
        if not header or any(not b for b in atoms):
            raise ValueError(f'Truncated extxyz frame {index}')
        lattice = re.search(rb'Lattice="([^"]+)"', header)
        if lattice is None:
            raise ValueError(f'Missing lattice in frame {index}')
        cell = np.array(lattice[1].split(), dtype=float).reshape(3, 3)
        fields = re.search(rb'Properties=([^\s]+)', header)
        columns, offset = {}, 0
        if fields:
            fields = fields[1].decode().split(':')
            for k in range(0, len(fields), 3):
                columns[fields[k]] = offset
                offset += int(fields[k + 2])
        sc, pc = columns.get('species', columns.get('Z', 0)), columns.get('pos', 1)
        split = [b.split() for b in atoms]
        species = [SYMBOLS.get(v[sc].decode(), v[sc].decode()) for v in split]
        positions = np.array([v[pc:pc + 3] for v in split], dtype=float)
        if positions.shape != (n, 3) or not np.isfinite(positions).all():
            raise ValueError(f'Invalid positions in frame {index}')
        yield index, line + header + b''.join(atoms), species, positions, cell
        index += 1


def rdf(pos, cell):
    """Finite-N-corrected total number RDF; periodic pair list, not an N*N matrix."""
    if not np.isfinite(cell).all() or np.linalg.det(cell) <= 0:
        raise ValueError('Invalid cell')
    if not np.allclose(cell, np.diag(np.diag(cell)), atol=1e-9, rtol=0):
        # Current approved inputs are orthorhombic. Refuse a silent wrong MIC.
        raise ValueError('Non-orthorhombic cell requires an explicit MIC implementation')
    box = np.diag(cell)
    if np.min(box) <= 2 * EDGES[-1]:
        raise ValueError('Cell too small for the fixed 3.5 A RDF radius')
    p = pos % box
    pairs = cKDTree(p, boxsize=box).query_pairs(EDGES[-1], output_type='ndarray')
    if len(pairs):
        delta = p[pairs[:, 1]] - p[pairs[:, 0]]
        delta -= box * np.rint(delta / box)
        distances = np.sqrt(np.einsum('ij,ij->i', delta, delta))
        counts = np.histogram(distances, bins=EDGES)[0]
    else:
        counts = np.zeros(len(RADII))
    n = len(pos)
    expected = n * (n - 1) / 2 * SHELLS / np.linalg.det(cell)
    return counts / expected


def peak(curve):
    smooth = np.convolve(np.pad(curve, 2, mode='edge'), np.ones(5) / 5, mode='valid')
    indices = np.flatnonzero(WINDOW)
    k = indices[np.argmax(smooth[WINDOW])]
    if smooth[k] <= 0 or k in (indices[0], indices[-1]):
        return None, 'boundary_or_absent_peak'
    maxima, _ = find_peaks(smooth, height=.5 * smooth[k],
                          prominence=.1 * smooth[k], distance=3)
    maxima = maxima[WINDOW[maxima]]
    if len(maxima) != 1 or maxima[0] != k:
        return None, 'ambiguous_first_shell_peak'
    den = smooth[k - 1] - 2 * smooth[k] + smooth[k + 1]
    shift = .5 * (smooth[k - 1] - smooth[k + 1]) / den if den < 0 else 0
    return float(RADII[k] + np.clip(shift, -.5, .5) * BIN), 'clear'


def xdatcar_sample(path, selected):
    """Read fixed-cell XDATCAR; keep selected 0-based frames for the reference."""
    selected = set(selected)
    curves, indices = [], []
    with Path(path).open('rb') as f:
        f.readline()
        scale = float(f.readline())
        if scale <= 0:
            raise ValueError('Unexpected nonpositive VASP scaling factor')
        cell = np.array([f.readline().split() for _ in range(3)], dtype=float) * scale
        symbols = f.readline().decode().split()
        counts = list(map(int, f.readline().split()))
        key = comp_key(dict(zip(symbols, counts)))
        n = sum(counts)
        index = 0
        for marker in f:
            if not marker.strip():
                continue
            if not re.match(rb'(Direct|Cartesian)\s+configuration\s*=', marker):
                raise ValueError(f'Unexpected XDATCAR marker: {marker[:80]!r}')
            lines = [f.readline() for _ in range(n)]
            if any(not b for b in lines):
                raise ValueError('Truncated XDATCAR coordinate block')
            if index in selected:
                coords = np.fromstring(b''.join(lines).decode('ascii'), sep=' ').reshape(n, 3)
                pos = coords @ cell if marker.startswith(b'Direct') else coords * scale
                curves.append(rdf(pos, cell))
                indices.append(index)
            index += 1
    if len(curves) != len(selected):
        raise ValueError('Missing requested XDATCAR frames')
    return np.array(curves), indices, key, n, float(np.linalg.det(cell)), index


def reference_task(args):
    root, record = args
    path = Path(root) / record['path']
    if digest(path) != record['sha']:
        raise ValueError(f'Source changed: {path}')
    result = dict(record)
    if record['excluded_reason']:
        result.update(included=False)
        return result
    selected = np.linspace(record['frames'] - 1000, record['frames'] - 1, 200, dtype=int).tolist()
    curves, indices, key, n, volume, total = xdatcar_sample(path, selected)
    assert key == record['composition'] and total == record['frames']
    value, status = peak(curves.mean(0))
    first, _ = peak(curves[:100].mean(0))
    last, _ = peak(curves[100:].mean(0))
    drift = 100 * abs(first - last) / value if all(v is not None for v in [first, last, value]) else None
    reason = status if status != 'clear' else ('reference_peak_drift_ge5pct' if drift is None or drift >= 5 else '')
    result.update(included=not reason, excluded_reason=reason,
                  sampled_indices=indices, rdf_mean=curves.mean(0), peak_A=value,
                  first_half_peak_A=first, last_half_peak_A=last,
                  drift_pct=drift, n_atoms=n, volume_A3=volume)
    return result


def source_identity(path):
    p = path.split('#attempt=')[0]
    if '/outputs_v4/' in p:
        run = 'gese-v3-c6d320g_step400000_interp_w2.5/outputs_v4'
    elif '/local_eval/' in p:
        run = p.split('/local_eval/', 1)[1].split('/')[0]
    elif '/outputs/' in p and p.startswith('workspace_docs/supercomputer/'):
        run = '/'.join(p.split('/outputs/', 1)[1].split('/')[:2])
    else:
        run = str(Path(p).parent)
    method = 'flow_matching' if 'gese-v2-fm' in p else 'diffusion'
    return run, method


def source_task(args):
    root, source, candidates, groups = args
    path = Path(root) / source['path']
    if digest(path) != source['sha256']:
        raise ValueError(f'Source changed: {path}')
    wanted = {v['locator']: v for v in candidates}
    result = []

    def analyze(locator, raw, species, pos, cell):
        record = wanted[locator]
        key = comp_key(Counter(species))
        gid = geometry_id(species, pos, cell)
        assert gid == record['geometry'] and key == record['composition']
        assert len(species) == record['n']
        value, status = peak(rdf(pos, cell))
        refs = groups[key]
        errors = [100 * abs(value - r['peak_A']) / r['peak_A'] for r in refs] if value is not None else []
        worst = int(np.argmax(errors)) if errors else 0
        error = max(errors) if errors else None
        accepted = status == 'clear' and error is not None and error < 5
        run, method = source_identity(source['path'])
        row = dict(record, status='accepted' if accepted else (status if status != 'clear' else 'peak_error_ge5pct'),
                   first_peak_A=value, max_error_pct=error,
                   reference_group_ids=';'.join(r['id'] for r in refs),
                   error_reference_group_id=refs[worst]['id'], reference_peak_A=refs[worst]['peak_A'],
                   n_Ge=species.count('Ge'), n_Se=species.count('Se'),
                   volume_A3=float(np.linalg.det(cell)), source_id=source['id'],
                   source_run=run, generation_method=method,
                   raw_sha256=hashlib.sha256(raw).hexdigest())
        result.append((row, raw if accepted else None))

    if source['kind'] == 'sqlite':
        connection = sqlite3.connect(f'file:{path}?mode=ro', uri=True)
        for locator in sorted(wanted):
            row = connection.execute('select state,payload from attempts where id=?', (locator,)).fetchone()
            assert row and row[0] == 'complete' and row[1] is not None
            parsed = list(xyz_frames(io.BytesIO(zlib.decompress(row[1]))))
            assert len(parsed) == 1
            _, raw, species, pos, cell = parsed[0]
            analyze(locator, raw, species, pos, cell)
        connection.close()
    else:
        # Skip unused frames without parsing their arrays; scan the file once.
        with path.open('rb') as f:
            index = 0
            while line := f.readline():
                if not line.strip():
                    continue
                n = int(line)
                header = f.readline()
                body = [f.readline() for _ in range(n)]
                if any(not v for v in body):
                    raise ValueError('Truncated source extxyz')
                if index in wanted:
                    parsed = next(xyz_frames(io.BytesIO(line + header + b''.join(body))))
                    _, raw, species, pos, cell = parsed
                    analyze(index, raw, species, pos, cell)
                index += 1
                if len(result) == len(wanted):
                    break
    assert len(result) == len(candidates)
    return result


def write_csv(path, rows, fields):
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def reference_groups(references):
    by_group = defaultdict(list)
    for r in references:
        if r['included']:
            by_group[r['group_id']].append(r)
    groups = defaultdict(list)
    for gid, members in sorted(by_group.items()):
        key = members[0]['composition']
        value, status = peak(np.mean([r['rdf_mean'] for r in members], axis=0))
        if status != 'clear':
            raise ValueError(f'Unclear pooled reference peak: {gid}')
        group = dict(id=gid, peak_A=value, reference_ids=[r['id'] for r in members])
        groups[key].append(group)
        for r in members:
            r['group_peak_A'] = value
    return dict(groups)


def build(root, target, manifest, workers):
    if (target / 'structures.csv').exists() or (target / 'generated').exists():
        raise ValueError('Output already contains collected data; choose a new --output')
    target.mkdir(parents=True, exist_ok=True)
    log('Reference curves: fixed 200 frames from the final 1000 steps of each trajectory')
    with ProcessPoolExecutor(max_workers=workers) as pool:
        references = list(pool.map(reference_task, [(str(root), r) for r in manifest['references']]))
    groups = reference_groups(references)
    if set(groups) != {c['composition'] for c in manifest['candidates']}:
        raise ValueError('Missing a fixed composition reference group')
    log(f'References ready: {sum(r["included"] for r in references)} trajectories, {len(groups)} compositions')
    for key, values in sorted(groups.items()):
        log(f'{comp_name(key)} DFT peaks: ' + ', '.join(f'{v["id"]}={v["peak_A"]:.4f} A' for v in values))
    for r in references:
        if not r['included']:
            log(f'Reference excluded: {r["path"]}: {r["excluded_reason"]}')
            continue
        destination = Path('dft_reference') / comp_name(r['composition']) / ('XDATCAR_' + r['id'])
        (target / destination).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / r['path'], target / destination)
        r['output_file'] = str(destination)
    by_source = defaultdict(list)
    for candidate in manifest['candidates']:
        by_source[candidate['source_id']].append(candidate)
    sources = {s['id']: s for s in manifest['sources']}
    tasks = [(str(root), sources[sid], candidates, groups) for sid, candidates in sorted(by_source.items())]
    accepted, all_results, handles, indices = [], [], {}, Counter()
    t = time.time()
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for i, rows in enumerate(pool.map(source_task, tasks), 1):
                for row, raw in rows:
                    all_results.append({k: v for k, v in row.items() if k != 'raw_sha256'})
                    if raw is None:
                        continue
                    rel = Path('generated') / comp_name(row['composition']) / f'atoms_{row["n"]}.extxyz'
                    if rel not in handles:
                        (target / rel).parent.mkdir(parents=True, exist_ok=True)
                        handles[rel] = (target / rel).open('wb')
                    row['output_file'], row['output_frame'] = str(rel), indices[str(rel)]
                    row['structure_id'] = 'G_' + row['geometry'][:20]
                    row['source_file'] = sources[row['source_id']]['path']
                    row['source_sha256'] = sources[row['source_id']]['sha256']
                    row['source_locator_kind'] = 'attempt_id' if sources[row['source_id']]['kind'] == 'sqlite' else 'frame_index_0based'
                    handles[rel].write(raw)
                    indices[str(rel)] += 1
                    accepted.append(row)
                if i % 100 == 0 or i == len(tasks):
                    log(f'Candidate source files {i}/{len(tasks)}; checked={len(all_results)} accepted={len(accepted)}; {time.time()-t:.1f}s')
    finally:
        for f in handles.values():
            f.close()
    fields = ['structure_id', 'output_file', 'output_frame', 'n', 'n_Ge', 'n_Se', 'composition',
              'volume_A3', 'first_peak_A', 'reference_group_ids', 'error_reference_group_id',
              'reference_peak_A', 'max_error_pct', 'source_category', 'generation_method',
              'source_run', 'source_file', 'source_locator_kind', 'locator', 'source_sha256',
              'raw_sha256', 'geometry']
    write_csv(target / 'structures.csv', accepted, fields)
    ref_fields = ['id', 'composition', 'group_id', 'temperature_K', 'state_label', 'included',
                  'excluded_reason', 'output_file', 'path', 'sha', 'frames', 'n_atoms', 'volume_A3',
                  'sampled_indices', 'peak_A', 'group_peak_A', 'first_half_peak_A',
                  'last_half_peak_A', 'drift_pct']
    csv_refs = []
    for r in references:
        copy = {k: v for k, v in r.items() if k != 'rdf_mean'}
        if 'sampled_indices' in copy:
            copy['sampled_indices'] = json.dumps(copy['sampled_indices'], separators=(',', ':'))
        csv_refs.append(copy)
    write_csv(target / 'references.csv', csv_refs, ref_fields)
    manifest['outcomes'] = all_results
    manifest['summary'] = dict(candidates=len(all_results), accepted=len(accepted),
                               reasons=dict(Counter(r['status'] for r in all_results)),
                               accepted_by_atoms=dict(Counter(r['n'] for r in accepted)),
                               accepted_by_composition=dict(Counter(r['composition'] for r in accepted)),
                               included_references=sum(r['included'] for r in references),
                               reference_groups=groups,
                               output_sha256={p: digest(target / p) for p in indices})
    (target / 'collection_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, separators=(',', ':')) + '\n')
    if HERE != target:
        shutil.copy2(Path(__file__).resolve(), target / 'collect_dataset.py')
    verified = verify(target, recompute=True)
    write_readme(target, manifest, references, verified)
    log(f'COMPLETE: {len(accepted)} accepted structures, {len(indices)} generated files. {target}')


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def verify(target, recompute=True):
    manifest = json.loads((target / 'collection_manifest.json').read_text())
    rows, refs = read_csv(target / 'structures.csv'), read_csv(target / 'references.csv')
    groups = manifest['summary']['reference_groups']
    group_lookup = {g['id']: (key, g['peak_A']) for key, values in groups.items() for g in values}
    for r in refs:
        if r['included'] == 'True':
            assert digest(target / r['output_file']) == r['sha']
    grouped = defaultdict(dict)
    seen = set()
    for row in rows:
        assert row['geometry'] not in seen
        seen.add(row['geometry'])
        assert float(row['max_error_pct']) < 5
        assert all(group_lookup[gid][0] == row['composition'] for gid in row['reference_group_ids'].split(';'))
        grouped[row['output_file']][int(row['output_frame'])] = row
    count = 0
    for filename, expected in sorted(grouped.items()):
        assert digest(target / filename) == manifest['summary']['output_sha256'][filename]
        with (target / filename).open('rb') as f:
            found = 0
            for index, raw, species, pos, cell in xyz_frames(f):
                row = expected[index]
                assert hashlib.sha256(raw).hexdigest() == row['raw_sha256']
                assert geometry_id(species, pos, cell) == row['geometry']
                assert comp_key(Counter(species)) == row['composition']
                assert len(species) == int(row['n'])
                assert species.count('Ge') == int(row['n_Ge']) and species.count('Se') == int(row['n_Se'])
                if recompute:
                    value, status = peak(rdf(pos, cell))
                    assert status == 'clear' and abs(value - float(row['first_peak_A'])) < 1e-10
                    errors = [100 * abs(value - group_lookup[gid][1]) / group_lookup[gid][1]
                              for gid in row['reference_group_ids'].split(';')]
                    assert max(errors) < 5 and abs(max(errors) - float(row['max_error_pct'])) < 1e-10
                found += 1
            assert found == len(expected)
            count += found
    assert count == len(rows) == manifest['summary']['accepted']
    log(f'VERIFIED {count} structures: exact composition pairing, original frame bytes, geometry, peak errors and reference copy hashes')
    return dict(structures=count, generated_files=len(grouped), reference_files=sum(r['included'] == 'True' for r in refs))


def write_readme(target, manifest, references, verified):
    s = manifest['summary']
    outcomes = manifest['outcomes']
    structure_rows = read_csv(target / 'structures.csv')
    methods = Counter(r['generation_method'] for r in structure_rows)
    endpoints = Counter(r['composition'] for r in structure_rows if r['composition'] in ['Ge:1', 'Se:1'])
    lines = ['# 生成结构与DFT参考：第一峰位配对数据集', '',
             '本地收集于2026-11-16。生成结构、DFT参考分别存放；按真实元素数的约分比例严格对应，允许原子数不同。', '',
             f'共检查 **{s["candidates"]:,}** 个去重候选，收集 **{s["accepted"]:,}** 个第一峰位误差严格小于5%的结构。', '',
             f'其中扩散生成 **{methods["diffusion"]:,}** 帧，Flow Matching **{methods["flow_matching"]:,}** 帧；可按structures.csv中的generation_method列区分。纯Ge {endpoints["Ge:1"]:,}帧、纯Se {endpoints["Se:1"]:,}帧，其余Ge–Se合金 {len(structure_rows)-sum(endpoints.values()):,}帧。', '',
             '| 原子数 | 候选 | 收集 |', '|---:|---:|---:|']
    for n in [100, 300, 1000, 30000]:
        lines.append(f'| {n:,} | {sum(r["n"]==n for r in outcomes):,} | {sum(r["n"]==n and r["status"]=="accepted" for r in outcomes):,} |')
    lines += ['', '| 组分 | 候选 | 收集 |', '|---|---:|---:|']
    keys = sorted({r['composition'] for r in outcomes}, key=lambda k: sum(int(v) for e,v in (p.split(':') for p in k.split(';')) if e=='Ge') / sum(int(p.split(':')[1]) for p in k.split(';')))
    for key in keys:
        lines.append(f'| {comp_name(key)} | {sum(r["composition"]==key for r in outcomes):,} | {s["accepted_by_composition"].get(key,0):,} |')
    lines += ['', '## 查看与使用', '',
              '- `generated/<组分>/atoms_<原子数>.extxyz`：通过验收的生成结构，按组分与尺寸合并；保留原始每帧字节。',
              '- `dft_reference/<组分>/XDATCAR_<ID>`：独立的原始DFT参考轨迹。',
              '- `structures.csv`：每行对应一个已收集结构；`output_frame`和源文件帧索引均从0开始，数据库来源的`locator`是attempt记录ID。',
              '- `references.csv`：参考轨迹、实际组分、已知温度、使用帧号、独立峰位和参考组峰位。`included=False`的轨迹未用于基准、未复制。',
              '- `collection_manifest.json`：机器读取的固定候选范围、源文件哈希、排除原因、参数和结果；一般无需人工浏览。', '',
              '## 验收定义', '',
              '总数目g(r)，0.05 Å bin，按晶胞体积和有限N原子对数归一；第一配位壳2.0–3.3 Å，统一5点滑动平均，主峰三点插值。计算邻近原子对，未创建大胞的N×N矩阵。边界峰、缺峰和明显多峰暂不收集。', '',
              '`error_pct = 100 * abs(r_generated - r_DFT) / r_DFT < 5`。该指标是第一峰位误差，不代表整条g(r)曲线误差或全部结构性质的准确性。未按SOAP等额外压缩多样性。', '',
              '参考用每条轨迹最后1000步内均匀抽出的200帧，完整2000帧原轨迹全部保留。相同组分的新300 K多seed轨迹按轨迹等权平均g(r)建立一个参考组；旧轨迹缺少充分条件记录时各自独立成组，不与新300 K轨迹混合。生成结构必须对其组分的全部参考组均通过5%；索引的`max_error_pct`保存最差误差，`error_reference_group_id`与`reference_peak_A`指产生该最差误差的组。参考组在筛选前确定，不逐帧挑更有利的参考。', '',
              '参考稳定性使用后半段前100/后100个抽样帧的均值峰位比较，漂移≥5%或无法明确定位者不作为基准。旧轨迹温度未查明时留空并标记legacy，不自动标作300 K。', '',
              f'实际包含DFT参考{verified["reference_files"]}条；下列参考暂未纳入：']
    for r in references:
        if not r['included']:
            lines.append(f'- `{r["path"]}`：{r["excluded_reason"]}。')
    lines += ['', '## 来源与数量口径', '',
              '候选来自本地新旧生成模型、归档版本、消融与Flow Matching等分支，来源在索引中保留，不将这些结果全部当作当前主模型的表现。训练POSCAR、MLP弛豫轨迹、DFT优化中间步、骨架截取和未训练冒烟测试未混作独立生成样本。', '',
              '原工作区逐帧初步去重候选93,733个，其中14,335个与DFT严格对应；无对应者不进入本包。纯Ge/Se也计入总量；参考轨迹帧数不计为生成数量。', '',
              '筛选结果：']
    for reason, count in sorted(s['reasons'].items()):
        lines.append(f'- `{reason}`：{count:,}。')
    lines += ['', '## 验证与复现', '',
              '依赖Python≥3.10、NumPy、SciPy。使用现有chaldiff环境即可：', '',
              '```bash',
              '/home/ysj/anaconda3/envs/chaldiff/bin/python -B collect_dataset.py verify',
              '# 有原始工作区时，按固定清单在另一个新目录完整重建：',
              '/home/ysj/anaconda3/envs/chaldiff/bin/python -B collect_dataset.py build \\',
              '  --source-root /home/ysj/amor_diff/amor_gen/graphite-main \\',
              '  --output /tmp/diffusion_rdf_dataset_rebuild --workers 4',
              '```', '',
              'verify无需原始工作区，逐帧重算峰位、核对严格组分对应和字节/几何哈希。build核对原始源文件哈希，使用固定候选范围和参考规则重建；拒绝覆盖已有生成数据。完整参考基准重算可使用verify --references。', '',
              f'收集后已核对{verified["structures"]:,}帧与索引一致，所有原始帧字节、参考文件哈希以及逐帧峰位复算通过。原工作区文件和public_release未改动，未上传。', '']
    (target / 'README.md').write_text('\n'.join(lines))


def verify_reference_curves(target):
    rows = read_csv(target / 'references.csv')
    pooled = defaultdict(list)
    for r in rows:
        if r['included'] != 'True':
            continue
        curves, _, key, _, _, _ = xdatcar_sample(target / r['output_file'], json.loads(r['sampled_indices']))
        assert key == r['composition']
        value, status = peak(curves.mean(0))
        assert status == 'clear' and abs(value - float(r['peak_A'])) < 1e-10
        pooled[r['group_id']].append(curves.mean(0))
    expected = {r['group_id']: float(r['group_peak_A']) for r in rows if r['included'] == 'True'}
    for gid, curves in pooled.items():
        value, status = peak(np.mean(curves, axis=0))
        assert status == 'clear' and abs(value - expected[gid]) < 1e-10
    log(f'VERIFIED DFT reference curves and all {len(pooled)} fixed groups from copied raw XDATCAR')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['build', 'verify'])
    parser.add_argument('--source-root', type=Path)
    parser.add_argument('--output', type=Path, default=HERE)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--references', action='store_true', help='Recompute all DFT reference curves during verify')
    args = parser.parse_args()
    if args.command == 'verify':
        verify(args.output)
        if args.references:
            verify_reference_curves(args.output)
    else:
        manifest = json.loads((HERE / 'collection_manifest.json').read_text())
        root = (args.source_root or Path(manifest['source_root'])).resolve()
        build(root, args.output.resolve(), manifest, args.workers)


if __name__ == '__main__':
    main()
