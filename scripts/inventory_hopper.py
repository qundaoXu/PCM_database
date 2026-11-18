#!/usr/bin/env python3
"""Read-only Hopper discovery/counting; remote script uses only Python stdlib.

Run locally with the pcm_database environment. Reuses an existing SSH control
socket. Saves full machine-readable results locally, never writes on Hopper.
"""
import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REMOTE = r'''
import os, sys, json, re, math, collections, mmap, time, hashlib, base64
roots = ['/scratch/t0937929/MLIP_benchmark', '/scratch/t0937929/deeph',
         '/scratch/t0937929/myproject/1-GST225-aging']
skip_dirs = {'.git', '__pycache__', 'node_modules', 'DeepH-dock', 'site-packages',
             '.cache', '.venv', 'venv'}
if config['mode'] == 'discover':
    files, excluded, errors = [], [], []
    def walk_error(error):
        errors.append(str(error))
    for root in roots:
        for parent, dirs, names in os.walk(root, onerror=walk_error):
            for d in list(dirs):
                if d in skip_dirs:
                    excluded.append(os.path.join(parent, d)); dirs.remove(d)
            for name in names:
                ext = os.path.splitext(name)[1].lower()
                if not (ext in {'.xyz', '.extxyz', '.dump', '.lammpstrj', '.traj', '.lammps'}
                        or 'xdatcar' in name.lower() or name in {'POSCAR','CONTCAR'}
                        or name.lower().startswith(('dump.', 'trajectory.', 'traj.'))):
                    continue
                path = os.path.join(parent, name)
                try:
                    stat = os.stat(path)
                    with open(path, 'rb') as f:
                        head = f.read(1500)
                except OSError as e:
                    errors.append(str(e)); continue
                if b'ITEM: TIMESTEP' in head[:100]: kind = 'lammps-dump-text'
                elif 'xdatcar' in name.lower(): kind = 'vasp-xdatcar'
                elif name in {'POSCAR', 'CONTCAR'}: kind = 'vasp'
                elif ext in {'.xyz', '.extxyz'} and re.match(rb'\s*\d+\s*\n', head): kind = 'xyz'
                elif head.startswith(b'- of Ulm'): kind = 'ase-traj'
                else: continue
                files.append(dict(path=path, bytes=stat.st_size, mtime=stat.st_mtime,
                                  format=kind, header=head.decode('utf-8','replace')))
    print(json.dumps(dict(roots=roots, files=files, excluded_directories=excluded,
                          errors=errors, mode='read_only_discovery', scanned_at_unix=time.time())))
elif config['mode']=='samples':
    samples=[]
    for r in config['files']:
        with open(r['path'],'rb') as f:
            header=[f.readline() for _ in range(9)]
            if header[0].strip()!=b'ITEM: TIMESTEP': raise ValueError('Not a dump')
            n=int(header[3]); lines=[f.readline() for _ in range(n)]
            if any(not line for line in lines): raise ValueError('Incomplete first frame')
        samples.append(dict(path=r['path'],first_frame_base64=base64.b64encode(b''.join(header+lines)).decode()))
    print(json.dumps(dict(files=samples,mode='read_only_first_frame_samples')))
else:
    results=[]
    for record in config['files']:
        path, kind = record['path'], record['format']
        row=dict(record); row.pop('header',None)
        before=os.stat(path)
        try:
            with open(path,'rb') as f:
                if kind=='lammps-dump-text':
                    count=0; natoms=set(); steps=[]; previous=None; strides=collections.Counter()
                    with mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as mm:
                        for marker in re.finditer(rb'(?m)^ITEM: TIMESTEP\r?\n',mm):
                            mm.seek(marker.end()); step=int(mm.readline()); steps.append(step)
                            if previous is not None: strides[step-previous]+=1
                            previous=step
                            if mm.readline().strip()!=b'ITEM: NUMBER OF ATOMS': raise ValueError('No atom-count marker')
                            n=int(mm.readline()); natoms.add(n)
                            if not mm.readline().startswith(b'ITEM: BOX BOUNDS'): raise ValueError('No box bounds')
                            for _ in range(3): mm.readline()
                            if not mm.readline().startswith(b'ITEM: ATOMS'): raise ValueError('No atom columns')
                            # Validate block boundary and completeness without parsing all coordinates.
                            for _ in range(n):
                                line=mm.readline()
                                if not line or line.startswith(b'ITEM:'): raise ValueError('Incomplete atom block')
                            end=mm.tell()
                            nextline=mm.readline()
                            if nextline and nextline.strip()!=b'ITEM: TIMESTEP': raise ValueError('Unexpected block boundary')
                            count+=1
                    row.update(frames=count,atom_counts=sorted(natoms),first_timestep=steps[0] if steps else None,
                               last_timestep=steps[-1] if steps else None, timestep_differences=dict(strides),
                               count_method='all_timestep_markers_and_complete_coordinate_blocks')
                elif kind=='xyz':
                    count=0; natoms=set(); species_counts=collections.Counter()
                    while True:
                        line=f.readline()
                        if not line: break
                        if not line.strip(): continue
                        n=int(line); natoms.add(n); header=f.readline()
                        if not header: raise ValueError('Missing xyz header')
                        counts=collections.Counter()
                        for _ in range(n):
                            atom=f.readline().split()
                            if len(atom)<4: raise ValueError('Incomplete xyz atom block')
                            counts[atom[0].decode()]+=1
                        key=';'.join('%s:%d'%(s,k) for s,k in sorted(counts.items()))
                        species_counts[key]+=1; count+=1
                    row.update(frames=count,atom_counts=sorted(natoms),species_counts=dict(species_counts),
                               count_method='all_complete_xyz_blocks')
                elif kind in {'vasp','vasp-xdatcar'}:
                    head=[f.readline() for _ in range(7)]
                    n=sum(map(int,head[6].split())); symbols=head[5].decode().split()
                    row.update(atom_counts=[n],species_counts=dict(zip(symbols,map(int,head[6].split()))))
                    if kind=='vasp': row.update(frames=1,count_method='vasp_header_only')
                    else:
                        count=0
                        while True:
                            line=f.readline()
                            if not line: break
                            if not line.strip(): continue
                            if b'configuration=' not in line: raise ValueError('Unexpected XDATCAR block')
                            for _ in range(n):
                                if not f.readline(): raise ValueError('Incomplete XDATCAR')
                            count+=1
                        row.update(frames=count,count_method='all_complete_xdatcar_blocks')
                else:
                    row.update(frames=None,count_method='unsupported_without_remote_ASE')
            row['status']='counted' if row.get('frames') is not None else 'pending'
        except Exception as e:
            row.update(status='error',error=str(e),frames=None)
        after=os.stat(path)
        with open(path,'rb') as hash_file:
            h=hashlib.sha256()
            for block in iter(lambda:hash_file.read(4*1024*1024),b''):
                h.update(block)
            row['sha256']=h.hexdigest()
        row['changed_during_scan']=(before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns)
        if row['changed_during_scan']: row['status']='active_file_snapshot'
        results.append(row)
        print('counted %d/%d %s'%(len(results),len(config['files']),path),file=sys.stderr,flush=True)
    print(json.dumps(dict(files=results,mode='read_only_complete_block_count',
                         note='Counts do not establish phase, independence, model accuracy, or uniqueness',
                         scanned_at_unix=time.time())))
'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['discover', 'count', 'samples'])
    p.add_argument('--selection', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--socket', default='/private/tmp/pcm_database_hopper_20261117.sock')
    args = p.parse_args()
    config = dict(mode=args.mode)
    if args.mode in {'count', 'samples'}:
        config['files'] = json.loads(args.selection.read_text())['files']
    script = 'config = ' + repr(config) + '\n' + REMOTE
    result = subprocess.run(['ssh', '-S', args.socket, '-o', 'BatchMode=yes',
                             'hopper.nus.edu.sg', 'python3 -'], input=script,
                            text=True, stdout=subprocess.PIPE, check=True)
    data = json.loads(result.stdout)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(args.output)
    print(f'Saved {len(data["files"])} files to {args.output}')


if __name__ == '__main__':
    main()
