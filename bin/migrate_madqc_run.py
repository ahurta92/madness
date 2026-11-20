#!/usr/bin/env python3
"""Relocate a stopped SCF/response development run to the shared root layout.

Stop every solver using this directory first. Moves use same-filesystem rename;
archives and evidence are never rewritten. Run the same command after interruption.
This is a one-time development utility, not a solver compatibility mode.
"""
import argparse
import json
import os
from pathlib import Path
import tempfile

MARKER = '.madqc-layout-migration.json'
OLD_SCF = Path('task_0/moldft')
OLD_RESPONSE = Path('task_1/molresponse')


def destination(source):
    p = Path(source)
    if p.is_absolute() or '..' in p.parts:
        raise ValueError('invalid migration source')
    if p.parent == OLD_SCF:
        return Path(p.name) if p.name in ('checkpoints', 'results.json', 'timings.json') else Path('analysis/ground_state')/p.name
    if p.parent in (OLD_RESPONSE, OLD_RESPONSE/'analysis') and p.name != 'analysis':
        return Path('analysis/response')/p.name
    raise ValueError('unrecognized migration source: '+str(p))


def reject_links(root):
    if any(p.is_symlink() for p in root.rglob('*')):
        raise ValueError('symlinks in the run directory require manual relocation')


def prepare(root):
    reject_links(root)
    old = root/OLD_SCF
    if not (old/'checkpoints/run.json').is_file():
        raise ValueError('no previous SCF/response shared checkpoint store found')
    if {p.name for p in root.glob('task_*')} - {'task_0', 'task_1'}:
        raise ValueError('extra workflow tasks require manual relocation')
    for task, engine in [('task_0', 'moldft'), ('task_1', 'molresponse')]:
        directory = root/task
        if directory.exists() and {p.name for p in directory.iterdir()} != {engine}:
            raise ValueError('unexpected contents in '+str(directory))
    for name in ['checkpoints', 'results.json', 'timings.json', 'analysis']:
        if (root/name).exists():
            raise ValueError('destination already exists: '+str(root/name))
    moves = [str(p.relative_to(root)) for p in sorted(old.iterdir())]
    response = root/OLD_RESPONSE
    if response.exists():
        moves += [str(p.relative_to(root)) for p in sorted(response.iterdir()) if p.name != 'analysis']
        if (response/'analysis').exists():
            moves += [str(p.relative_to(root)) for p in sorted((response/'analysis').iterdir())]
    targets = [destination(p) for p in moves]
    if len(targets) != len(set(targets)):
        raise ValueError('conflicting response diagnostics require manual relocation')
    return {'version': 1, 'moves': moves}


def sync_directory(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_plan(root, plan):
    # A durable marker stops the new solver throughout relocation. No archive
    # moves happen until this complete plan has been atomically published.
    fd, name = tempfile.mkstemp(prefix='.madqc-layout-', suffix='.tmp', dir=root)
    temp = Path(name)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(plan, stream, indent=2)
            stream.flush()
            os.fsync(stream.fileno())
        # link is atomic and refuses to overwrite another relocation plan.
        os.link(temp, root/MARKER)
        sync_directory(root)
    finally:
        temp.unlink()


def migrate(root):
    root = Path(root).resolve(strict=True)
    reject_links(root)
    for store in [root/OLD_SCF/'checkpoints', root/'checkpoints']:
        if (store/'.store-lock').exists():
            raise ValueError('active or stale checkpoint lock; stop writers and resolve the lock first')
    marker = root/MARKER
    plan = json.loads(marker.read_text()) if marker.exists() else prepare(root)
    if plan.get('version') != 1 or not isinstance(plan.get('moves'), list) or not plan['moves']:
        raise ValueError('invalid relocation plan')
    pairs = [(root/p, root/destination(p)) for p in plan['moves']]
    if len({s for s, _ in pairs}) != len(pairs) or len({d for _, d in pairs}) != len(pairs):
        raise ValueError('duplicate relocation entries')
    for source, target in pairs:
        if source.exists() == target.exists():
            raise ValueError('expected exactly one source or destination: '+str(source))
    if not marker.exists():
        write_plan(root, plan)
    for source, target in pairs:
        if source.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                raise ValueError('destination appeared during relocation: '+str(target))
            source.rename(target)
            sync_directory(source.parent)
            sync_directory(target.parent)
    # Never recursively remove anything: unexpected files leave the marker and
    # a recoverable error, with every completed rename still in the plan.
    for relative in [OLD_RESPONSE/'analysis', OLD_RESPONSE, Path('task_1'), OLD_SCF, Path('task_0')]:
        directory = root/relative
        if directory.exists():
            directory.rmdir()
            sync_directory(directory.parent)
    marker.unlink()
    sync_directory(root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    args = parser.parse_args()
    try:
        migrate(args.run_directory)
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.exit(1, f'Relocation stopped: {error}\n')
    print(f'Relocated {args.run_directory}. Restart with the original madqc command.')
    print('Historical analysis may contain old path text; native source identities are unchanged.')


if __name__ == '__main__':
    main()
