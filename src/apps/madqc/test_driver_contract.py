#!/usr/bin/env python3
"""One public driver: old named input cannot start work or publish files."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', required=True)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--mpiexec', default='mpiexec')
    args = parser.parse_args()
    env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1', OMPI_MCA_hwloc_base_binding_policy='none')
    launch = [] if args.ranks == 1 else [args.mpiexec, '-n', str(args.ranks)]
    binary = str(Path(args.binary).resolve())
    with tempfile.TemporaryDirectory(prefix='madqc-driver-') as scratch:
        work = Path(scratch)
        (work/'named.in').write_text('''molecule helium
 units atomic
 He 0 0 0
end
moldft ground
 molecule helium
 xc hf
 protocol [1e-4]
 maxiter 1
 requested_properties [energy]
end
''')
        for deck, extra in [('named.in', []), ('absent.in', ['--stage=ground']),
                            ('absent.in', ['--from-store=previous'])]:
            result = subprocess.run(launch+[binary,deck,'--prefix=forbidden']+extra,
                                    cwd=work,env=env,capture_output=True,text=True,timeout=60)
            assert result.returncode == 1, result.stdout+result.stderr
            assert not (work/'forbidden').exists(), 'retired interface started work or published output'
            assert 'unknown input' in result.stdout.lower(), result.stdout+result.stderr
        for name, contents in [
                ('uppercase.in', 'DFT\n xc hf\nend\n'),
                ('missing_end.in', 'io\nresponse\n dipole.directions z\nend\n')]:
            (work/name).write_text(contents)
            result = subprocess.run(launch+[binary,name,'--geometry=he','--prefix=forbidden'],
                cwd=work,env=env,capture_output=True,text=True,timeout=60)
            assert result.returncode == 1, result.stdout+result.stderr
            assert not (work/'forbidden').exists(), 'malformed group input started computation'
            assert 'unknown input' in result.stdout.lower(), result.stdout+result.stderr
        for key in ['dconv', 'excited.guess_max_iter']:
            rejected = subprocess.run(launch+[binary,'--wf=response','--geometry=he',
                '--prefix=forbidden','--response='+key+'=1','--response=maxiter=0'],
                cwd=work,env=env,capture_output=True,text=True,timeout=30)
            assert rejected.returncode == 1, rejected.stdout+rejected.stderr
            assert 'unknown input' in (rejected.stdout+rejected.stderr).lower()
            assert not (work/'forbidden').exists()
        # A geometry load failure on rank zero must release every rank.
        bad = subprocess.run(launch+[binary,'--geometry=missing_geometry_for_test',
            '--prefix=forbidden'],cwd=work,env=env,capture_output=True,text=True,timeout=20)
        assert bad.returncode == 1, bad.stdout+bad.stderr
        assert not (work/'forbidden').exists()
        # Preserve CLI-only calls and workflow alias precedence. The last --w
        # must not override the higher-precedence --workflow selector.
        result = subprocess.run(launch+[binary,'--workflow=scf','--w=not_a_workflow',
            '--geometry=he','--prefix=CaseRun','--dft=protocol=[1e-4];maxiter=0;convergence.criteria=[bsh_residual];convergence.bsh_residual=1000'],
            cwd=work,env=env,capture_output=True,text=True,timeout=90)
        assert result.returncode == 0, result.stdout+result.stderr
        assert (work/'CaseRun/results.json').exists()
        # JSON/file input keeps CLI override precedence and case-sensitive paths.
        (work/'MixedCase.json').write_text(json.dumps({'dft': {
            'maxiter': 0, 'protocol': [1e-4],
            'convergence.criteria': ['bsh_residual'],
            'convergence.bsh_residual': 0}}))
        result = subprocess.run(launch+[binary,'MixedCase.json','--wf=scf',
            '--geometry=he','--prefix=JsonRun','--dft=convergence.bsh_residual=1000'],
            cwd=work,env=env,capture_output=True,text=True,timeout=90)
        assert result.returncode == 0, result.stdout+result.stderr
        assert (work/'JsonRun/results.json').exists()
        (work/'nested.in').write_text('geometry\nmolecule\n He 0 0 0\nend\nend\n')
        result = subprocess.run(launch+[binary,'nested.in','--wf=scf','--prefix=nested',
            '--dft=protocol=[1e-4];maxiter=0;convergence.criteria=[bsh_residual];convergence.bsh_residual=1000'],
            cwd=work,env=env,capture_output=True,text=True,timeout=90)
        assert result.returncode == 0, result.stdout+result.stderr
        assert (work/'nested/results.json').exists()
        result = subprocess.run(launch+[binary,'--restart_info=missing-checkpoint'],
            cwd=work,env=env,capture_output=True,text=True,timeout=30)
        assert result.returncode == 1, result.stdout+result.stderr
    print('madqc driver contract passed')


if __name__ == '__main__':
    main()
