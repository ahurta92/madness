"""A named test's work directory, removed when the test passes (D-108).

work_directory makes the directory where the tests always made it (a
mkdtemp in the current directory, the build tree under ctest). When the
script ends without an uncaught exception (a failed assert, a failed
subprocess) the directory is removed; otherwise it is kept and its path
printed, so a failure can be inspected. Left behind, the directories of
a full named suite reach hundreds of GB.
"""
import atexit
import shutil
import sys
import tempfile
from pathlib import Path


def work_directory(prefix):
    work = Path(tempfile.mkdtemp(prefix=prefix, dir=Path.cwd()))
    failed = []
    previous = sys.excepthook

    def excepthook(*exc_info):
        failed.append(exc_info[0])
        previous(*exc_info)

    def finish():
        if failed:
            print(f'test failed: work directory kept at {work}', file=sys.stderr, flush=True)
        else:
            shutil.rmtree(work)

    sys.excepthook = excepthook
    atexit.register(finish)
    return work
