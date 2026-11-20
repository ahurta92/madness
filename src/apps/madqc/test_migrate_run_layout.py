#!/usr/bin/env python3
"""Offline layout migration preserves bytes and can resume an interrupted rename."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('migration', Path(__file__).resolve().parents[3]/'bin/migrate_madqc_run.py')
migration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(migration)

class MigrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)/'run'
        self.old = self.root/'task_0/moldft'
        (self.old/'checkpoints').mkdir(parents=True)
        (self.old/'checkpoints/run.json').write_text('{"version":2}')
        (self.old/'checkpoints/native.00000').write_bytes(b'original native data')
        (self.old/'results.json').write_text('{"original":true}')
        (self.old/'scf.dat').write_text('scf analysis')
        self.response = self.root/'task_1/molresponse'
        (self.response/'analysis').mkdir(parents=True)
        (self.response/'analysis/history.jsonl').write_text('metrics')
        (self.response/'response_metadata.json').write_text('{}')

    def test_relocation(self):
        inode = (self.old/'checkpoints/native.00000').stat().st_ino
        migration.migrate(self.root)
        self.assertEqual((self.root/'checkpoints/native.00000').stat().st_ino, inode)
        self.assertEqual((self.root/'checkpoints/native.00000').read_bytes(), b'original native data')
        self.assertEqual((self.root/'results.json').read_text(), '{"original":true}')
        self.assertEqual((self.root/'analysis/ground_state/scf.dat').read_text(), 'scf analysis')
        self.assertEqual((self.root/'analysis/response/history.jsonl').read_text(), 'metrics')
        self.assertFalse(list(self.root.glob('task_*')))
        self.assertFalse((self.root/migration.MARKER).exists())

    def test_no_overwrite(self):
        (self.root/'results.json').write_text('unrelated')
        with self.assertRaises(ValueError): migration.migrate(self.root)
        self.assertEqual((self.root/'results.json').read_text(), 'unrelated')
        self.assertTrue((self.old/'checkpoints/native.00000').exists())
        self.assertFalse((self.root/migration.MARKER).exists())

    def test_symlink_rejected(self):
        (self.old/'external').symlink_to(self.root.parent, target_is_directory=True)
        with self.assertRaises(ValueError): migration.migrate(self.root)
        self.assertFalse((self.root/migration.MARKER).exists())

    def test_locked_or_extra_task_rejected(self):
        lock = self.old/'checkpoints/.store-lock'
        lock.mkdir()
        with self.assertRaises(ValueError): migration.migrate(self.root)
        lock.rmdir()
        (self.root/'task_2').mkdir()
        with self.assertRaises(ValueError): migration.migrate(self.root)
        self.assertFalse((self.root/migration.MARKER).exists())

    def test_resume_after_rename(self):
        plan = migration.prepare(self.root)
        migration.write_plan(self.root, plan)
        src = self.root/plan['moves'][0]
        dst = self.root/migration.destination(plan['moves'][0])
        dst.parent.mkdir(parents=True, exist_ok=True)
        src.rename(dst)  # killed after a rename, before completion
        migration.migrate(self.root)
        self.assertEqual((self.root/'checkpoints/native.00000').read_bytes(), b'original native data')
        self.assertFalse((self.root/migration.MARKER).exists())

    def test_invalid_resume_plan(self):
        (self.root/migration.MARKER).write_text(json.dumps({'version':1, 'moves':['../unrelated']}))
        with self.assertRaises(ValueError): migration.migrate(self.root)
        self.assertTrue((self.old/'checkpoints/native.00000').exists())

    def test_duplicate_response_filename_rejected_before_changes(self):
        (self.response/'history.jsonl').write_text('other history')
        with self.assertRaises(ValueError): migration.migrate(self.root)
        self.assertFalse((self.root/migration.MARKER).exists())

if __name__ == '__main__': unittest.main()
