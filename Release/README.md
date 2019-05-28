# Release

This folder includes necessary files for public releases:

1. `InputTables`: version specific input files.

    *Note: Don't modify these files as they are automatically prepared via the test procedure. Chnages to these files are should be made in the `Test/BaseRun/{version}` folder.*

2. `pack_release.py`: script to prepare the release archive.
3. `Makefile`: makefile to conveniently call the above packing script.