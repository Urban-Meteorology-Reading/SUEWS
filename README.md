# SUEWS

This is a public repo for SUEWS source and documentation. 

## Documentation

* Documentation site: <https://suews-docs.readthedocs.io/>

* Documentation source: `docs` folder in this repo

## Developer Note

### Branch

#### Central curated branches
These branches are regularly curated by admin members with specific purposes and set with triggers for automatic deployment (via MS Azure Pipeline) in the [*releases* page](https://github.com/Urban-Meteorology-Reading/SUEWS/releases) named **Latest Release Test**:

* `master`:
  * the main branch that keeps milestone and stable features.
  * `push` is restricted to admin members.
* `PublicRelease`:
  * used for public releases: compiled binaries are archived and published to the public.
  * version information (e.g., 2018c) will be added as a tag once published and a formal release will be available via the `Releases` tab.
  * `push` is restricted to admin members.
* `test-dev`:
  * used for the fast development and test
  * new features should be incorporated into this branch for testing before merging into `master`

#### Feature branches
Branches for specific features are NOT maintained by admin members: they are proposed and maintained by related developers. Being without central curation, developers are strongly suggested to following certain rules:

* using the latest `test-dev` as the basis to conduct feature development for easier merging at a later time.
* naming the feature branch as `test-dev-{feature}` so certain branch maintenance rules can be easily applied by curators whenever necessary.
* testing the features as thorough as possible; details refer to the [Test](###Test) section.

#### General merging workflow

`test-dev-{feature}` --[feature ready]--> `test-dev` --[central test passed]--> `master` --[release ready]--> `PublicRelease`



### Manual

* Please keep the development changes documented in the [Documentation repo](https://github.com/Urban-Meteorology-Reading/SUEWS-Docs).
* The SUEWS docs are mainly written in [RST](http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html) but [Markdown](https://guides.github.com/features/mastering-markdown/) is also acceptable.
* Your changes to docs will be reviewed and merged into public release if appropriate.

### Test

Whenever changes are made, please run `make test` in the repo root to check if your changes are working or not.
If any error, please resolve it or justify that the test is giving false alarm.

#### Tests and purposes
`make test` will perform the following checks:

- if single-grid-multi-year run could be successful: to check if multi-year run is functional;
- if multi-grid-multi-year run could be successful: to check if multi-grid run is functional;
- if test run results could match those from the standard run (runs under `BaseRun`): to check if any non-functional changes would break the current code;
- if all physics schemes are working: to check possible invalid physics schemes.

#### Workflow
The test workflow is as follows (details refer to the Makefile `test` recipe and related python code):

1. clean existing build and rebuild the code;
2. create a temporary directory as working directory to perform checks;
3. copy the rebuilt `SUEWS_{version}` binary to the temporary folder;
4. copy the version specific input files under `Release/InputTables/{version}` to the temporary folder (see below for its preparation);
5. run python code to perform the above checks and write out test results to the console:
   1. if all tests are successful, the code is generally good to go;
   2. if any test failed, we NEED to look into the reasons for the failure and resolve them before any further feature-added operations.

#### Preparation of tests

1. Prepare a base run:
   - under `Test/BaseRun`, create a folder named with version/feature info (e.g., `2019a`);
   - perform a simulation to produce example output files, which will later be used as standard run to verify the correct code functionalities.

   *Note: all the above input files will be automatically copied under `Release/InputTables` with explicit version/feature (e.g., `Release/InputTables/2019a`) and later archived in public releases for users; so carefully construct test data to include in the input files.*
2. Configure test namelist file `Test/code/BTS_config.nml`:

   - `name_exe`: the SUEWS binary name that will be used for testing.
   - `dir_exe`: the directory to copy `name_exe`.
   - `dir_input`: the directory to copy input files; suggested to be `Release/InputTables/{version}`.
   - `dir_baserun`: the base run against which to test identity in results.

### Questions

* Please [raise issues](https://github.com/Urban-Meteorology-Reading/SUEWS/issues/new) for questions in the development so our progress can be well managed.
