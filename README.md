# SUEWS

This is a private repo for SUEWS development.

## Documentation

* Documentation site: <https://suews-docs.readthedocs.io/>

* Documentation repo: <https://github.com/Urban-Meteorology-Reading/SUEWS-Docs>

## Developer Note

### Branch

#### Central curated branches
These branches are regularly curated by admin members with specific purposes and set with triggers for automatic deployment (via MS Azure Pipeline) in the [*releases* page](https://github.com/Urban-Meteorology-Reading/SUEWS/releases) named **Latest Release Test**:

* `master`:
  * the main branch that keeps milestone and stable features.
  * `push` is restricted to admin members.
* `PublicRelease`:
  * used for public releases: compiled binaries are archived and published to the public
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

### Questions

* Please [raise issues](https://github.com/Urban-Meteorology-Reading/SUEWS/issues/new) for questions in the development so our progress can be well managed.
