# SUEWS

This is a private repo for SUEWS development.

## Documentation

* Documentation site: https://suews-docs.readthedocs.io/

* Documentation repo: https://github.com/Urban-Meteorology-Reading/SUEWS-Docs


## Developer Note

### Branch

These branches are designed with specific purposes and set with triggers for automatic deployment (via MS Azure Pipeline) in the [*releases* page](https://github.com/Urban-Meteorology-Reading/SUEWS/releases) named **Latest Release Test**:

* `master`: the main branch that keeps milestone and stable features
* `PublicRelease`: the branch used for public releases that are archived and published to the public
* `test-dev`: the fast development and test branch used for new development and test

### Manual

* Please keep the development changes documented in the [Documentation repo](https://github.com/Urban-Meteorology-Reading/SUEWS-Docs).
* The SUEWS docs are mainly written in [RST](http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html) but [Markdown](https://guides.github.com/features/mastering-markdown/) is also acceptable.
* Your changes to docs will be reviewd and merged into public release if appropriate.

### Questions

* Plesae [raise issues](https://github.com/Urban-Meteorology-Reading/SUEWS/issues/new) for questions in the deveopment so our progress can be well managed.
