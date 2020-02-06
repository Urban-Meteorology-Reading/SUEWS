---
name: release checklist
about: Checklist for regular release. ONLY USEFUL FOR DEVELOPERS!
title: 'Checklist for V202Yx'
labels: docs
assignees: ''

---
## New Features
<!--planned new features-->

## Fixes
<!--planned fixes-->

## Code Health
- [ ] update release identifier: [`progname` in `suews_ctrl_const.f95`](../../SUEWS-SourceCode/suews_ctrl_const.f95)
- [ ] update test identifier: [`name_exe` in `BTS_config.nml`](../../Test/code/BTS_config.nml)
- [ ] run simple test:
  - [ ] set up test folder for the new release: `Test/BaseRun/$release_ver`
  - [ ] use `make test` in the root folder for testing

**Docs**
- [ ] update [scheme descriptions](../../docs/source/parameterisations-and-sub-models.rst) if any.
- [ ] update [option descriptions](../../docs/source/input_files/RunControl/RunControl.rst) and relevant csv files if any.


**Additional Remarks**
<!-- Other remarks about this release -->
