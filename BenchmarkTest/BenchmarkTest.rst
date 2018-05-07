this tool is mostly designed for development to check if new develpments can
produce a program (i.e. a new version of SUEWS) meet certain expected behaviours.

current, the following tests will be performed:
  - if a multi-grid run can produce the same results as runs with individual grids?
    so we can know the multi-grid function performs well .

  - if a multi-year run can successfully run?
    so we can know the program run multi-year simulations.

  - if the new release can produce the identical results as the base run?
    so we can know the function of a new release stay the same as the previous one given we are conducting under-the-hood work.
