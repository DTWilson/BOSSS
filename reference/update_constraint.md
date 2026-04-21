# Update problem constraints

Update problem constraints

## Usage

``` r
update_constraint(
  problem,
  number = NULL,
  name = NULL,
  nom = NULL,
  delta = NULL
)
```

## Arguments

- problem:

  BOSSS problem

- number:

  index of the constraint(s) to be updated

- name:

  name of the constraint(s) to be updated

- nom:

  new nominal value(s)

- delta:

  new acceptance threshold(s)

## Value

an updated BOSSS problem object.
