# Velocity_Fields

Compile

```gcc calculate_velocity_fields.c -lm -O3 -o make_fields```

The code currently anticipates files of structure:
```velocity_data/XV_data_t=00003800```
but could easily be embedded in the molecular dynamics code itself to operate more efficiently on future runs.
The code generates field data organized by time
```velocity_data/field_00003800```
so each file contains a uniform x-y grid (2D marked out by j,k)
that has an estimate of the velocity of the particles at that spot, a measure of the system curl, and enstropy, in the format:

``` x_grid, y_grid, vx_field[j][k],  vy_field[j][k], curl_field[j][k], enstrophy_field[j][k]```

once this has run, you can apply the python file to make some plots.  modify these as you see fit.

```python plot_fields.py```

where the code expects files of time
```velocity_data/field_00003800```
```velocity_data/XV_data_t=00003800```

and is currently hardcoded to plot time ```9950```
