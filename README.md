# MOS Capacitor Simulator
For Numerical Methods for Device Modelling course.

This simulator uses Scharfetter-Gummel Discretization to simulate a MOS Capacitor and output capacitance and charge accumulated at each bias. This is still a purely academic project meant for learning numerical methods and for qualitative clarity.

## How to run
- `make all`
- `./build/moscv <config file> <output file>`
- plot column 1 and 6. for V vs Cap.

## Note
- If cofig file does not exist it will be created.
- Not all parameters can be changed through config file as of now. timesteps, bias points, etc are hardcoded in program. Only material parameters are in config file.


