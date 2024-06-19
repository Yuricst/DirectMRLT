# DirectMRLT

Many revolution low-thrust (MRLT) trajectory optimization with direct collocation using ICLOCS on MATLAB


For initializing the repository, run

```bash
git clone --recurse-submodules https://github.com/Yuricst/DirectMRLT.git
git submodule update --init --recursive  # if --recurse-submodules flag was forgotten
```

### Submodule dependencies:

- [ICLOCS](https://github.com/ImperialCollegeLondon/ICLOCS)

### Toolboxes:

- **Optimization** : needed if using fmincon
- **Signal Processing Toolbox** : needed if using hp-grids (recommended)


### Setting up IPOPT

- Windows: [OPTI](https://github.com/jonathancurrie/OPTI)
- OSX: [mexIPOPT](https://github.com/ebertolazzi/mexIPOPT)


## Dev notes

- [ ] MEE dynamics
    - [x] Two-body
    - [ ] Spherical harmonics
    - [ ] SRP
- [ ] Eclipse path constraints