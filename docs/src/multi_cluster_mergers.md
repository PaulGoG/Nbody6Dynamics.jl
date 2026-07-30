# Multi-Cluster Merger Simulations

## Background

Cluster merger simulations model the gravitational encounter and coalescence of two or more star clusters. This is relevant to:

- **Nuclear star cluster formation** via inspiral of globular clusters toward a galactic centre
- **Hierarchical cluster assembly** in young massive cluster complexes (e.g. R136, Westerlund 1)
- **Merger remnant identification** (e.g. NGC 1851, Terzan 5 -- suspected merger products with multiple stellar populations)
- **Intermediate-mass black hole (IMBH) formation** via runaway collisions during core mergers

## Physics of Cluster Mergers

### Orbital parameters

Two clusters on a Keplerian orbit are characterised by:

| Parameter | Symbol | Meaning |
|:----------|:-------|:--------|
| Apocentre distance | ``d_{\rm apo}`` | Maximum separation |
| Eccentricity | ``e`` | Orbital shape (0 = circular, 1 = parabolic) |
| Semi-major axis | ``a = d_{\rm apo}/(1+e)`` | Orbit size |
| Mass ratio | ``q = M_2/M_1 \leq 1`` | Relative cluster masses |

The **mutual escape velocity** at separation ``d`` is:

```math
v_{\rm esc} = \sqrt{\frac{2G(M_1 + M_2)}{d}}
```

- ``v_{\rm rel} < v_{\rm esc}``: bound orbit, merger inevitable (sub-virial)
- ``v_{\rm rel} = v_{\rm esc}``: parabolic, merger likely for close passages
- ``v_{\rm rel} > v_{\rm esc}``: hyperbolic fly-by, merger only via strong dissipation

### Typical physical scales

| Quantity | Typical range | Notes |
|:---------|:-------------|:------|
| Initial separation | 5--20 ``r_{\rm hm}`` | Must exceed tidal radii of both clusters |
| Relative velocity | 1--10 km/s | Sub-virial to mildly parabolic in young complexes |
| Merger timescale (``q \sim 1``) | 50--200 Myr | 2--5 orbital periods from first pericentre |
| Merger timescale (``q \sim 0.1``) | 500 Myr -- several Gyr | Dynamical friction dominated |
| "Independent" separation | ``> 10\, r_{\rm hm}`` | Tidal perturbation negligible |

### Density profiles

| Model | Properties | Use case |
|:------|:-----------|:---------|
| **King** (recommended) | Finite tidal radius, flat core, concentration ``c = \log_{10}(r_t/r_c)`` | Realistic GC mergers; natural truncation prevents density overlap |
| **Plummer** | Infinite extent (truncated in practice), single scale radius | Quick parameter surveys, analytically tractable |
| **Wilson** | Like King but with different outer profile | Alternative when King truncation is too sharp |
| **limepy family** | Generalises King (``g=1``), Wilson (``g=2``), Woolley (``g=0``) | Gold standard for modern IC generation |

**King models are preferred** for merger simulations because their finite tidal radius prevents artificial density overlap between subclusters at the initial separation.

### Known issues with merger ICs

**1. Virial ratio**: Two independently virialised clusters (``Q = 0.5`` each) are NOT in virial equilibrium as a combined system. The mutual gravitational potential energy must be accounted for. The Nbody6++ `SCALE` subroutine must NOT rescale velocities for the combined system.

**Solution**: Use `KZ(22)=2` (external particle input), which sets `LSCALE=.FALSE.` -- no centre-of-mass correction or velocity rescaling is applied.

**2. Tidal truncation at initial separation**: Each cluster feels the tidal field of the other. Stars beyond the instantaneous Jacobi radius:

```math
r_J \approx d \left(\frac{M_{\rm self}}{3\,M_{\rm other}}\right)^{1/3}
```

will be immediately unbound. Best practice: truncate each King model at ``\min(r_t, r_J)`` to avoid transient mass loss and energy errors.

**3. Mass segregation**: Primordial mass segregation in each subcluster is disrupted during the merger. The merged remnant re-establishes mass segregation on its own (longer) half-mass relaxation timescale.

## Nbody6++ Support

### Built-in: Two Plummers (`KZ(5)=2`)

Nbody6++GPU has native support for binary Plummer encounters via `KZ(5)=2` in `setup.F`. Parameters are read from the `&INSETUP` block:

| Parameter | Meaning | Constraints |
|:----------|:--------|:------------|
| `APO` | Apocentre distance [NB units] | Semi = APO/(1+ECC), clamped to [2, 50] |
| `ECC` | Eccentricity | [0, 0.999] |
| `N2` | Particle count in second cluster | ``\leq N`` |
| `SCALE` | Size ratio of second cluster | [0.2, 5.0] |

**Limitations**:
- Only 2 clusters
- Both are Plummer profiles (no King model option)
- Second cluster is subsampled from the first (same mass function realisation)
- No independent density profiles or concentrations

### External particle input (`KZ(22)=2`)

The recommended approach for custom merger ICs. Write a `dat.10` file with one line per particle:

```
MASS  X  Y  Z  VX  VY  VZ
```

All values in N-body units (``G=1``, ``M_{\rm total}=1``). Set `KZ(22)=2` in the `.inp` file. With `KZ(22)=2`:
- `LSCALE=.FALSE.`: no velocity rescaling, no CM correction
- Masses are normalised to sum to 1
- Positions and velocities are used as-is

For astrophysical units (``M_\odot``, pc, km/s), use `KZ(22)=10` instead.

### Tidal field options (`KZ(14)`)

| `KZ(14)` | Field type | Parameters |
|:----------|:-----------|:-----------|
| 0 | Isolated (no tidal field) | -- |
| 1 | Oort constants (solar neighbourhood) | built-in |
| 2 | Point-mass galaxy | `GMG` [``M_\odot``], `RG0` [kpc] |
| 3 | Disk + halo | `GMG`, `DISK`, `VCIRC`, etc. |
| 5 | Milky Way potential (Bovy 2015) | `RG(1:3)`, `VG(1:3)` |

## Proposed Implementation in Nbody6Setup.jl

### New module: `InitialConditions`

A Julia-native initial condition generator that produces `dat.10` files for Nbody6++ with `KZ(22)=2`. This fills a gap -- no existing tool provides a turnkey multi-cluster merger IC generator.

### Cluster models to implement

1. **Plummer model**: Inverse CDF sampling (analytic). Single parameter: scale radius ``a``.
2. **King model**: Eddington inversion of the lowered isothermal distribution function. Parameters: central potential ``W_0`` (or concentration ``c``), tidal radius ``r_t``.
3. **IMF sampling**: Kroupa (2001) broken power law as default, with configurable slopes and mass limits.

### Orbital setup

Given ``N_{\rm clusters} \geq 2`` clusters with masses ``M_i``, positions ``\mathbf{r}_i``, and velocities ``\mathbf{v}_i``:

1. Generate each cluster independently in its own centre-of-mass frame
2. Optionally apply Jacobi truncation at the mutual tidal radius
3. Compute two-body orbital velocities from Kepler's equation (for 2-cluster case)
4. For ``N > 2`` clusters: user specifies positions and velocities directly, or places them on a hierarchical (nested two-body) orbital configuration
5. Apply CM offsets to each cluster's particles
6. Combine into a single particle array, normalise to N-body units
7. Write `dat.10` + generate matching `.inp` file with `KZ(22)=2`

### Configuration

A new TOML section `[merger]` or a standalone merger config:

```toml
[merger]
n_clusters = 2

[merger.cluster1]
model = "king"          # "king", "plummer"
N = 50000
W0 = 6.0               # King concentration (W0), ignored for Plummer
mass_total = 1e5        # Solar masses
rbar = 2.0              # Half-mass radius [pc]
imf = "kroupa"
body1 = 100.0           # Upper mass limit
bodyn = 0.08            # Lower mass limit
metallicity = 0.001

[merger.cluster2]
model = "king"
N = 50000
W0 = 4.0
mass_total = 5e4
rbar = 3.0
imf = "kroupa"
body1 = 100.0
bodyn = 0.08
metallicity = 0.001

[merger.orbit]
apocentre = 15.0        # pc (physical units)
eccentricity = 0.7
# OR for N>2 clusters, specify positions/velocities directly:
# positions = [[0,0,0], [15,0,0], [-10,5,0]]
# velocities = [[0,0,0], [0,-2,0], [0,1,0]]  # km/s

[merger.output]
format = "nbody"        # "nbody" (KZ22=2) or "astro" (KZ22=10)
truncate_jacobi = true  # Truncate at mutual Jacobi radius
```

### Output

1. `dat.10` -- particle data file for Nbody6++
2. `.inp` -- matching input file with `KZ(22)=2`, correct `N`, `RBAR`, `ZMBAR`
3. Summary log with physical parameters, expected merger timescale, virial ratios of each subcluster

### Diagnostics and post-processing extensions

The existing post-processing pipeline handles merger simulations without modification -- energy, snapshots, HR diagrams, and animations all work on the combined particle set. However, useful additions would include:

- **Subcluster identification**: Track which particles belong to which original cluster (via particle NAME ranges) and plot their spatial separation over time
- **Merger detection**: Identify the merger epoch as the time when the two density peaks merge into one
- **Lagrangian radii per subcluster**: Track half-mass radii of each original cluster component separately

## References

- Aarseth, S.J. (2003). *Gravitational N-Body Simulations*. Cambridge University Press.
- Fujii, M.S. et al. (2012). "The formation of young dense star clusters through mergers." ApJ, 753, 85.
- Gieles, M. & Zocchi, A. (2015). "A family of lowered isothermal models." MNRAS, 454, 576.
- Kuepper, A.H.W. et al. (2011). "McLuster -- A tool to make star clusters." MNRAS, 417, 2300.
- Livernois, A.R. et al. (2022). "Modelling star cluster formation: Mergers." MNRAS, 513, 6095.
- de Oliveira, M.R. et al. (2000). "Final stages of N-body star cluster encounters." MNRAS, 311, 589.
- Kroupa, P. (2001). "On the variation of the initial mass function." MNRAS, 322, 231.
- King, I.R. (1966). "The structure of star clusters. III." AJ, 71, 64.
- Bovy, J. (2015). "galpy: A Python library for galactic dynamics." ApJS, 216, 29.
