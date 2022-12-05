# Energy consumption for Gaussian Filter

## Benchmarking

To run benchmarks:

```sh
# 1. Run make to build project
make all
# 2. Script to run performance benchmarks
./scripts/bench.sh <host> <image>
# 3. Script to run energy benchmarks
./scripts/benchEnergy.sh <host> <image>
```

-------------

## Performance

### i7-9750H

1. Gaussian filter with Golden SSE and AVX versions

![Stats](./stat/i7-9750H/image.i.NGC6744.SB31944.cont.taylor.0.restored.conv.fits-performance.png)

### i7-9750H - Work-station

1. Gaussian filter with Golden SSE and AVX versions

![Stats](./stat/i7-9750H-workstation/image.i.NGC6744.SB31944.cont.taylor.0.restored.conv.fits-performance.png)

### Ryzen-3700X - Work-station

1. Gaussian filter with Golden SSE and AVX versions

![Stats](./stat/Ryzen-3700X/image.i.NGC6744.SB31944.cont.taylor.0.restored.conv.fits-performance.png)

### Apple Silicon: M1 Pro

1. Gaussian filter with Golden SSE and NEON versions

![Stats](./stat/m1pro/image.i.NGC6744.SB31944.cont.taylor.0.restored.conv.fits-performance.png)


## Energy

### i7-9750H-Workstation

1. Gaussian filter with Golden SSE and AVX versions

![Stats](./stat/i7-9750H-workstation/energy/image.i.NGC6744.SB31944.cont.taylor.0.restored.conv.fits-performance.png)

### Ryzen-3700X - Work-station

1. Gaussian filter with Golden SSE and AVX versions

![Stats](./stat/Ryzen-3700X/energy/image.i.NGC6744.SB31944.cont.taylor.0.restored.conv.fits-performance.png)
