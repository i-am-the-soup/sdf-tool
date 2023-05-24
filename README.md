# sdf-tool
This little doodad takes in black-or-white image sequences and smoothly blends them together through some SDF magic.
# Building
Currently requires a relatively new development zig build (tested on Linux with zig version `0.11.0-dev.2969+855493bb8`)
No fuss needed, just run `zig build -Doptimize=ReleaseFast`.
# Usage
```sdf-tool <folder with input images> <path to output image> [number of blending steps (defaults to 20)]```

e.g.
`sdf-tool test-input/ sdf.png 50`
# Caveats
Input folder may only contain greyscale image files of matching dimensions.
The threading used is very rudimentary at the moment, with hard-coded thread counts that might not be optimal for other machines.
