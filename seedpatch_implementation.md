# Seed-patch implementation: stage 1

The existing Sir-style sequence remains Compute_Field_Bounds, Calculate_dPhi,
Compute_Wedges, Filter_Hits, followed by the new Form_Seed_Patches function.
Sir's source file itself is not changed by this implementation.

## Input and units

Default input: hits_truth1000.csv. An optional first command-line argument selects
another CSV. Columns are located by name. Source coordinates must be millimetres;
Filter_Hits divides x/y/z by 10 before applying the centimetre luminous bounds.
Momentum and angles are not scaled. Filtered phi remains in degrees.
The existing geometry is retained, including outer z bound 49.078 cm (the later
discussion's 49.708 example did not change that configuration).

Event/event_id, when present, separates collisions. Without it, the file is
explicitly treated as ONE event labelled 0; do not supply a merged multi-event
file without event identifiers. Particle IDs are output metadata only.
Invalid input stops the run with a row error. Outputs from a failed run may be
partial and must not be treated as complete.

## First-patch selection

1. Store filtered hits by event, wedge, and physical layer (1 through 4).
2. Sort each list by z, with hit_id as a deterministic tie-breaker.
3. Select the highest 16 consecutive hits in each endpoint layer, 1 and 4.
4. For layers 2 and 3, interpolate endpoint minimum and maximum coordinates:
   alpha = (Radius[layer]-Radius[1])/(Radius[4]-Radius[1]).
   z = (1-alpha)*z1 + alpha*z4.
5. Try a 16-hit consecutive window whose upper endpoint is the first available
   hit at or above the required upper bound (or index 15 if more preceding hits
   are needed). Accept it as covering only if its lower bound covers the required
   minimum. If no such covering window exists, select the 16 hits ending at the
   last hit at or below the required upper bound and clip the patch accordingly.
6. Start with the rectangle defined by the endpoint superpoints in (z1,z4).
   Clip against both intermediate-layer intervals using Clip_Patch.
7. Skip groups with insufficient selectable hits or zero-area intersections.

The polygon describes the shared straight-line parameter region. Its 64 stored
hits are candidates, not a claim that they all belong to one particle or one line.
This stage does not additionally clip by the extrapolated beam intercept.

## Functions and output

- Filter_Hits: header parsing, unit conversion, filtering, and in-memory grouping.
- Select_Superpoint: a 16-hit window ending at/below a target z.
- Clip_Patch: intersection of a convex polygon with one linear inequality.
- Form_Seed_Patches: first-patch selection, geometry, and CSV output.

seedpatch_hits.csv contains 64 rows per accepted patch, hit metadata, centimetre
coordinates, and superpoint bounds. seedpatch_corners.csv contains one row per
polygon vertex, with z1_cm and z4_cm. Patch IDs start at zero.

Only ONE first patch per event/wedge is implemented. Advancing through subsequent
patches, full coverage, and complementary patches remain for the next discussion.

## Build and run

From a Visual Studio developer terminal:

    cl /std:c++17 /EHsc /W4 /O2 seedpatch.cpp /Fe:seedpatch.exe
    .\seedpatch.exe

Or with GCC:

    g++ -std=c++17 -O2 -Wall -Wextra seedpatch.cpp -o seedpatch.exe
    .\seedpatch.exe

Outputs are written to the current working directory. Rebuild before running;
an executable already in the project directory may belong to an older source.

## Validation

MSVC C++17 /W4 build succeeded. The supplied dataset produced 128 first patches,
512 superpoints, 8192 hit rows, and 671 polygon vertices. Every superpoint had
16 hits, and every exported vertex satisfied all four layer constraints within
1e-9 cm. Validation outputs and executable were kept in a temporary directory.
