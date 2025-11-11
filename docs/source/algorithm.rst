Algorithms
=============



The **uraster** package includes several algorithms to convert existing raster datasets into the unstructured mesh-based datasets.


The package conduct this operation using several steps:

1. Read the user provided configuration, which includes the source mesh, the raster dataset to sample from, and the output target mesh
2. Check the validity of the input datasets. This step includes checking both the source mesh and the raster dataset to ensure they are compatible for sampling. This step also includes the sub-steps:
    - Check the raster datasets, build the raster mesh if necessary.
    - Check the source mesh, and build the mesh topology, which will be used for the visualization later.
3. Sample the raster dataset using the source mesh. The specific sampling method can be configured by the user. The default method is nearest interpolation.
4. Write the sampled data into the target mesh file. The output format can be configured by the user. The default format is GeoParquet.


