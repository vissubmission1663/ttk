## Supplementary code for IEEE VIS 2023 submission #1663

This repository contains the implementation of the methods described in
VIS submission 1663 "ExTreeM: Scalable Augmented Merge Tree Computation
via Extremum Graphs". 

The implementation is based on [TTK - The Topology
Toolkit](topology-tool-kit.github.io) and hence this repository is a
fork of the [TTK repository](github.com/topology-tool-kit). The
additions to TTK we make towards our implementation can be easily viewed
by
[comparing](https://github.com/topology-tool-kit/ttk/compare/dev...vissubmission1663:ttk:dev)
the two repositories.

Our core algorithms are implemented in TTK's base layer and as VTK
filters, available through a Paraview plugin. We provide a Paraview
Python script to reproduce a figure from the paper. The build is
configured such that only the minmal subset of TTK required to run the
test case is compiled, resulting in a much smaller-than-typical build
overhead (approx. 170 files).

### Running the benchmarks

As a fork of TTK, our implementation can be compiled following TTK's
installation instructions. We provide several easier options, however,
to allow reproduction of the test case.

#### Visual Studio Code – Development Container

The repository contains a `.devcontainer` description from which [Visual
Studio Code](code.visualstudio.com) can directly generate a pre-made
compilation environment running in a container:

- Install Docker Desktop and ensure the "Remote Development" extension for Visual Studio Code is installed
- Clone the repository and open the top-level directory in Visual Studio Code
- Run the Visual Studio command "Reopen in Container"

Once the folder is reopened in the container, the code can be configured
and compiled by selecting the "Benchmark" CMake preset and then clicking
"Build" (or pressing F7).

To run the test case, press F5, or open a Terminal in the Visual Code
shell, and use 
```
./run-benchmark.sh
```
to run the test case. This will output a single image into the `output`
folder, which can be viewed by clicking on it in the sidebar on the left.

To vary the number of threads, the `OMP_NUM_THREADS` environment
variable can be used, e.g.
```
OMP_NUM_THREADS=4 ./run-benchmark.sh
```

Notes:
- The build process should work identically in a local (i.e.,
  non-container) environment that can successfully build Paraview and TTK.

- It is feasible to open the repository to inspect and compile the code
  in a **Github Codespace**. To do so, navigate to the repository start
  page on Github and press '.', or activate a codespace from the "Code"
  button. Due to memory limitations, it is unfortunately not possible to
  run the test case in the codespace.

#### Build and run Docker container

We provide a top-level [`Dockerfile`](Dockerfile) to create a Docker
image that contains a compiled and runnable version of our
implementation. It is based on the publically available [TTK Docker images](https://hub.docker.com/u/topologytoolkit).

To build the container, run 
```
docker build -t vis1663 .
```
from the repository top-level directory. The resulting image runs the
test case when started:
```
docker run -it --rm vis1663 
```
The image output is written into the container's `/output` directory and
can thus be made accessible outside using e.g.
```
docker run -it --rm -v $(pwd)/output:/output vis1663
```

To vary the number of threads, the `OMP_NUM_THREADS` environment
variable can be set via
```
docker run -it --rm -e OMP_NUM_THREADS=4 vis1663
```

#### Singularity Image

For environments where Docker cannot be used, e.g. on supercomputers, we
provide a [Singularity](https://sylabs.io/singularity/) definition file
equivalent to the Dockerfile described above.

To create the Singularity image, run
```
sudo singularity build ../vis1663.sif ttk.def
```
This will compile and install the implementation. The test case can then
be executed using the command 
```
singularity run ../vis1663.sif
```
or, alternatively, using the longer command 
```
singularity run -B /tmp/output:/output --env OMP_NUM_THREADS=8 ../vis1663.sif
```
that also writes the output image `test.png` to `/tmp` and limits the
number of threads using the `OMP_NUM_THREADS` environment variable.

Note: Due to restrictions on many HPC systems, it may be necessary to
build the image file vis1663.sif on a host where sudo rights are
available. The image file can then be transported to the HPC environment
and run there.