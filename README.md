# Locally Injective Parametrization with Arbitrary Fixed Boundaries

(Unofficial) Implementation of the paper: Locally Injective Parametrization with Arbitrary Fixed Boundaries

Ofir Weber
Bar Ilan University, Israel
Denis Zorin
New York University, USA

## What does it do

We present an algorithm for mapping a triangle mesh, which is homeomorphic to a disk, to a planar domain with arbitrary fixed boundaries. The algorithm is guaranteed to produce a globally bijective map when the boundary is fixed to a shape that does not self-intersect. 
Meaning: given a model with disk topology, and a target 2D polygon (weakly self-overlapping or non-self-intersect), generate an locally injective/globally bijective map.

Here are some examples using the cut from the dataset (http://vcg.isti.cnr.it/Publications/2014/MPZ14/)

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `param_bin` binary.

## Run

From within the `build` directory just issue:

    ./param_bin --in path_to_data

## Dependencies

The only dependencies are stl, eigen, my fork of [libigl](https://github.com/hankstag/libigl.git) (branch `dev`) and
the dependencies of the `igl::opengl::glfw::Viewer`.

We recommend you to install libigl using git via:

    git clone https://github.com/libigl/libigl.git
    cd libigl/
    git submodule update --init --recursive
    cd ..

If you have installed libigl at `/path/to/libigl/` then a good place to clone
this library is `/path/to/libigl-example-project/`.
