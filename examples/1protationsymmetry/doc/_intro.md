# Rotation-symmetric pressure distribution

__In this example, you will learn how to__

* solve a rotation-symmetric problem one-dimensionally
* perform a convergence test against an analytical solution
* apply the `Rotational Extrusion` filters in [ParaView](https://www.paraview.org/) for a two-dimensional visualization of the one-dimensional results


__Result__. With the `Rotational Extrusion` and the `Warp By Scalar` filters in [ParaView](https://www.paraview.org/),
the pressure distribution of this example looks as shown in the following picture:

<figure>
    <center>
        <img src="img/result.png" alt="Rotation-symmetric pressure distribution" width="60%"/>
        <figcaption> <b> Fig.1 </b> - Rotation-symmetric pressure distribution on a disc (warped to 3D). </figcaption>
    </center>
</figure>


__Table of contents__. This description is structured as follows:

[[_TOC_]]


## Problem setup

We consider a single-phase problem that leads to a rotation-symmetric pressure distribution.
The following figure illustrates the setup:

<figure>
    <center>
        <img src="img/setup.svg" alt="Rotation-symmetric setup" width="60%"/>
        <figcaption> <b> Fig.2 </b> - Setup for the rotation-symmetric problem. The pressure boundary conditions are shown by the colored lines and the simulation domain is depicted in grey.</figcaption>
    </center>
</figure>

This could, for example, represent a cross section of an injection/extraction well in a homogeneous
and isotropic porous medium, where the well with radius $`r_1`$ is cut out and the
injection/extraction pressure $`p_1`$ is prescribed as a Dirichlet boundary condition. At the outer
radius $`r_2`$, we set the pressure $`p_2`$. In the polar coordinates $`r`$ and $`\varphi`$, the
solution to this problem is independent of the angular coordinate $`\varphi`$ and can be reduced to
a one-dimensional problem in the radial coordinate $`r`$. Therefore, in this example, we want to
solve the problem on a one-dimensional computational domain as illustrated by the orange line in
the above figure.

## Mathematical model

In this example we are using the single-phase model of DuMuX, which considers Darcy's law to relate
the Darcy velocity $`\textbf u`$ to gradients of the pressure $`p`$. In the case of rotational
symmetry, the mass balance equation for the fluid phase can be transformed using polar coordinates:

```math
-\frac{1}{r} \frac{\partial}{\partial r} \left( r  \frac{\varrho k}{\mu} \frac{\partial p}{\partial r} \right) = 0.
```
Here, $`k`$ is the permeability of the porous medium, $`\mu`$ is the dynamic viscosity of the
fluid, $`\phi`$ is the porosity, and $`\varrho`$ is the fluid density.

## Discretization

We employ a finite-volume scheme to spatially discretize the mass balance equation shown above.
The discrete equation describing mass conservation inside a control volume $`K`$ is obtained
by integration and reads:

```math
    \sum_{\sigma \in \mathcal{S}_K} | \sigma | \left( \varrho u \right)_\sigma
    = 0,
```

where $`\sigma`$ are the faces of the control volume such that
$`\bigcup_{\sigma \in \mathcal{S}_K} \sigma \equiv \partial K`$ and where the notation
$`( \cdot )_\sigma`$ was used to denote quantities evaluated for a
face $`\sigma`$. The area of a face is denoted with $`| \sigma |`$.

DuMuX provides the classes `RotationSymmetricSubControlVolume` and
`RotationSymmetricSubControlVolumeFace`, which implement one-dimensional control
volumes and faces, that, in the computations of volumes and areas, take into account
the extrusion about the rotation axes of symmetry. This will be discussed in part 1
of the documentation.

# Implementation & Post processing
