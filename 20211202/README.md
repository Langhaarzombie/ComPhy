# First meshing exercise

This should briefly explain how to reproduce the screenshots and the struggles I had on the way.

## Installation

### Programs used

For this exercise we were using *gmsh* for creating a mesh on a geometry, *FEniCS* for translating the mesh network into
a .xml format, *python* to convert that .xml to a .pvd file and finally *Paraview* to look at the geometry and the mesh
we have created.

#### Installing gmsh

You can technically just download the sources from their website.
For Arch based systems, you can download the application from the AUR repository.
There are two versions though and you should pick the *gmsh-bin* package as the other one will just build the sources and throw an error after about an hour.

This is all you need to run gmsh. At point of writing I was using version 4.8.4.

#### Installing FEniCS

You can either use the provided VirtualBox (from the university) or go the quicker way and use the Docker image maintained by the FEniCS people
[here](https://fenics.readthedocs.io/projects/containers/en/latest/).
This website should clarify everything you need to run the container.
No problems there for me.

#### Installing Paraview

Technically there is also a Paraview Container (web based).
That one however only works with Nvidia graphics cards and I had too little time to find out if there was a serious option for Intel or AMD ones.

Again for Arch based systems you might be tempted to download the applciation from the community or the AUR repository.
This however leads to a conflict with *vtk* (or so) for me.
*vtk* is used by OpenCascade which again is used by gmsh (I think).
So you might be able to resolve that conflict by running *gmsh* in a container as well.

In the end I downloaded the binary from their website and used just that.

## Doing the exercise

### gmsh

```
Create mesh.geo and check for the correct representation of your geometry in gmsh.
Then create a 3D mesh.
```

So first things first, the syntax for gmsh is quite horrible.
However, we got a pretty good tutorial during the lecture and that is really all you need.
The only problem that occured on my machine was a *Segmentation faul*.
I don't know why.
Simply running the command again solved it (for some time).

In the end, one error that also occured during the lecture, happend with my script.
I had the problem that my Volume 100 would not have any elements.
After trying and googling everything, I found an [issue](https://gitlab.onelab.info/gmsh/gmsh/-/issues/686) that helped me out.
Deleting the volumes (as done in the *mesh.geo* file) made the script run error free in the end.

### FEniCS

```
Transform the .msh file to .xml.
```

Open a bash line in the running Docker container and execute the given command.

### Pyhton

```
Transform .xml to .pvd.
```

Still in the Docker container (as you need the dolfin dependency that comes with FEniCS), you run the script that was provided in the PDF.
Note however, that you do not need to stuff with the facets.

### Paraview

```
Create the screenshots and check if everything went well.
```

As discussed in the lecture you can make some cuts, to see if everything worked.
If there is anything missing (like for me the coil due to the error described above),
go back to the .geo script and check for errors.

