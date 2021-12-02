from dolfin import *

mesh = Mesh("mesh.xml")
domains = MeshFunction("size_t", mesh, "mesh_physical_region.xml")
File("domains.pvd") << domains
