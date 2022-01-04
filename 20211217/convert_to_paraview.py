from dolfin import *

mesh = Mesh("coil.xml")
domains = MeshFunction("size_t", mesh, "coil_physical_region.xml")
File("domains.pvd") << domains
