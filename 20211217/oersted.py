from dolfin import *

class Box(SubDomain):
    def inside(self, x, on_boundary):
        return near(abs(x[0]), 150) or near(abs(x[1]), 150) or near(abs(x[2]), 150)

mesh = Mesh("coil.xml")

volume_domains = MeshFunction("size_t", mesh, "coil_physical_region.xml")
dx = Measure("dx", mesh, subdomain_data=volume_domains)

facet_domains = MeshFunction("size_t", mesh, "coil_facet_region.xml")
facet_domains.set_all(0)

box = Box()
box.mark(facet_domains, 1)

ds = Measure("ds", mesh, subdomain_data=facet_domains)

V = FunctionSpace(mesh, "CG", 1)

H_x = TrialFunction(V)
H_y = TrialFunction(V)
H_z = TrialFunction(V)
v = TestFunction(V)

# LHS
a_x = dot(grad(H_x), grad(v))*dx
a_y = dot(grad(H_y), grad(v))*dx
a_z = dot(grad(H_z), grad(v))*dx

# RHS
j = Expression(("-x[1]", "x[0]", "0.0"), degree=1, domain=mesh)
cj = curl(j)
L_x = dot(cj[0], v)*dx(5)
L_y = dot(cj[1], v)*dx(5)
L_z = dot(cj[2], v)*dx(5)

# Assemble System
A_x, b_x = assemble_system(a_x, L_x, DirichletBC(V, Constant(0), box))
A_y, b_y = assemble_system(a_y, L_y, DirichletBC(V, Constant(0), box))
A_z, b_z = assemble_system(a_z, L_z, DirichletBC(V, Constant(0), box))

# Solve it
H_x = Function(V)
H_y = Function(V)
H_z = Function(V)
solve(A_x, H_x.vector(), b_x)
solve(A_y, H_y.vector(), b_y)
solve(A_z, H_z.vector(), b_z)

H = Expression(("A", "B", "C"), A=H_x, B=H_y, C=H_z, degree=1, domain=mesh)

VV = VectorFunctionSpace(mesh, "CG", 1)
H = project(H, VV)

File("oersted.pvd") << H

