from dolfin import *

# Read mesh
mesh = Mesh("bone.xml")

# Define facet domains
class Contact1(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], 0.)

class Contact2(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], 10.)

# Functions for different boundary conditions
def dirichlet(c1, c2):
  u0 = Constant(0)
  uL = Constant(1)
  return [DirichletBC(V, u0, contact1), DirichletBC(V, uL, contact2)]

def dirichlet_neumann(c1):
  u0 = Constant(0)
  return [DirichletBC(V, u0, contact1)]

facet_domains = MeshFunction("size_t", mesh, 2)
facet_domains.set_all(0)

boundary = DomainBoundary()
boundary.mark(facet_domains, 1)

contact1 = Contact1()
contact1.mark(facet_domains, 2)

contact2 = Contact2()
contact2.mark(facet_domains, 3)

# Create measures and function spaces
# now you can integrate e.g. over contact 1 with "ds(2)"
ds = Measure('ds', mesh, subdomain_data=facet_domains)

V  = FunctionSpace(mesh, "CG", 1)

# Create coefficients
u = TrialFunction(V)
v = TestFunction(V)

rho = Constant(1.0)

# Define forms and boundary conditions
# Bilinear form
a  = rho*dot(grad(u), grad(v))*dx
# Linear form
f = Constant(0)
L_DD  = f*v*dx
L_DN  = f*v*dx + Constant(1.0)*v*ds(3)
L_NN  = f*v*dx + Constant(0.0)*v*ds(2) + Constant(1.0)*v*ds(3)

# Assemble
A_DD, b_DD = assemble_system(a, L_DD, dirichlet(contact1, contact2))
A_DN, b_DN = assemble_system(a, L_DN, dirichlet_neumann(contact1))
A_NN, b_NN = assemble_system(a, L_NN)

# Solve system
u_DD = Function(V)
u_DN = Function(V)
u_NN = Function(V)
solve(A_DD, u_DD.vector(), b_DD)
solve(A_DN, u_DN.vector(), b_DN)
solve(A_NN, u_NN.vector(), b_NN)

VV = VectorFunctionSpace(mesh, "CG", 1)
j_DD = project(grad(u_DD), VV)
j_DN = project(grad(u_DN), VV)
j_NN = project(grad(u_NN), VV)

File("j_DD.pvd") << j_DD
File("j_DN.pvd") << j_DN
File("j_NN.pvd") << j_NN
