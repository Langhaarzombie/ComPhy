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
uC1 = Constant(0)
uC2 = Constant(1)
L  = f*v*dx + uC1*v*ds(2) + uC2*v*ds(3)
# Dirichlet boundary conditions
#  bcs = dirichlet(contact1, contact2)
#  bcs = dirichlet_neumann(contact1)

# Assemble
#  A, b = assemble_system(a, L, bcs)
A, b = assemble_system(a, L)

# Solve system
u = Function(V)
solve(A, u.vector(), b)

VV = VectorFunctionSpace(mesh, "CG", 1)
j = project(grad(u), VV)

File("j.pvd") << j
