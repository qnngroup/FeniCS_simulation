from dolfin import *
from numpy import array

# Simple mesh/space for testing:
mesh = UnitIntervalMesh(3)
tree = mesh.bounding_box_tree()
V = FunctionSpace(mesh,"CG",1)

# Evaluate the i-th basis function at point x:
def basis_func(i,x):
    cell_index = tree.compute_first_entity_collision(Point(*x))
    cell_global_dofs = V.dofmap().cell_dofs(cell_index)
    for local_dof in range(0,len(cell_global_dofs)):
        if(i==cell_global_dofs[local_dof]):
            cell = Cell(mesh, cell_index)
            values = array([0,])
            return V.element().evaluate_basis(local_dof,x,
                                              cell.get_vertex_coordinates(),
                                              cell.orientation())[0]
    # If none of this cell's shape functions map to the i-th basis function,
    # then the i-th basis function is zero at x.
    return 0.0

# Build a Python list of all basis functions, as requested:
list_of_basis_funcs = []
for i in range(0,V.dim()):
    list_of_basis_funcs += [(lambda i: lambda x : basis_func(i,x))(i),]

# Try evaluating all basis functions at a point:
for i in range(0,V.dim()):
    print(list_of_basis_funcs[i](array([1/3+1/12,])))
