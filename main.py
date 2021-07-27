from Geometry import Geometry
#from Simulation import Simulation
import numpy as np 
from fenics import *
from mshr import *
import matplotlib.pyplot as plt


def intersect_surface(outline_points, point_star, tol):
    x_star, y_star = point_star
    is_on = False
    for i in range(len(outline_points[:, 0]) -1):
        if abs(outline_points[i, 0] - outline_points[i+1, 0]) < tol:
            if abs(x_star - outline_points[i, 0]) < tol and ((y_star <= outline_points[i,1] and y_star >= outline_points[i+1,1]) or (y_star >= outline_points[i,1] and y_star <= outline_points[i+1,1])):
                is_on = True
        elif abs(outline_points[i, 1] - outline_points[i+1, 1]) < tol:
            if abs(y_star - outline_points[i, 1]) < tol and ((x_star <= outline_points[i,0] and x_star >= outline_points[i+1,0]) or (x_star >= outline_points[i,0] and x_star <= outline_points[i+1,0])):
                is_on = True
        elif (x_star <= outline_points[i,0] and x_star >= outline_points[i+1,0]) or (x_star >= outline_points[i,0] and x_star <= outline_points[i+1,0]):
            m = (outline_points[i,1] - outline_points[i+1,1])/(outline_points[i,0] - outline_points[i+1,0])
            b = outline_points[i,1] - m*outline_points[i,0]
            if abs(y_star - m*x_star - b) < tol:
                is_on = True
    return is_on
    
def run_electric_sim(emitter_voltage, collector_voltage, gate_voltage, function_space, trial_fxn, test_fxn, filename):
    a = inner(grad(trial_fxn), grad(test_fxn))*dx
    L = Constant(0.0)*test_fxn*dx
    out = DirichletBC(function_space, Constant(0), out_boundary)
    emitter = DirichletBC(function_space, Constant(emitter_voltage), on_emitter)
    collector = DirichletBC(function_space, Constant(collector_voltage), on_collector)
    gate = DirichletBC(function_space, Constant(gate_voltage), on_gate)
    bcs = [out, emitter, collector, gate]
    A, b = assemble_system(a, L, bcs)
    solver = KrylovSolver('cg', 'ilu')
    INITIAL_GUESS_ZERO = True
    u = Function(function_space)
    U = u.vector()
    solver.solve(A, U, b)
    plt.figure()
    plot(u)
    plt.savefig("test{}.png".format(filename))
    plt.close()
    print(u(Point(-9, 0)))
    print(u(Point(0, 0)))
    print(u(Point(9, 0)))
    return u

if __name__=="__main__":
    x_gap = 20                 # gap between emitter and collector in the x direction 
    y_gap = 50                # gap between 2 gates in the y direction
    w = 100                    # width of emitter and collector
    wg = 80                    # width of gates
    l = 100                     # length of emitter
    lg = 100                   # length of gate
    h = 50                     # height of triangle
    hg = 50                     # height of gate triangle 

    emitter = np.array([[-x_gap/2, 0], [-x_gap/2-h, w/2], [-x_gap/2-h-l, w/2], [-x_gap/2-h-l, -w/2], [-x_gap/2-h, -w/2], [-x_gap/2, 0]])
    collector = np.array([[x_gap/2, 0], [x_gap/2+h, -w/2], [x_gap/2+h+l, -w/2], [x_gap/2+h+l, w/2], [x_gap/2+h, w/2], [x_gap/2, 0]])
    gate = np.array([[0, y_gap/2], [wg/2, y_gap/2+hg], [wg/2,  y_gap/2+hg+lg], [-wg/2,  y_gap/2+hg+lg], [-wg/2, y_gap/2+hg], [0, y_gap/2]])
    geom_object = Geometry(emitter, collector, gate, round_corner=True)
    geom_object.draw_geom()
    Vc = [15]
    Ve = [0]
    Vg = [10]
    tol = 1e-2
    minx = min([np.min(geom_object.emitter_points[:, 0]), np.min(geom_object.collector_points[:, 0]), np.min(geom_object.gate_points[:, 0])]) - 10
    miny = min([np.min(geom_object.emitter_points[:, 1]), np.min(geom_object.collector_points[:, 1]), np.min(geom_object.gate_points[:, 1])]) - 10
    maxx = max([np.max(geom_object.emitter_points[:, 0]), np.max(geom_object.collector_points[:, 0]), np.max(geom_object.gate_points[:, 0])]) + 10
    maxy = max([np.max(geom_object.emitter_points[:, 1]), np.max(geom_object.collector_points[:, 1]), np.max(geom_object.gate_points[:, 1])]) +10
    emitter_coords = [Point(geom_object.emitter_points[i, 0], geom_object.emitter_points[i, 1]) for i in range(len(geom_object.emitter_points))]
    collector_coords = [Point(geom_object.collector_points[i, 0], geom_object.collector_points[i, 1]) for i in range(len(geom_object.collector_points))]
    gate_coords = [Point(geom_object.gate_points[i, 0], geom_object.gate_points[i, 1]) for i in range(len(geom_object.gate_points))]
    emitter = Polygon(emitter_coords)
    collector = Polygon(collector_coords)
    gate = Polygon(gate_coords)
    domain = Rectangle(Point(minx, miny), Point(maxx, maxy)) - emitter - collector - gate
    mesh = generate_mesh(domain, 20)
    function_space =  FunctionSpace(mesh, 'CG', 1)
    trial_fxn = TrialFunction(function_space )
    test_fxn = TestFunction(function_space )
    electric_solution = None
    collector_current = np.zeros((len(Vc), len(Vg)))
    gate_current = np.zeros((len(Vc), len(Vg)))
    num_iter = 0

    def out_boundary(x, on_boundary):
        return on_boundary and (abs(x[0] - minx) <= tol or abs(x[0] - maxx) <= tol or abs(x[1] - miny) <= tol or abs(x[1] - maxy) <= tol) 

    def on_emitter(x, on_boundary):
        return on_boundary and intersect_surface(geom_object.emitter_points, x, tol)

    def on_collector(x, on_boundary):
        return on_boundary and intersect_surface(geom_object.collector_points, x, tol)

    def on_gate(x, on_boundary):
        return on_boundary and intersect_surface(geom_object.gate_points, x, tol)

    xposition = []
    yposition = []

    for Ve_i in Ve:
        for Vc_j in Vc:
            for Vg_k in Vg: 
                U = run_electric_sim(Ve_i, Vc_j, Vg_k, function_space, trial_fxn, test_fxn, "2")

  
 
    #simulation = Simulation(emitter, collector, gate, [0], [5], [10], 4.9)
    #simulation.run_electric_sim(0, 5, 10)
    #simulation.save_solutions("test1")