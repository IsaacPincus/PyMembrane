#import the code
import pymembrane as mb
import os
import numpy as np
from pprint import pprint
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scienceplots
plt.style.use(['science'])

# ### Parse arguments
# ## Now we want to have X snapshots every X steps each
# parser = argparse.ArgumentParser(description="Please provide: snapshots and run_steps")
# ## Add arguments for snapshots and run_steps
# parser.add_argument("--snapshots", type=int, required=True, help="Number of snapshots")
# parser.add_argument("--max_iter", type=int, required=True, help="Number of iteration steps")

# user_args = parser.parse_args()
# # Access the parsed arguments
# snapshots = user_args.snapshots
# max_iter = user_args.max_iter

os.chdir('/home/ipincus/fork_pymembrane/PyMembrane/docs/examples/07_particle/BD/')

## Now we want to have x snapshots every y steps each
snapshots = 50
run_steps = 3000
# snapshots = 5
# run_steps = 500

vertex_file = 'minimized_vertices_RBC_shape.dat'
face_file = 'minimized_faces_RBC_shape.dat'
# vertex_file = 'minimized_vertices.dat'
# face_file = 'minimized_faces.dat'
# vertex_file = 'vertices_sphere_N1024.inp'
# face_file = 'faces_sphere_N1024.inp'
# vertex_file = 'vertices_ellipse.inp'
# face_file = 'faces_ellipse.inp'

#create a system 
box = mb.Box(40,40,40)
system = mb.System(box)
system.read_mesh_from_files(files={'vertices':vertex_file, 'faces':face_file})
compute = system.compute

#save the mesh to display
#create dumper
dump = system.dumper
dump.vtk("initial mesh", False)

vertices = system.vertices
current_max = 0.0
for vertex in vertices:
    z_val = vertex.r.z
    if z_val>current_max:
        current_max = z_val
print("current max is: {}".format(current_max))

vertices = system.vertices
p1 = np.array([-0.41, 1, 0.66])
p2 = np.array([0.30, -1.12, -0.51])
dist = 0.4
for vertex in vertices:
    x_val = vertex.r.x
    y_val = vertex.r.y
    z_val = vertex.r.z
    p = np.array([x_val, y_val, z_val])
    r_p1_p = np.linalg.norm(p1-p)
    r_p2_p = np.linalg.norm(p2-p)
    # print(r_p2_p)
    if np.isclose(r_p1_p, 0.0, rtol=0, atol=dist):
        vertex.type = 5
    if np.isclose(r_p2_p, 0.0, rtol=0, atol=dist):
        vertex.type = 6
    # print(vertex.id, vertex.type)

system.vertices = vertices
for v in system.vertices:
    # print(v.id, v.type)
    if v.type != 1:
        print(v.id, v.type)

#add the evolver class where the potentials and integrators are added
evolver = mb.Evolver(system)

# tether force on vertices
Fm = 1.0
u_p1p2 = (p2-p1)/np.linalg.norm(p2-p1)
evolver.add_force("Mesh>Tether", {
    "Fx": {"5": str(-1.0*Fm*u_p1p2[0]), "6":str(1.0*Fm*u_p1p2[0])},
    "Fy": {"5": str(-1.0*Fm*u_p1p2[1]), "6":str(1.0*Fm*u_p1p2[1])},
    "Fz": {"5": str(-1.0*Fm*u_p1p2[2]), "6":str(1.0*Fm*u_p1p2[2])}
})

# add particle and force
system.add_particle(-0.35,-0.4,0.5, 0.46)
# system.add_particle(0.0,2.0,0.0, 0.5)
# system.add_particle(2.0,0.0,0.0, 0.5)
evolver.add_force("Mesh>Particle", {
    "epsilon": {"1": str(1.0)},
    "sigma": {"1": str(0.02)},
    "phi": {"1": str(1.0)}
})

## Compute the initial energy
energy = compute.energy(evolver)
print("[Initial] energy = ", energy)

# # first we need to know the edge length to move it appropriate:
# edge_lengths = compute.edge_lengths()
# avg_edge_length= np.mean(edge_lengths)
# print("[Initial] avg_edge_length = ", avg_edge_length)
# k = str(5.0)
# l0 = str(avg_edge_length)
# evolver.add_force("Mesh>Harmonic", {"k":{"0":k}, 
#                                     "l0":{"0":l0}})

edge_lengths = compute.edge_lengths()
avg_edge_length= np.mean(edge_lengths)
l0 = avg_edge_length
print("max edge length: " + str(np.max(edge_lengths)) + "  min edge length: " + str(np.min(edge_lengths)) + "  avg edge length" + str(avg_edge_length))
ks = 3.0
m = 2.0
lmax = 5.0*np.max(edge_lengths)
kp = l0**(m+1)*ks/(1-(l0/lmax)**2)
print(kp)
evolver.add_force("Mesh>FENEPOW", {"ks":{"0":str(ks)}, 
                                    "lmax":{"0":str(lmax)},
                                    "m":{"0":str(m)},
                                    "kp":{"0":str(kp)}})

# bending potential
kappa = str(2.0)
# evolver.add_force("Mesh>Bending>Dihedral", {"kappa":{"0":kappa}})
evolver.add_force("Mesh>Bending>Helfrich", {"kappaH":{"1":kappa},
                                            "H0":{"1":str(0)},
                                            "kappaG":{"1":str(0)}})

# check the total bending energy matches analytical results
initial_area = compute.area()
initial_volume = compute.volume()
nu_v = initial_volume/(4*np.pi/3*(initial_area/(4*np.pi))**(3/2))
print("Reduced volume: " + str(nu_v))

#add the Brownian integrator
evolver.add_integrator("Mesh>Brownian>vertex>move", {"seed":"202208"})

dt = str(1e-5)
evolver.set_time_step(dt)

# constant area potential
initial_area = compute.area()
target_area_per_face = 12.0/len(system.faces)
kappa_al = 100/target_area_per_face
print("target area per face:")
print(str(target_area_per_face))
evolver.add_force("Mesh>Constant Area", {
    "sigma": {"1": str(kappa_al)},
    "target_area": {"1": str(target_area_per_face)}
})

# global area potential
kappa_ag = 500
target_global_area = 12.0
print("target global area:")
print(str(target_global_area))
evolver.add_force("Mesh>Constant Global Area", {
    "kappa_ag": str(kappa_ag),
    "target_area": str(target_global_area)
})

# global volume potential
initial_volume = compute.volume()
kappa_v = 1000
target_global_volume = 2.5
print("target global volume:")
print(str(target_global_volume))
evolver.add_force("Mesh>Constant Global Volume", {
    "kappa_v": str(kappa_v),
    "target_volume": str(target_global_volume)
})

## Compute the initial area and volume
print("[Initial] volume:{} area:{}".format(initial_volume, initial_area))

## Compute the initial energy
energy = compute.energy(evolver)
print("[Initial] energy = ", energy)

dump.vtk("BD_snapshot_t0")
for snapshot in range(1, snapshots):
    for temperature in [1e-3, 1e-5]:
        evolver.set_global_temperature(str(temperature))
        evolver.evolveMD(steps=run_steps)
    dump.vtk("BD_snapshot_t" + str(snapshot*run_steps))
    ## Compute the current initial energy
    energy = compute.energy(evolver)
    volume = compute.volume()
    area = compute.area()
    edge_lengths = compute.edge_lengths()
    avg_edge_length= np.mean(edge_lengths)
    print("[Current] energy:{} volume:{} area:{} edge length:{}".format(energy, volume, area, avg_edge_length))
    if snapshot==3:
        dt = str(1e-5)
        evolver.set_time_step(dt)
    #     evolver.delete_force("Mesh>Bending>Helfrich")
    #     kappa = str(3)
    #     evolver.add_force("Mesh>Bending>Helfrich", {"kappaH":{"1":kappa},
    #                                                 "H0":{"1":str(0)},
    #                                                 "kappaG":{"1":str(0)}})
    # if snapshot==10:
    #     dt = str(1e-3)
    #     evolver.set_time_step(dt)

dump.txt("minimized")


