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

os.chdir('/home/ipincus/fork_pymembrane/PyMembrane/docs/examples/06_tet_tests/BD/')

# set up the vertices
vertex_file = 'vertices.inp'
face_file = 'faces.inp'

#create a system 
box = mb.Box(40,40,40)
system = mb.System(box)
system.read_mesh_from_files(files={'vertices':vertex_file, 'faces':face_file})

vertices = system.vertices
for vertex in vertices:
    print(str(vertex.id) + ' ' + str(vertex.r))

faces = system.faces
for face in faces:
    print(str(face.id) + ' ' + str(face.v1) + ' ' + str(face.v2) + ' ' + str(face.v3))

#save the mesh to display
#create dumper
dump = system.dumper
dump.vtk("initial mesh", False)

#add the evolver class where the potentials and integrators are added
evolver = mb.Evolver(system)
evolver.add_integrator("Mesh>Brownian>vertex>move", {"seed":"5"})
evolver.set_time_step(str(0.2))

# compute module
compute = system.compute

# global volume potential
initial_volume = compute.volume()
kappa_v = 1
target_global_volume = 0.5
print("Init volume:{} Target volume:{} ".format(initial_volume, target_global_volume))
evolver.add_force("Mesh>Constant Global Volume", {
    "kappa_v": str(kappa_v),
    "target_volume": str(target_global_volume)
})

# # global area potential
# initial_area = compute.area()
# kappa_ag = 1
# target_global_area = 4
# print("Init area:{} Target area:{} ".format(initial_area, target_global_area))
# evolver.add_force("Mesh>Constant Global Area", {
#     "kappa_ag": str(kappa_ag),
#     "target_area": str(target_global_area)
# })

# # constant area potential
# initial_area_per_face = compute.area()/len(system.faces)
# target_area_per_face = 1
# kappa_al = 1
# print("Init area per face:{} Target area per face:{} ".format(initial_area_per_face, target_area_per_face))
# evolver.add_force("Mesh>Constant Area", {
#     "sigma": {"0": str(kappa_al)},
#     "target_area": {"0": str(target_area_per_face)}
# })

forces = compute.forces(evolver)
for force in forces:
    print(forces)

dump.vtk("mesh_t0")

run_steps = 1
evolver.set_global_temperature(str(0.0))
for t in range(1,10):
    evolver.evolveMD(steps=run_steps)
    dump.vtk("mesh_t" + str(t*run_steps))
    # area = compute.area()
    # print("Init area:{} Target area:{} Current area:{} ".format(initial_area, target_global_area, area))
    volume = compute.volume()
    print("Init volume:{} Target volume:{} Current volume:{} ".format(initial_volume, target_global_volume, volume))

# volume = compute.volume()
# print("Init volume:{} Target volume:{} Current volume:{} ".format(initial_volume, target_global_volume, volume))

# area = compute.area()
# print("Init area:{} Target area:{} Current area:{} ".format(initial_area, target_global_area, area))

# # first we need to know the edge length to move it appropriate:
# 
# edge_lengths = compute.edge_lengths()
# avg_edge_length= np.mean(edge_lengths)
# print("[Initial] avg_edge_length = ", avg_edge_length)
# k = str(400.0)
# l0 = str(avg_edge_length)
# # evolver.add_force("Mesh>Harmonic", {"k":{"0":k}, 
# #                                     "l0":{"0":l0}})

# # limit potential
# lmin = str(0.5*np.min(edge_lengths))
# lmax = str(2.0*np.max(edge_lengths))
# evolver.add_force("Mesh>Limit", {"lmin":{"0":lmin}, 
#                                  "lmax":{"0":lmax}})

# # bending potential
# kappa = str(2.0)
# # evolver.add_force("Mesh>Bending>Dihedral", {"kappa":{"0":kappa}})
# evolver.add_force("Mesh>Bending>Helfrich", {"kappaH":{"0":kappa},
#                                             "H0":{"0":str(0)},
#                                             "kappaG":{"0":str(0)}})

# #add the Brownian integrator
# evolver.add_integrator("Mesh>Brownian>vertex>move", {"seed":"202208"})

# dt = str(1e-4)
# 

# ## then we want to run the simulation for a temperature 1e-4
# temperature = str(1e-4)
# evolver.set_global_temperature(temperature)

# # constant area potential
# initial_area_per_face = compute.area()
# target_area_per_face = initial_area_per_face/len(system.faces)
# kappa_al = 1e2/target_area_per_face
# print("target area per face:")
# print(str(target_area_per_face))
# evolver.add_force("Mesh>Constant Area", {
#     "sigma": {"1": str(kappa_al)},
#     "target_area": {"1": str(target_area_per_face)}
# })

# # global area potential
# kappa_ag = 1e2
# target_global_area = initial_area_per_face
# print("target global area:")
# print(str(target_global_area))
# evolver.add_force("Mesh>Constant Global Area", {
#     "kappa_ag": str(kappa_ag),
#     "target_area": str(target_global_area)
# })

# # global volume potential
# initial_volume = compute.volume()
# kappa_v = 5e2
# target_global_volume = initial_volume*0.9
# print("target global volume:")
# print(str(target_global_volume))
# evolver.add_force("Mesh>Constant Global Volume", {
#     "kappa_v": str(kappa_v),
#     "target_volume": str(target_global_volume)
# })

# # Randomize the mesh to avoid initial freezing
# rnd_move = 0.5*np.std(edge_lengths)/6.0
# vertices = system.vertices
# for vertex in vertices:
#     vertex.r.x += np.random.uniform(-rnd_move, rnd_move)
#     vertex.r.y += np.random.uniform(-rnd_move, rnd_move)
#     vertex.r.z += np.random.uniform(-rnd_move, rnd_move)
# system.vertices = vertices

# ## Compute the initial energy
# min_energy = snapshots*[None]
# min_energy[0] = compute.energy(evolver)['edges']/system.Numedges
# print("[Initial] energy:{} volume:{} area:{}".format(min_energy[0], initial_volume, initial_area_per_face))

# ## Compute the initial energy
# energy = compute.energy(evolver)
# print("[Initial] energy = ", energy)

# dump.vtk("BD_snapshot_t0")
# for snapshot in range(1, snapshots):
#     for temperature in [1e-6]:
#         evolver.set_global_temperature(str(temperature))
#         evolver.evolveMD(steps=run_steps)
#     dump.vtk("BD_snapshot_t" + str(snapshot*run_steps))
#     ## Compute the current initial energy
#     energy = compute.energy(evolver)
#     volume = compute.volume()
#     area = compute.area()
#     print("[Current] energy:{} volume:{} area:{}".format(energy, volume, area))

# # Compute the final volume
# energy = compute.energy(evolver)
# print("[Final] energy:{} volume:{}".format(min_energy[snapshots-1], compute.volume()))
# final_volume = compute.volume()
# print("volume difference: {}".format(final_volume-initial_volume))
