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

os.chdir('/home/ipincus/fork_pymembrane/PyMembrane/docs/examples/08_flat_plane/testing_repulsion/')
os.makedirs('results', exist_ok=True)

## Now we want to have x snapshots every y steps each
# snapshots = 50
# run_steps = 2000
snapshots = 50
run_steps = 1

vertex_file = 'vertices_small.inp'
face_file = 'faces_small.inp'

#create a system 
box = mb.Box(50,50,50)
print(box)
system = mb.System(box)
system.read_mesh_from_files(files={'vertices':vertex_file, 'faces':face_file})
compute = system.compute

#save the mesh to display
#create dumper
dump = system.dumper
dump.vtk("results/initial mesh", False)

#add the evolver class where the potentials and integrators are added
evolver = mb.Evolver(system)

# add particle and force
system.add_particle(0.0,0.0,1.4, 0.4)
# system.add_particle(0.0,2.0,0.0, 0.5)
# system.add_particle(2.0,0.0,0.0, 0.5)
evolver.add_force("Mesh>Particle", {
    "epsilon": {"1": str(1.0)},
    "sigma": {"1": str(1.0)},
    "phi": {"1": str(1.0)}
})

# edge_lengths = compute.edge_lengths()
# avg_edge_length= np.mean(edge_lengths)
# l0 = avg_edge_length
# print("max edge length: " + str(np.max(edge_lengths)) + "  min edge length: " + str(np.min(edge_lengths)))
# ks = 5.0
# m = 2.0
# lmax = 3.0*np.max(edge_lengths)
# kp = l0**(m+1)*ks*lmax**2/(1-(l0/lmax)**2)
# print(kp)
# evolver.add_force("Mesh>FENEPOW", {"ks":{"0":str(ks)}, 
#                                     "lmax":{"0":str(lmax)},
#                                     "m":{"0":str(m)},
#                                     "kp":{"0":str(kp)}})

## Compute the initial energy
energy = compute.energy(evolver)
print("[Initial] energy = ", energy)

# # bending potential
# kappa = str(2.0)
# # evolver.add_force("Mesh>Bending>Dihedral", {"kappa":{"0":kappa}})
# evolver.add_force("Mesh>Bending>Helfrich", {"kappaH":{"1":kappa},
#                                             "H0":{"1":str(0)},
#                                             "kappaG":{"1":str(0)}})

#add the Brownian integrator
evolver.add_integrator("Mesh>Brownian>vertex>move", {"seed":"202208"})

dt = str(1e-4)
evolver.set_time_step(dt)

# constant area potential
initial_area = compute.area()
target_area_per_face = initial_area/len(system.faces)
kappa_al = 100/target_area_per_face
print("target area per face:")
print(str(target_area_per_face))
evolver.add_force("Mesh>Constant Area", {
    "sigma": {"1": str(kappa_al)},
    "target_area": {"1": str(target_area_per_face)}
})

# # global area potential
# kappa_ag = 100
# target_global_area = initial_area
# print("target global area:")
# print(str(target_global_area))
# evolver.add_force("Mesh>Constant Global Area", {
#     "kappa_ag": str(kappa_ag),
#     "target_area": str(target_global_area)
# })

## Compute the initial area and volume
print("[Initial] area:{}".format(initial_area))

## Compute the initial energy
energy = compute.energy(evolver)
print("[Initial] energy = ", energy)

particle_positions = []
dump.vtk("results/BD_snapshot_t0")
for snapshot in range(1, snapshots):
    for temperature in [0.0]:
        evolver.set_global_temperature(str(temperature))
        evolver.evolveMD(steps=run_steps)
    dump.vtk("results/BD_snapshot_t" + str(snapshot*run_steps))
    ## Compute the current initial energy
    energy = compute.energy(evolver)
    area = compute.area()
    print("[Current] energy:{} area:{}".format(energy, area))
    particle_positions.append(system.get_particle_position(0))
    print("Particle center: {}".format(particle_positions[-1]))

dump.txt("results/minimized")


