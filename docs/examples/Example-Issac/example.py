#import the code
import pymembrane as mb
import numpy as np
from pprint import pprint

#create a system 
box = mb.Box(100.0, 100.0, 100.0)

system = mb.System(box)

#check if the box is loaded correctly
print(system.box)

#read the mesh
N = 3 #pentagon size
vertex_file = 'vertices_N' + str(N) + '.inp'
face_file = 'faces_N' + str(N) + '.inp'
system.read_mesh_from_files(files={'vertices':vertex_file, 'faces':face_file})

#Print the types:
print('Vertices:')
for v in system.vertices:
    print(v.id, v.type)

print('Faces:')
for f in system.faces:
    print(f.id, f.type)

print('Edges:') 
for e in system.edges:
    print(e.id, e.type)
    
#now if you want to change the type of a vertex, you can do it like this:

vertices = system.vertices

#change the type of the first vertex to 10
for v in vertices:
    v.type = 10

system.vertices = vertices

#Print the types:
print('Vertices:')
for v in system.vertices:
    print(v.id, v.type)
    
#same for the other object