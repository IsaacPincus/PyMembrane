import open3d as o3d
import numpy as np

mesh = o3d.io.read_triangle_mesh('plane_10x10.ply')

vertices = np.asarray(mesh.vertices)
triangles = np.asarray(mesh.triangles)

with open('vertices.inp', 'w') as file:
    for index, vertex in enumerate(vertices):
        file.write(str(index) + ' ' + ' '.join(map(str, vertex)) + ' 1 \n')


with open('faces.inp', 'w') as file:
    for index, face in enumerate(triangles):
        file.write(str(index) + ' ' + ' '.join(map(str, face)) + ' 1 1 \n')