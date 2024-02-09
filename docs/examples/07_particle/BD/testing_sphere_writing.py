import os
import numpy as np
import vtk

p_radius = 1.0
particle_positions = np.array([[1,0,0],[2,0,0]])

# Function to create a sphere source with varying position over time
def create_time_varying_sphere(radius, centers):
    spheres = []

    for current_center in centers:

        # Create a sphere source for the current time step
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(radius)
        sphere.SetCenter(current_center)
        sphere.Update()

        spheres.append(sphere.GetOutput())

    return spheres

# Function to write a sequence of spheres to a VTK file
def write_time_varying_spheres_to_vtk(spheres, filename):
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(filename)

    # Create a multi-block dataset to store spheres at different time steps
    multi_block = vtk.vtkMultiBlockDataSet()

    for i, sphere in enumerate(spheres):
        # Create a block for each time step
        block = vtk.vtkPolyData()
        block.ShallowCopy(sphere)
        multi_block.SetBlock(i, block)

    # Write the multi-block dataset to the VTK file
    writer.SetInputData(multi_block)
    writer.Write()

# Example usage: Create time-varying spheres, write them to a VTK file
spheres_output = create_time_varying_sphere(p_radius, particle_positions)
write_time_varying_spheres_to_vtk(spheres_output, "results/time_varying_spheres.vtm")
