##############################################################################################################
####  The program visualizes seismic catalog data using RBF Kernel resampling with CPU parallelism.
####  Author: He Pei; 2024.02.25
##############################################################################################################

import sys
import os
import math
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import vtk
from joblib import Parallel, delayed

import CONST  # Import user-defined constants



# Constants
STEP = 0.01  # Grid resolution
Z_OFFSET = 0.01 / CONST.angle_to_kilometers  # Z-axis offset for visualization

# Create output directory for results if it doesn't exist
output_files_path = r'.\resultData\RBF'
os.makedirs(output_files_path, exist_ok=True)



def gaussian_kernel(distance, amplitude, sigma):
    """Compute the Gaussian kernel value based on distance, amplitude, and sigma."""
    return amplitude * np.exp(-0.5 * (distance ** 2) / (sigma ** 2))

def process_point(y_coordinate, x_values, data, data_kdTree):
    """Process each point to calculate RBF value, frequency, and maximum magnitude."""
    print(f"{y_coordinate}  --->  {CONST.latitude_up}")
    
    results = []
    for x in x_values:
        # Convert geodetic coordinates to ECEF
        xyz0 = CONST.geodetic_to_ecef_km(x, y_coordinate, 0.0)
        max_magnitude = 0.0
        rbf_value = 0.0
        
        # Find points within a radius of 4.0
        indices = data_kdTree.query_ball_point([x, y_coordinate], r=4.0)
        
        for index in indices:
            point = data[index]
            xyz1 = CONST.geodetic_to_ecef_km(point[0], point[1], 0.0)
            sigma = math.pow(2.0, point[-1])  # Use magnitude to compute sigma
            distance = np.linalg.norm(np.array(xyz1) - np.array(xyz0))  # Calculate Euclidean distance
            
            # Calculate Gaussian weight if within range
            if distance < 4.0 * sigma:
                amplitude = point[-1] / 1000.0
                gaussian_value = round(gaussian_kernel(distance, amplitude, sigma), 5)
            else:
                gaussian_value = 0.0
            
            rbf_value += gaussian_value
            max_magnitude = max(max_magnitude, point[2])  # Update maximum magnitude if needed
        
        results.append([x, y_coordinate, round(rbf_value, 8), len(indices), max_magnitude])
    
    return results

def main():
    """
    Main function to visualize seismic data.
    
    Usage:
      python src/rbf_Mcpu_vis.py rawData/RBF_csv/CENC_year1980-year1990_minMag1.csv
      
    Inputs:
      - catalog_file: Path to the seismic catalog CSV file

    Outputs:
      - CSV file containing the results of the RBF computation data
      - VTK file containing the visualization data
      - Optional PNG screenshot when 'h' key is pressed during visualization
    """
    
    # Check for command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python src/rbf_Mcpu_vis.py rawData/RBF_csv/CENC_year1980-year1990_minMag1.csv")
        sys.exit()
        
    half_step = round(STEP * 0.5, 5)
    catalog_filename = os.path.basename(sys.argv[1])
    
    # Read data from CSV file
    df = pd.read_csv(sys.argv[1])
    data = np.array(df[['lon', 'lat', 'mag']].values.tolist())
    data_kdTree = cKDTree(data[:, :2])  # Create KD tree for spatial queries

    # Generate grid points
    x_values = [round(x0, 5) for x0 in np.arange(CONST.longitude_left, CONST.longitude_right + half_step, STEP)]
    y_values = [round(y0, 5) for y0 in np.arange(CONST.latitude_down, CONST.latitude_up + half_step, STEP)]

    # Parallel processing of points
    results = Parallel(n_jobs=-1)(delayed(process_point)(yj, x_values, data, data_kdTree) for yj in y_values)

    # Write results to CSV file
    csv_filename = os.path.join(output_files_path, f"RBF_{catalog_filename}")
    with open(csv_filename, "w") as file:
        file.write("lon,lat,rbf_value,freq,max_mag\n")
        for result in results:
            for r in result:
                file.write(','.join(map(str, r)) + '\n')

    # Visualization
    colors = vtk.vtkNamedColors()
    ugrid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()

    # Create VTK arrays for RBF value, frequency, and max magnitude
    rbf_array = vtk.vtkFloatArray()
    rbf_array.SetNumberOfComponents(1)
    rbf_array.SetName('RBF value')
    
    freq_array = vtk.vtkFloatArray()
    freq_array.SetNumberOfComponents(1)
    freq_array.SetName('frequency')

    max_mag_array = vtk.vtkFloatArray()
    max_mag_array.SetNumberOfComponents(1)
    max_mag_array.SetName('MAX magnitude')
    
    # Insert points and associated data into VTK structures
    for lin in results:
        for line in lin:
            points.InsertNextPoint(line[0], line[1], line[2] + Z_OFFSET)
            rbf_array.InsertNextValue(line[2])
            freq_array.InsertNextValue(line[3])
            max_mag_array.InsertNextValue(line[4])
    
    # Create the unstructured grid cells
    x_num = len(x_values)
    y_num = len(y_values)
    for j in range(y_num - 1):
        for i in range(x_num - 1):
            id0 = i + j * x_num
            id1 = id0 + 1
            id2 = id1 + x_num
            id3 = id0 + x_num
            
            # Create triangles for the grid
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id3, id0, id1])
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id3, id1, id2])

    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(freq_array)
    ugrid.GetPointData().AddArray(max_mag_array)
    ugrid.GetPointData().AddArray(rbf_array)

    # Write VTK data to file
    vtk_filename = os.path.join(output_files_path, f"RBF_{catalog_filename.replace('csv', 'vtu')}")
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(vtk_filename)
    writer.SetDataModeToBinary()
    writer.Write()

    # Create a color transfer function for RBF values
    color_transfer_function = vtk.vtkColorTransferFunction()
    color_transfer_function.SetColorSpaceToHSV()  # Set to HSV for better control

    # Define color mapping points
    color_points = [
        (0.01,     0.95, 0.95, 0.95),
        (0.02,     0.0, 0.0, 0.5625),
        (0.342222, 0.0, 0.0, 1.0),
        (1.10159,  0.0, 1.0, 1.0),
        (1.48127,  0.5, 1.0, 0.5),
        (1.86095,  1.0, 1.0, 0.0),
        (2.62032,  1.0, 0.0, 0.0),
        (3,        0.5, 0.0, 0.0)
    ]
    
    for point in color_points:
        color_transfer_function.AddRGBPoint(*point)

    # Create a mapper and actor for the visualization
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)
    mapper.SetLookupTable(color_transfer_function)
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("RBF value")
    mapper.SetColorModeToMapScalars()
    mapper.SetScalarRange(1, 3)  # Set scalar range for coloring

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('Silver'))
    actor.GetProperty().SetPointSize(2)

    # Create a scalar bar for visualization
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_transfer_function)
    scalar_bar.SetTitle("RBF")
    scalar_bar.GetLabelTextProperty().SetFontSize(20)
    scalar_bar.SetOrientationToVertical()
    scalar_bar.SetPosition(0.95, 0.1)
    scalar_bar.SetWidth(0.03)
    scalar_bar.SetHeight(0.8)

    # Set properties for annotation labels
    label_text_property = scalar_bar.GetAnnotationTextProperty()
    label_text_property.SetFontSize(20)
    label_text_property.SetColor(0, 0, 0)  # Black color
    label_text_property.SetJustificationToCentered()
    label_text_property.SetVerticalJustificationToCentered()

    # Set properties for the title
    title_text_property = scalar_bar.GetTitleTextProperty()
    title_text_property.SetFontSize(50)
    title_text_property.SetColor(0, 0, 0)  # Black color

    # Setup the visualization pipeline
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor(scalar_bar)
    renderer.SetBackground(colors.GetColor3d('Lavender'))

    # Create render window and interactor
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(1920, 1080)
    render_window.AddRenderer(renderer)

    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Add keyboard event for saving screenshot
    def save_screenshot_callback(obj, event):
        render_window.SetWindowName("Seismic Data Visualization")
        screenshot_filter = vtk.vtkWindowToImageFilter()
        screenshot_filter.SetInput(render_window)
        screenshot_filter.Update()
        
        png_writer = vtk.vtkPNGWriter()
        png_writer.SetFileName(os.path.join(output_files_path, "screenshot.png"))
        png_writer.SetInputData(screenshot_filter.GetOutput())
        png_writer.Write()
    
    render_window_interactor.AddObserver("KeyPressEvent", save_screenshot_callback)

    # Start rendering
    render_window.Render()
    render_window_interactor.Start()

if __name__ == '__main__':
    main()
