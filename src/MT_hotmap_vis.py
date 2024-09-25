################################################################################
####  The program is used to plot M-T-F hotmap diagrams of seismic catalog data.
####  Author: He Pei; 2024.02.25                                                 
################################################################################

import vtk
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from datetime import datetime
import os
import sys

# Add the current directory to the system path for importing constants
sys.path.append(r".")
import CONST  # User-defined constants



# Define constants
TIME_FORMAT = "%Y.%m.%d %H:%M:%S.%f"
STEP_MONTHS = 12
STEP_MAG = 0.1

# Create output directory for results if it doesn't exist
output_files_path = r'.\resultData\MT_hotmap'
os.makedirs(output_files_path, exist_ok=True)  # Create directory if it doesn't exist



def seconds_between_times(time1, time2):
    """Calculate the number of months between two time strings."""
    time1_dt = datetime.strptime(time1, TIME_FORMAT)
    time2_dt = datetime.strptime(time2, TIME_FORMAT)
    return abs((time2_dt - time1_dt).total_seconds()) / (3600.0 * 24 * 30 * STEP_MONTHS)

def main():
    """
    Main function to visualize seismic data.
    
    Usage:
      python src/MT_hotmap_vis.py rawData/CENC_minMag1_1970-2023.csv

    Inputs:
      - catalog_file: Path to the seismic catalog CSV file

    Outputs:
      - VTK file containing the visualization data
      - Optional PNG screenshot when 'h' key is pressed during visualization
    """

    # Check for command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python src/MT_hotmap_vis.py rawData/CENC_minMag1_1970-2023.csv")
        sys.exit()

    catalog_filename = os.path.basename(sys.argv[1])

    # Initialize VTK structures
    colors = vtk.vtkNamedColors()
    ugrid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()

    # Initialize data arrays for maximum magnitude and frequency
    max_mag = vtk.vtkFloatArray()
    max_mag.SetName('MAX_Magnitude')

    freq = vtk.vtkFloatArray()
    freq.SetName('frequency')

    # Read seismic catalog data
    df = pd.read_csv(sys.argv[1])
    time = df['time'].to_list()
    mag = df['mag']

    # Determine time range for the plot
    year0 = datetime.strptime(time[0], TIME_FORMAT).year
    year1 = datetime.strptime(time[-1], TIME_FORMAT).year

    t0 = f"{year0}.01.01 00:00:00.00"
    x0 = 0.0
    x1 = seconds_between_times(t0, f"{year1 + 1}.01.01 00:00:00.00")
    y0 = 0.0
    y1 = max(mag) + 0.5
    
    k = (year1 + 1 - year0) / (x1 - x0)
    print(f"{1 / k},  {-1.0 / k * year0} \n")

    # Calculate points
    xy = [[seconds_between_times(t0, time[i]), mag[i]] for i in range(len(mag))]
    xy_kdTree = cKDTree(xy)

    # Create grid points
    half_mag_step = STEP_MAG * 0.5
    x_step = STEP_MAG
    half_x_step = x_step * 0.5
    kdTree_radius = half_x_step + half_mag_step

    x = np.round(np.arange(x0, x1 + half_x_step, x_step), 5)
    y = np.round(np.arange(y0, y1 + half_mag_step, STEP_MAG), 5)
    max_frequency = -float('inf')

    # Populate points and calculate frequency
    for yj in y:
        for xi in x:
            p_list = []
            mag0 = 0.0
            indices = xy_kdTree.query_ball_point([xi, yj], r=kdTree_radius)

            for index in indices:
                p = xy[index]
                if (p[0] >= xi - half_x_step) and (p[0] < xi + half_x_step) and (p[1] >= yj - half_mag_step) and (p[1] < yj + half_mag_step):
                    p_list.append(index)
                    mag0 = max(mag0, p[1])

            freq0 = len(p_list) + 1
            
            points.InsertNextPoint(xi, yj, 0.0)
            freq.InsertNextValue(freq0)
            max_mag.InsertNextValue(mag0)
            
            max_frequency = max(max_frequency, freq0)

    # Construct the unstructured grid
    x_num = len(x)
    y_num = len(y)
    
    for j in range(y_num - 1):
        for i in range(x_num - 1):
            id0 = i + j * x_num
            id1 = id0 + 1
            id2 = id1 + x_num
            id3 = id0 + x_num
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id3, id0, id1])
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id3, id1, id2])

    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(max_mag)
    ugrid.GetPointData().AddArray(freq)

    # Setup the writer for saving the VTK data
    output_VTK_filename = os.path.join(output_files_path, f"MT_hotmap_{catalog_filename.replace('csv', 'vtu')}")
    writer1 = vtk.vtkXMLUnstructuredGridWriter()
    writer1.SetInputData(ugrid)
    writer1.SetFileName(output_VTK_filename)
    writer1.SetDataModeToBinary()
    writer1.Write()  # Use Write instead of Update for writer to save file

    # Create a color transfer function for jet color mapping
    color_transfer_function = vtk.vtkColorTransferFunction()
    color_transfer_function.SetColorSpaceToHSV()  # Optional: set to HSV color space for better control

    # Add color points for 'jet' color mapping
    color_points = [
        (1, 1, 1, 1),                 # White
        (2, 0, 0, 0.5625),            # Dark Blue
        (77.7777, 0, 0, 1),           # Blue
        (253.27, 0, 1, 1),            # Cyan
        (314.016, 0.5, 1, 0.5),       # Green
        (428.762, 1, 1, 0),           # Yellow
        (604.254, 1, 0, 0),           # Red
        (max_frequency, 0.5, 0, 0)    # Dark Red
    ]
    
    for point in color_points:
        color_transfer_function.AddRGBPoint(*point)

    # Create a mapper and actor
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)
    mapper.SetLookupTable(color_transfer_function)
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("frequency")
    mapper.SetColorModeToMapScalars()
    mapper.SetScalarRange(1, max_frequency)  # Use the maximum magnitude for coloring

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('Silver'))
    actor.GetProperty().SetPointSize(2)

    # Create and add color bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_transfer_function)
    scalar_bar.SetTitle("frequency")
    scalar_bar.GetLabelTextProperty().SetFontSize(20)
    scalar_bar.SetOrientationToVertical()
    scalar_bar.SetPosition(0.95, 0.1)
    scalar_bar.SetWidth(0.03)
    scalar_bar.SetHeight(0.8)

    # Adjust annotation properties
    label_text_property = scalar_bar.GetAnnotationTextProperty()
    label_text_property.SetFontSize(20)
    label_text_property.SetColor(0, 0, 0)  # Black color
    label_text_property.SetJustificationToCentered()
    label_text_property.SetVerticalJustificationToCentered()

    # Adjust title properties
    title_text_property = scalar_bar.GetTitleTextProperty()
    title_text_property.SetFontSize(50)
    title_text_property.SetColor(0, 0, 0)  # Black color

    # Setup visualization pipeline
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor2D(scalar_bar)
    renderer.SetBackground(colors.GetColor3d('White'))
    renderer.ResetCamera()

    render_window = vtk.vtkRenderWindow()
    render_window.SetWindowName('M-T Hotmap diagrams of the Earthquake Catalog')
    render_window.SetSize(200, 800)  # Set render window size
    render_window.AddRenderer(renderer)

    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Define key press event callback function
    def save_screenshot(obj, event):
        key = obj.GetKeySym()
        if key == 'h':
            window_to_image_filter = vtk.vtkWindowToImageFilter()
            window_to_image_filter.SetInput(render_window)
            window_to_image_filter.SetScale(2)  # Increase resolution
            window_to_image_filter.Update()
            
            output_PNG_filename = os.path.join(  output_files_path, f"MT_hotmap_{catalog_filename.replace('csv', 'png')}"  )
            writer2 = vtk.vtkPNGWriter()
            writer2.SetFileName(output_PNG_filename)
            writer2.SetInputConnection(window_to_image_filter.GetOutputPort())
            writer2.Write()

            print(f"Screenshot saved as {output_PNG_filename}.")

    # Bind the key press event for saving the image
    render_window_interactor.AddObserver("KeyPressEvent", save_screenshot)

    # Start rendering and interaction
    render_window.Render()
    print("Press 'h' to save a screenshot.")
    render_window_interactor.Start()

if __name__ == "__main__":
    main()
