#######################################################################
####  The program is used to plot M-T diagrams of seismic catalog data.
####  Author: He Pei; 2024.02.25                                       
#######################################################################

import vtk
import pandas as pd
from datetime import datetime
import sys
import os

# Add the current directory to the system path for importing constants
sys.path.append(r".")
import CONST  # User-defined constants

# Define constants
TIME_FORMAT = "%Y.%m.%d %H:%M:%S.%f"
STEP_MONTHS = 12

# Create output directory for results if it doesn't exist
output_files_path = r'.\resultData\MT'
os.makedirs(output_files_path, exist_ok=True)  # Create directory if it doesn't exist

def main():
    """
    Main function to visualize seismic data.
    
    Usage:
      python src/MT_plot.py rawData/CENC_catalog_1970-2020.csv

    Inputs:
      - catalog_file: Path to the seismic catalog CSV file

    Outputs:
      - VTK file containing the visualization data
      - Optional PNG screenshot when 'h' key is pressed during visualization
    """

    # Check for command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python src/MT_plot.py rawData/CENC_catalog_1970-2020.csv")
        sys.exit()

    catalog_filename = os.path.basename(sys.argv[1])  # Extract filename from argument

    # Create polydata for geometric data
    lines_polydata = vtk.vtkPolyData()
    named_colors = vtk.vtkNamedColors()

    points = vtk.vtkPoints()
    cell_lines = vtk.vtkCellArray()
    magnitude = vtk.vtkFloatArray()
    magnitude.SetNumberOfComponents(1)
    magnitude.SetName('magnitude')

    def seconds_between_times(time1, time2):
        """Calculate the number of months between two time strings."""
        time1_dt = datetime.strptime(time1, TIME_FORMAT)
        time2_dt = datetime.strptime(time2, TIME_FORMAT)
        return abs((time2_dt - time1_dt).total_seconds()) / (3600.0 * 24 * 30 * STEP_MONTHS)

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

    # Calculate scaling factor for the plot
    k = (year1 + 1 - year0) / (x1 - x0)
    print(f"{1 / k},  {-1.0 / k * year0}\n")

    # Populate points and lines for visualization
    for i in range(len(mag)):
        poly_line = vtk.vtkLine()
        
        x = seconds_between_times(t0, time[i])
        points.InsertNextPoint(x, y0, 0.0)  # Start point
        points.InsertNextPoint(x, mag[i], 0.0)  # End point
        
        magnitude.InsertNextValue(mag[i])
        
        poly_line.GetPointIds().SetId(0, 2 * i)
        poly_line.GetPointIds().SetId(1, 2 * i + 1)
        
        cell_lines.InsertNextCell(poly_line)

    lines_polydata.SetPoints(points)
    lines_polydata.SetLines(cell_lines)
    lines_polydata.GetCellData().AddArray(magnitude)

    # Setup the writer for saving the VTK data
    output_VTK_filename = os.path.join(output_files_path, f"MT_{catalog_filename.replace('csv', 'vtu')}")
    writer1 = vtk.vtkXMLDataSetWriter()
    writer1.SetInputData(lines_polydata)
    writer1.SetFileName(output_VTK_filename)
    writer1.SetDataModeToBinary()
    writer1.Update()

    # Setup the visualization pipeline
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(lines_polydata)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetLineWidth(1.5)
    actor.GetProperty().SetColor(1.0, 0.0, 0.0)  # Set line color to red

    # Renderer and window setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(named_colors.GetColor3d("White"))
    renderer.ResetCamera()

    render_window = vtk.vtkRenderWindow()
    render_window.SetWindowName('M-T diagrams of the Earthquake Catalog')
    render_window.SetSize(100, 1_000)  # Set render window size
    render_window.AddRenderer(renderer)

    # Create render window interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Define key press event callback function
    def save_screenshot(obj, event):
        key = obj.GetKeySym()
        if key == 'h':
            # Create a vtkWindowToImageFilter object
            window_to_image_filter = vtk.vtkWindowToImageFilter()
            window_to_image_filter.SetInput(render_window)
            window_to_image_filter.SetScale(2)  # Increase resolution
            window_to_image_filter.Update()
            
            output_PNG_filename = os.path.join(
                output_files_path, f"MT_{catalog_filename.replace('csv', 'png')}")
            # Create a vtkPNGWriter object
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
