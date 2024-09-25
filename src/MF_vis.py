##########################################################################
####  The program is used to plot hotmap diagrams of seismic catalog data.
####  Author: He Pei; 2024.02.25                                          
##########################################################################

import sys
import os
import vtk
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm  # For progress bar

import CONST  # Import user-defined constants



# Configuration parameters
step = 0.02  # Grid step size
window_Width = 2  # Width of the search window
z0 = 5000.0 / CONST.angle_to_meters  # Maximum elevation adjusted for translation in Z direction
scale = 0.001  # Scale factor for frequency representation

# Output directory for results
output_files_path = r'.\resultData\MF'
if not os.path.exists(output_files_path):
    os.makedirs(output_files_path)



def main():
    """
    Main function to visualize seismic data.
    
    Usage:
      python src/MF_vis.py rawData/CENC_catalog_1970-2020.csv
      
    Inputs:
      - catalog_file: Path to the seismic catalog CSV file

    Outputs:
      - VTK file containing the visualization data
      - Optional PNG screenshot when 'h' key is pressed during visualization
    """
    
    # Check command-line arguments
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python src/MF_vis.py rawData/CENC_catalog_1970-2020.csv")
        sys.exit()
        
    catalog_filename = os.path.basename(sys.argv[1])

    half_step = round(step * 0.5, 5)
    search_radius = window_Width * step / 2.0
    core_radius = search_radius * 1.5

    # Initialize VTK structures
    points = vtk.vtkPoints()
    ugrid = vtk.vtkUnstructuredGrid()
    freq = vtk.vtkFloatArray()
    freq.SetNumberOfComponents(1)
    freq.SetName('frequency')

    max_mag = vtk.vtkFloatArray()
    max_mag.SetNumberOfComponents(1)
    max_mag.SetName('MAX Magnitude')
    
    # Read seismic data from CSV
    df = pd.read_csv(sys.argv[1])
    data = np.array(df[['lon', 'lat', 'mag']].values.tolist())
    data_kdTree = cKDTree(data[:, :2])  # KD-tree for spatial querying
    max_magnitude = df['mag'].max()

    # Define grid coordinates
    x = [round(x0, 5) for x0 in np.arange(CONST.longitude_left, CONST.longitude_right + half_step, step)]
    y = [round(y0, 5) for y0 in np.arange(CONST.latitude_down, CONST.latitude_up + half_step, step)]

    # Progress bar for processing
    for yj in tqdm(y, desc="Processing Rows"):
        for xi in x:
            freq0 = 0
            mag0 = 0.0
            
            # Query for nearby points
            indices = data_kdTree.query_ball_point([xi, yj], r=core_radius)
            if indices is not None:
                MF_list = []
                
                for index in indices:
                    p = data[index]
                    if (p[0] >= xi - search_radius) and (p[0] < xi + search_radius) and (p[1] >= yj - search_radius) and (p[1] < yj + search_radius):
                        MF_list.append(index)
                        if p[2] > mag0:
                            mag0 = p[2]
                
                freq0 = len(MF_list)        
            
            # Insert point data into VTK structures
            points.InsertNextPoint(xi, yj, scale * freq0 + z0)
            freq.InsertNextValue(freq0)
            max_mag.InsertNextValue(mag0)

    # Create unstructured grid from points
    x_num = len(x)
    y_num = len(y)

    for j in range(y_num - 1):
        for i in range(x_num - 1):
            id0 = i + j * x_num
            id1 = id0 + 1
            id2 = id1 + x_num
            id3 = id0 + x_num
            
            # Define triangles for grid
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id3, id0, id1])
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id3, id1, id2])

    # Add data arrays to the grid
    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(freq)
    ugrid.GetPointData().AddArray(max_mag)

    # Write data to '.VTK' file
    writer1 = vtk.vtkXMLUnstructuredGridWriter()
    output_VTK_filename = os.path.join(  output_files_path, f"M-F_{catalog_filename.replace('csv', 'vtu')}"  )
    writer1.SetFileName(output_VTK_filename)
    writer1.SetInputData(ugrid)
    writer1.SetDataModeToBinary()
    writer1.Write()
    
    # Define magnitude ranges and corresponding colors
    color_map = [
        [1.0, 3.0, (0.9, 0.9, 0.9)],                                           # 1.0 <= mag < 3.0
        [3.0, 5.0, (1.0, 0.0, 0.498039)],                                           # 3.0 <= mag < 5.0
        [5.0, max_magnitude if max_magnitude > 8.0 else 8.0, (1.0, 0.0, 0.0)]  # 5.0 <= mag <= max_magnitude
    ]
    
    # Create a Discrete Color Transfer Function
    color_transfer_func = vtk.vtkDiscretizableColorTransferFunction()
    color_transfer_func.DiscretizeOn()  # Enable discretization
    color_transfer_func.SetNumberOfValues(len(color_map))
    
    # Add color segments and annotations to the color transfer function
    mag_mid = 0.0
    d_label = (max_magnitude - 1.0) / len(color_map) if max_magnitude > 8.0 else (8.0 - 1.0) / len(color_map)
    for mag_min, mag_max, color in color_map:
        color_transfer_func.AddRGBSegment(mag_min, *color, mag_max, *color)  # Add color segment
        mag_mid += d_label
        label = f"{mag_min}~{mag_max}"
        color_transfer_func.SetAnnotation(mag_mid, label)

    # Create mapper and actor for visualization
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)
    mapper.SetLookupTable(color_transfer_func)
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("MAX Magnitude")
    mapper.SetColorModeToMapScalars()
    mapper.SetScalarRange(1.0, max_magnitude if max_magnitude > 8.0 else 8.0)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create color bar (Scalar Bar)
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_transfer_func)
    scalar_bar.SetTitle("MAX_Magnitude")
    scalar_bar.SetTextPositionToSucceedScalarBar()  # Labels on the right
    scalar_bar.DrawTickLabelsOff()  # Disable tick labels
    scalar_bar.DrawAnnotationsOn()  # Enable annotations

    # Adjust annotation properties
    label_text_property = scalar_bar.GetAnnotationTextProperty()
    label_text_property.SetFontSize(20)
    label_text_property.SetColor(0, 0, 0)  # Black color
    label_text_property.SetJustificationToCentered()
    label_text_property.SetVerticalJustificationToCentered()
    
    # Adjust title properties
    title_text_property = scalar_bar.GetTitleTextProperty()
    title_text_property.SetFontSize(5)
    title_text_property.SetColor(0, 0, 0)  # Black color

    # Adjust Scalar Bar size and position
    scalar_bar.SetPosition(0.95, 0.1)
    scalar_bar.SetWidth(0.03)
    scalar_bar.SetHeight(0.8)

    # Create renderer and render window
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor2D(scalar_bar)
    renderer.SetBackground(1.0, 1.0, 1.0)  # White background
    renderer.ResetCamera()

    render_window = vtk.vtkRenderWindow()
    render_window.SetWindowName(f'M-F diagrams of Earthquake Catalog Scatter Visualization')
    render_window.SetSize(1600, 900)  # Set render window size
    render_window.AddRenderer(renderer)

    # Create render window interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Define key press event callback function
    def save_screenshot(obj, event):
        key = obj.GetKeySym()
        if key == 'h':
            window_to_image_filter = vtk.vtkWindowToImageFilter()
            window_to_image_filter.SetInput(render_window)
            window_to_image_filter.SetScale(2)  # Increase resolution
            window_to_image_filter.ReadFrontBufferOff()
            window_to_image_filter.Update()

            output_PNG_filename = os.path.join(output_files_path, f"{catalog_filename.replace('csv', 'png')}")
            png_writer = vtk.vtkPNGWriter()
            png_writer.SetFileName(output_PNG_filename)
            png_writer.SetInputConnection(window_to_image_filter.GetOutputPort())
            png_writer.Write()
            print(f"Screenshot saved as {output_PNG_filename}.")

    # Add key press event listener
    render_window_interactor.AddObserver("KeyPressEvent", save_screenshot)

    # Start the visualization
    render_window.Render()
    render_window_interactor.Start()

if __name__ == "__main__":
    main()
