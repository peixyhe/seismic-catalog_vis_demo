###############################################################################
# The program visualizes seismic catalog data as spheres in 2D or 3D space.
# Author: He Pei, at 110 office room
# Date: 2024.02.25
###############################################################################

import sys
import os
import csv
import vtk
import pandas as pd

# Import user-defined constants
import CONST  # Ensure CONST.py is in the same directory or in Python's PATH



# Maximum elevation adjusted for translation in the Z direction
z0 = 5000.0 / CONST.angle_to_meters  # max(elevation) = 7169.34 meters
sphere_scale = 0.35

output_files_path = r'.\resultData\catalog_scatter'
if not os.path.exists(output_files_path):
    os.makedirs(output_files_path)



def main():
    """
    Main function to visualize seismic data.
    Usage:
      python src/catalog_scatter_vis.py 2D rawData/CENC_catalog_1970-2020.csv
      or
      python src/catalog_scatter_vis.py 3D rawData/CENC_catalog_1970-2020.csv

    Inputs:
      - mode: '2D' or '3D' visualization
      - catalog_file: Path to the seismic catalog CSV file

    Outputs:
      - VTK file containing the visualization data
      - VTK file containing the visualization data
      - Optional PNG screenshot when 'h' key is pressed during visualization
    """

    # Check command-line arguments
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python src/catalog_scatter_vis.py 2D rawData/CENC_catalog_1970-2020.csv")
        print("  or")
        print("  python src/catalog_scatter_vis.py 3D rawData/CENC_catalog_1970-2020.csv")
        sys.exit()

    mode = sys.argv[1]
    catalog_filename = os.path.basename(sys.argv[2])

    # Read catalog dataset
    df = pd.read_csv(sys.argv[2])
    try:
        data = df[['lon', 'lat', 'dep', 'mag']].values
    except KeyError:
        try:
            data = df[['longitude', 'latitude', 'depth', 'magnitude']].values
        except KeyError:
            print("The input of earthquake catalog dataset ERROR:")
            print(f"  {catalog_filename} must be a standard CSV file with headers:")
            print("  'lon','lat','dep','mag' or 'longitude','latitude','depth','magnitude'")
            sys.exit()

    # Filter catalog dataset based on minimum magnitude
    data = data[data[:, 3] >= CONST.MIN_MAGNITUDE]
    if data.size == 0:
        print("No data points meet the minimum magnitude requirement.")
        sys.exit()

    # Get maximum magnitude for color mapping
    max_magnitude = data[:, 3].max()
    max_radius = CONST.Get_mag2radius(max_magnitude)

    # Create VTK point data and arrays
    vtk_points = vtk.vtkPoints()
    
    diameter_array = vtk.vtkFloatArray()
    diameter_array.SetName("diameter")
    
    depth_array = vtk.vtkFloatArray()
    depth_array.SetName("depth")
    
    magnitude_array = vtk.vtkFloatArray()
    magnitude_array.SetName("magnitude")

    # Define output CSV filename
    output_CSV_filename = os.path.join(
        output_files_path, f"{mode}_catalog_minMAG{CONST.MIN_MAGNITUDE}_scatter_{catalog_filename}")

    # Populate VTK point data and write to CSV
    with open(output_CSV_filename, 'w', newline='', encoding='utf-8') as fw:
        writer0 = csv.writer(fw)
        # Write header
        writer0.writerow(['longitude', 'latitude', 'depth', 'magnitude', 'radius', 'z0_2d', 'z_3d'])
        for lon, lat, depth, magnitude in data:
            z = -depth / CONST.angle_to_kilometers
            if mode == "2D":
                z_coord = z0
            elif mode == "3D":
                z_coord = z + z0
            else:
                print("Error: The first parameter must be '2D' or '3D'.")
                sys.exit()

            vtk_points.InsertNextPoint(lon, lat, z_coord)
            
            r0 = CONST.Get_mag2radius(magnitude)
            diameter = (r0 / max_radius) * sphere_scale
            
            diameter_array.InsertNextValue(diameter)
            depth_array.InsertNextValue(depth)
            magnitude_array.InsertNextValue(magnitude)

            writer0.writerow([lon, lat, depth, magnitude, r0, z0, z])

    # Create polydata object
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    
    # Add arrays to point data
    polydata.GetPointData().AddArray(diameter_array)
    polydata.GetPointData().AddArray(magnitude_array)
    polydata.GetPointData().AddArray(depth_array)

    # Create sphere source for glyphs
    sphere_source = vtk.vtkSphereSource()
    sphere_source.SetPhiResolution(10)
    sphere_source.SetThetaResolution(10)

    # Use vtkGlyph3D to place spheres at points
    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(polydata)
    glyph.SetSourceConnection(sphere_source.GetOutputPort())
    # Set the array to use for scaling (diameter)
    glyph.SetScaleModeToScaleByScalar()
    glyph.SetInputArrayToProcess(
        0,  # idx
        0,  # port
        0,  # connection
        vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,
        "diameter"  # The name of the array to use for scaling
    )
    glyph.SetScaleFactor(1.0)
    glyph.Update()

    # Write data to '.VTK' file
    writer1 = vtk.vtkXMLPolyDataWriter()
    output_VTK_filename = os.path.join(
        output_files_path, f"{mode}_catalog_minMAG{CONST.MIN_MAGNITUDE}_scatter_{catalog_filename.replace('csv', 'vtu')}")
    writer1.SetFileName(output_VTK_filename)
    writer1.SetInputConnection(glyph.GetOutputPort())
    writer1.SetDataModeToBinary()
    writer1.Write()

    # Define magnitude ranges and corresponding colors
    color_map = [
        [1.0, 4.0, (1.0, 0.0, 0.498039)],       # 1.0 <= mag < 4.0
        [4.0, 5.0, (0.0, 1.0, 1.0)],            # 4.0 <= mag < 5.0
        [5.0, 6.0, (0.666667, 0.666667, 1.0)],  # 5.0 <= mag < 6.0
        [6.0, 7.0, (0.0, 0.0, 1.0)],            # 6.0 <= mag < 7.0
        [7.0, max_magnitude if max_magnitude > 8.0 else 8.0, (1.0, 0.0, 0.0)],  # 7.0 <= mag <= max_magnitude
    ]

    # Create a Discrete Color Transfer Function
    color_transfer_func = vtk.vtkDiscretizableColorTransferFunction()
    color_transfer_func.DiscretizeOn()  # Enable discretization
    color_transfer_func.SetNumberOfValues(len(color_map))
    
    # Add color segments and annotations to the color transfer function
    mag_mid = 0.25
    if max_magnitude <= 8.0:
        d_label = (  8.0-1.0  ) / len(color_map)
    else:
        d_label = (  max_magnitude-1.0  ) / len(color_map)
    for mag_min, mag_max, color in color_map:
        # Add a segment from mag_min to mag_max with the specified color
        color_transfer_func.AddRGBSegment(mag_min, *color, mag_max, *color)
        # Calculate the midpoint for annotation
        mag_mid += d_label
        label = f"{mag_min}~{mag_max}"
        color_transfer_func.SetAnnotation(mag_mid, label)

    # Create mapper and actor for visualization
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.SetLookupTable(color_transfer_func)
    # Set the array to use for coloring (magnitude)
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("magnitude")
    mapper.SetColorModeToMapScalars()
    mapper.SetScalarRange(1.0, max_magnitude if max_magnitude > 8.0 else 8.0)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create color bar (Scalar Bar)
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_transfer_func)
    scalar_bar.SetTitle("Magnitude")
    # Position labels on the right side of the color bar
    scalar_bar.SetTextPositionToSucceedScalarBar()  # Labels on the right
    # Enable annotations and disable tick labels
    scalar_bar.DrawTickLabelsOff()    # Disable tick labels
    scalar_bar.DrawAnnotationsOn()    # Enable annotations
    # Adjust annotation properties
    label_text_property = scalar_bar.GetAnnotationTextProperty()
    label_text_property.SetFontSize(30)
    label_text_property.SetColor(0, 0, 0)  # Black color
    label_text_property.SetJustificationToCentered()
    label_text_property.SetVerticalJustificationToCentered()
    # Adjust title properties
    title_text_property = scalar_bar.GetTitleTextProperty()
    title_text_property.SetFontSize(40)
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
    render_window.SetWindowName(f'{mode} Earthquake Catalog Scatter Visualization')
    render_window.SetSize(1600, 900)  # Set render window size
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
                output_files_path, f"{mode}_catalog_minMAG{CONST.MIN_MAGNITUDE}_scatter_{catalog_filename.replace('csv', 'png')}")
            # Create a vtkPNGWriter object
            writer2 = vtk.vtkPNGWriter()
            writer2.SetFileName(output_PNG_filename)
            writer2.SetInputConnection(window_to_image_filter.GetOutputPort())
            writer2.Write()

            print(f"Screenshot saved as {output_PNG_filename}.")

    # Add key press event observer to the interactor
    render_window_interactor.AddObserver("KeyPressEvent", save_screenshot)

    # Start rendering and interaction
    render_window.Render()
    print("Press 'h' to save a screenshot.")
    render_window_interactor.Start()

if __name__ == "__main__":
    main()
