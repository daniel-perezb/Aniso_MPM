import cv2
import numpy as np
import open3d as o3d
import trimesh
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
import os
import configparser

def load_and_preprocess_mask(mask_path):
    """
    Load the mask image, normalize, and apply preprocessing steps.
    """
    # Load mask as grayscale
    mask = cv2.imread(mask_path, cv2.IMREAD_GRAYSCALE)
    if mask is None:
        raise FileNotFoundError(f"Mask image '{mask_path}' not found.")
    
    # Normalize to [0,1]
    mask = mask.astype(np.float32) / 255.0
    
    # Gaussian Blur 
    mask = cv2.GaussianBlur(mask, (5,5), 0)
    
    # Threshold to binarize
    _, mask = cv2.threshold(mask, 0.5, 1.0, cv2.THRESH_BINARY)
    
    # Remove noise
    kernel = np.ones((3,3), np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
    
    return mask

def extract_contours(mask):
    """
    Extract contours from the binary mask.
    """
    # Convert mask to uint8
    mask_uint8 = (mask * 255).astype(np.uint8)
    
    # Find contours using OpenCV
    contours, hierarchy = cv2.findContours(mask_uint8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Convert contours to list of (x, y) tuples
    polygons = []
    for cnt in contours:
        cnt = cnt.squeeze()
        if len(cnt.shape) < 2:
            continue  # Invalid contour
        polygon = [(point[0], point[1]) for point in cnt]
        polygons.append(polygon)
    
    return polygons

def create_extruded_mesh(polygons, extrusion_height):
    """
    Extrude the 2D polygons into a 3D mesh using Trimesh.
    """
    meshes = []
    for polygon in polygons:
        # Create Shapely polygon
        shapely_polygon = Polygon(polygon)
        if not shapely_polygon.is_valid:
            # Attempt to fix invalid polygons
            shapely_polygon = shapely_polygon.buffer(0)
            if not shapely_polygon.is_valid:
                print("Invalid polygon encountered and couldn't be fixed. Skipping.")
                continue
        if shapely_polygon.is_empty:
            print("Empty polygon encountered. Skipping.")
            continue
        
        # Handle MultiPolygons
        if isinstance(shapely_polygon, MultiPolygon):
            for poly in shapely_polygon:
                if poly.is_valid and not poly.is_empty:
                    mesh = trimesh.creation.extrude_polygon(poly, extrusion_height)
                    meshes.append(mesh)
        else:
            # Extrude polygon to create 3D mesh
            mesh = trimesh.creation.extrude_polygon(shapely_polygon, extrusion_height)
            meshes.append(mesh)
    
    if not meshes:
        raise ValueError("No valid meshes were created from the polygons.")
    
    # Combine all meshes into a single mesh
    combined_mesh = trimesh.util.concatenate(meshes)
    
    return combined_mesh

def convert_trimesh_to_open3d(trimesh_mesh):
    """
    Convert a Trimesh mesh to an Open3D mesh.
    """
    o3d_mesh = o3d.geometry.TriangleMesh()
    o3d_mesh.vertices = o3d.utility.Vector3dVector(trimesh_mesh.vertices)
    o3d_mesh.triangles = o3d.utility.Vector3iVector(trimesh_mesh.faces)
    o3d_mesh.compute_vertex_normals()
    return o3d_mesh

def load_parameters_from_file(parameter_file):
    """
    Load parameters from an external file.
    """
    config = configparser.ConfigParser()
    config.read(parameter_file)

    camera_distance = float(config['CameraParameters']['camera_distance'])
    focal_length_mm = float(config['CameraParameters']['focal_length_mm'])
    sensor_width_mm = float(config['CameraParameters']['sensor_width_mm'])
    sensor_height_mm = float(config['CameraParameters']['sensor_height_mm'])

    return camera_distance, focal_length_mm, sensor_width_mm, sensor_height_mm

def main():
    # Load parameters from file
    parameter_file = "../../Data/TetMesh/fleece_files/camera_params.txt"
    camera_distance, focal_length_mm, sensor_width_mm, sensor_height_mm = load_parameters_from_file(parameter_file)

    # Paths to your images
    mask_path = '../../Data/TetMesh/fleece_files/intial_mask.png'  # Replace with your mask path
    
    # Load and preprocess mask
    mask = load_and_preprocess_mask(mask_path)
    print("Mask loaded and preprocessed.")
    
    # Extract contours
    polygons = extract_contours(mask)
    print(f"Extracted {len(polygons)} polygon(s) from the mask.")
    
    if not polygons:
        print("No polygons found. Exiting.")
        return
    
    # Define extrusion height
    extrusion_height = 50  # Extrusion height
    
    # Create extruded mesh using Trimesh
    try:
        extruded_trimesh = create_extruded_mesh(polygons, extrusion_height)
        print("Extrusion completed using Trimesh.")
    except ValueError as e:
        print(f"Error during extrusion: {e}")
        return
    
    # Convert Trimesh mesh to Open3D mesh
    o3d_mesh = convert_trimesh_to_open3d(extruded_trimesh)
    print("Converted Trimesh mesh to Open3D mesh.")
    
    # Get image resolution in pixels
    image_height, image_width = mask.shape  
    
    # Convert focal length to meters
    focal_length_m = focal_length_mm / 1000.0  # Focal length in meters
    
    # Calculate real-world width and height based on camera distance and sensor size
    real_width = (camera_distance * sensor_width_mm) / focal_length_mm  
    real_height = (camera_distance * sensor_height_mm) / focal_length_mm  
    thickness = 0.05  # Desired thickness of the object in meters 
    
    # Scaling factors
    voxel_size_x = real_width / image_width 
    voxel_size_y = real_height / image_height 
    voxel_size_z = thickness / extrusion_height
 
    # Apply scaling to the mesh
    print("Applying scaling to the mesh...")
    vertices = np.asarray(o3d_mesh.vertices)
    vertices[:, 0] *= voxel_size_x
    vertices[:, 1] *= voxel_size_y
    vertices[:, 2] *= voxel_size_z
    o3d_mesh.vertices = o3d.utility.Vector3dVector(vertices)
    
    # Center the mesh at the origin
    print("Centering the mesh at the origin...")
    o3d_mesh.translate(-o3d_mesh.get_center())
    
    # Define rotation angles in degrees
    rotation_deg_x = -90.0
    rotation_deg_y = 180.0
    rotation_deg_z = 0.0

    # Convert rotation angles from degrees to radians
    rotation_rad_x = np.deg2rad(rotation_deg_x)
    rotation_rad_y = np.deg2rad(rotation_deg_y)
    rotation_rad_z = np.deg2rad(rotation_deg_z)
    
    # Create combined rotation matrix (Z -> Y -> X rotation)
    combined_rotation = o3d.geometry.get_rotation_matrix_from_xyz(
        [rotation_rad_x, rotation_rad_y, rotation_rad_z]
    )
    
    # Rotate the mesh
    print("Applying rotation to the mesh...")
    o3d_mesh.rotate(combined_rotation, center=(0, 0, 0))
    
    # Define target centroid coordinates (x, y, z) 
    # This is becauce of the centroid in the fleece.h for simulation
    target_centroid_x = 2.1
    target_centroid_y = 1.832
    target_centroid_z = 2.0
    target_centroid = np.array([target_centroid_x, target_centroid_y, target_centroid_z])
    
    # Compute current centroid after rotation
    current_centroid = o3d_mesh.get_center()
    print(f"Current centroid after rotation: {current_centroid}")
    
    # Compute translation vector
    translation_vector = target_centroid - current_centroid
    print(f"Translation vector: {translation_vector}")
 
    # Apply translation to the mesh
    print("Applying translation to the mesh...")
    o3d_mesh.translate(translation_vector)
    
    # Save the transformed mesh as a PLY file
    output_file = "../../Data/TetMesh/fleece_files/fleece_mesh.ply"
    success = o3d.io.write_triangle_mesh(output_file, o3d_mesh)
    if success:
        print(f"Transformed mesh saved as '{output_file}'")
    else:
        print(f"Failed to save the transformed mesh as '{output_file}'.")
 
    print("Process completed")

if __name__ == "__main__":
    main()

