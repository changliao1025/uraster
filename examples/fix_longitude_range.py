"""
Script to fix longitude values in GeoJSON files to ensure they fall within [-180, 180] range.
Values outside this range are wrapped (e.g., -180.1 becomes 179.9, 180.1 becomes -179.9).
"""

import json

INPUT_GEOJSON = "data/example_1/input/rhealpix_global_res3.geojson"
OUTPUT_GEOJSON = "data/example_1/input/rhealpix_global_res3_fixed.geojson"


def normalize_longitude(lon):
    """
    Normalize longitude to [-180, 180] range by wrapping.
    
    Args:
        lon (float): Longitude value
        
    Returns:
        float: Normalized longitude in [-180, 180] range
    """
    while lon < -180:
        lon += 360
    while lon > 180:
        lon -= 360
    return lon


def fix_coordinates(coords):
    """
    Recursively fix longitude values in coordinate arrays.
    Handles nested coordinate structures for different geometry types.
    
    Args:
        coords: Coordinate array (can be nested)
        
    Returns:
        Fixed coordinate array with normalized longitudes
    """
    if not coords:
        return coords
    
    # Check if this is a coordinate pair [lon, lat]
    if isinstance(coords[0], (int, float)):
        return [normalize_longitude(coords[0]), coords[1]]
    else:
        return [fix_coordinates(coord) for coord in coords]


def fix_geometry(geometry):
    """
    Fix longitude values in a Polygon or MultiPolygon geometry object.
    
    Args:
        geometry (dict): GeoJSON geometry object (Polygon or MultiPolygon)
        
    Returns:
        dict: Geometry object with fixed coordinates
    """
    if geometry is None:
        return geometry
    
    coords = geometry.get("coordinates")
    if coords is None:
        return geometry
    
    fixed_geometry = geometry.copy()
    fixed_geometry["coordinates"] = fix_coordinates(coords)
    return fixed_geometry


def fix_geojson(input_file, output_file):
    """
    Read a GeoJSON file, fix all longitude values, and write to output file.
    
    Args:
        input_file (str): Path to input GeoJSON file
        output_file (str): Path to output GeoJSON file
    """
    print(f"Reading GeoJSON file: {input_file}")
    
    with open(input_file, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    print(f"Processing {len(geojson_data.get('features', []))} features...")
    
    # Fix coordinates in all features
    fixed_features = []
    for idx, feature in enumerate(geojson_data.get("features", [])):
        if idx % 1000 == 0 and idx > 0:
            print(f"  Processed {idx} features...")
        
        fixed_feature = feature.copy()
        if "geometry" in feature:
            fixed_feature["geometry"] = fix_geometry(feature["geometry"])
        fixed_features.append(fixed_feature)
    
    geojson_data["features"] = fixed_features
    
    print(f"Writing fixed GeoJSON to: {output_file}")
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(geojson_data, f, indent=2)
    
    print("Done!")


# if __name__ == "__main__":
#     fix_geojson(INPUT_GEOJSON, OUTPUT_GEOJSON)

