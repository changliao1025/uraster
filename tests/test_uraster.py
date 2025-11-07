"""
Comprehensive tests for uraster package
"""

import unittest
import sys
import os
import tempfile
import numpy as np
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import uraster
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from uraster.classes.uraster import uraster
    from uraster.classes.sraster import sraster
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure uraster package is properly installed or in PYTHONPATH")


class TestUrasterConfiguration(unittest.TestCase):
    """Test uraster configuration and initialization"""

    def test_default_initialization(self):
        """Test uraster initialization with default config"""
        u = uraster()
        self.assertIsNotNone(u)
        self.assertEqual(u.iFlag_remap_method, 1)
        self.assertIsNone(u.sFilename_source_mesh)
        self.assertEqual(u.aFilename_source_raster, [])

    def test_custom_configuration(self):
        """Test uraster initialization with custom config"""
        config = {
            'iFlag_remap_method': 2,
            'sFilename_source_mesh': 'test_mesh.geojson',
            'aFilename_source_raster': ['test_raster.tif']
        }
        u = uraster(config)
        self.assertEqual(u.iFlag_remap_method, 2)
        self.assertEqual(u.sFilename_source_mesh, 'test_mesh.geojson')
        self.assertEqual(u.aFilename_source_raster, ['test_raster.tif'])

    def test_invalid_remap_method(self):
        """Test handling of invalid remap method"""
        config = {'iFlag_remap_method': 99}
        u = uraster(config)
        # Should default to 1 for invalid method
        self.assertEqual(u.iFlag_remap_method, 1)


class TestSrasterConfiguration(unittest.TestCase):
    """Test sraster configuration and initialization"""

    def test_sraster_initialization_no_file(self):
        """Test sraster initialization without file"""
        s = sraster()
        self.assertIsNotNone(s)
        self.assertIsNone(s.sFilename)

    def test_sraster_initialization_nonexistent_file(self):
        """Test sraster initialization with nonexistent file"""
        with self.assertRaises(FileNotFoundError):
            sraster('nonexistent_file.tif')

    def test_sraster_mesh_filename_generation(self):
        """Test mesh filename generation"""
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
            tmp_path = tmp.name

        try:
            s = sraster(tmp_path)
            expected_mesh = tmp_path.replace('.tif', '_mesh.geojson')
            self.assertEqual(s.sFilename_mesh, expected_mesh)
        finally:
            os.unlink(tmp_path)


class TestUrasterValidation(unittest.TestCase):
    """Test uraster validation methods"""

    def setUp(self):
        self.uraster_instance = uraster()

    def test_check_raster_files_empty_list(self):
        """Test raster file validation with empty list"""
        result = self.uraster_instance.check_raster_files([])
        self.assertIsNone(result)

    def test_check_raster_files_invalid_type(self):
        """Test raster file validation with invalid type"""
        result = self.uraster_instance.check_raster_files("not_a_list")
        self.assertIsNone(result)

    def test_check_raster_files_nonexistent_file(self):
        """Test raster file validation with nonexistent file"""
        result = self.uraster_instance.check_raster_files(['nonexistent.tif'])
        self.assertIsNone(result)

    def test_check_mesh_file_no_filename(self):
        """Test mesh file validation with no filename"""
        result = self.uraster_instance.check_mesh_file()
        self.assertIsNone(result)

    def test_check_mesh_file_nonexistent(self):
        """Test mesh file validation with nonexistent file"""
        self.uraster_instance.sFilename_source_mesh = 'nonexistent.geojson'
        result = self.uraster_instance.check_mesh_file()
        self.assertIsNone(result)


class TestUrasterUtilities(unittest.TestCase):
    """Test uraster utility methods"""

    def setUp(self):
        self.uraster_instance = uraster()

    def test_geometry_type_name_conversion(self):
        """Test geometry type name conversion"""
        # Mock ogr constants for testing
        with patch('uraster.classes.uraster.ogr') as mock_ogr:
            mock_ogr.wkbPolygon = 3
            mock_ogr.wkbPoint = 1
            mock_ogr.wkbUnknown = 0
            mock_ogr.wkbLineString = 2
            mock_ogr.wkbMultiPoint = 4
            mock_ogr.wkbMultiLineString = 5
            mock_ogr.wkbMultiPolygon = 6
            mock_ogr.wkbGeometryCollection = 7

            # Test known geometry types
            self.assertEqual(
                self.uraster_instance._get_geometry_type_name(3),
                "wkbPolygon"
            )
            self.assertEqual(
                self.uraster_instance._get_geometry_type_name(1),
                "wkbPoint"
            )

            # Test unknown geometry type
            result = self.uraster_instance._get_geometry_type_name(999)
            self.assertIn("Unknown geometry type", result)

    def test_cleanup_method(self):
        """Test cleanup method"""
        # Should not raise any exceptions
        self.uraster_instance.cleanup()

    def test_print_methods(self):
        """Test print methods don't crash"""
        # These methods should not crash even with empty data
        with patch('builtins.print'):
            self.uraster_instance.print_mesh_info()
            self.uraster_instance.print_raster_info()


class TestUrasterReporting(unittest.TestCase):
    """Test uraster reporting methods"""

    def setUp(self):
        self.uraster_instance = uraster()

    def test_report_inputs(self):
        """Test input reporting"""
        with patch.object(self.uraster_instance, 'print_raster_info') as mock_raster, \
             patch.object(self.uraster_instance, 'print_mesh_info') as mock_mesh:
            self.uraster_instance.report_inputs()
            mock_raster.assert_called_once()
            mock_mesh.assert_called_once()

    def test_report_inputs_with_gpu_info(self):
        """Test input reporting with GPU info"""
        with patch.object(self.uraster_instance, 'print_raster_info'), \
             patch.object(self.uraster_instance, 'print_mesh_info'), \
             patch('builtins.print'):
            # Test with geovista available
            with patch('uraster.classes.uraster.gvreport') as mock_gv:
                mock_gv.Report.return_value = "GPU Info"
                self.uraster_instance.report_inputs(iFlag_show_gpu_info=True)

    def test_report_outputs_existing_file(self):
        """Test output reporting with existing file"""
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp.write(b"test data")
            tmp_path = tmp.name

        try:
            with patch('uraster.classes.uraster.logger') as mock_logger:
                self.uraster_instance.report_outputs(tmp_path)
                mock_logger.info.assert_called()
        finally:
            os.unlink(tmp_path)

    def test_report_outputs_nonexistent_file(self):
        """Test output reporting with nonexistent file"""
        with patch('uraster.classes.uraster.logger') as mock_logger:
            self.uraster_instance.report_outputs('nonexistent.txt')
            mock_logger.warning.assert_called()


class TestUrasterVisualization(unittest.TestCase):
    """Test uraster visualization methods"""

    def setUp(self):
        self.uraster_instance = uraster()

    @patch('uraster.classes.uraster._visual')
    def test_visualize_source_mesh(self, mock_visual):
        """Test source mesh visualization"""
        mock_visual.visualize_source_mesh.return_value = True

        result = self.uraster_instance.visualize_source_mesh(
            sFilename_out='test.png',
            dLongitude_focus_in=0.0,
            dLatitude_focus_in=0.0
        )

        mock_visual.visualize_source_mesh.assert_called_once()
        self.assertTrue(result)

    @patch('uraster.classes.uraster._visual')
    def test_visualize_target_mesh(self, mock_visual):
        """Test target mesh visualization"""
        mock_visual.visualize_target_mesh.return_value = True

        result = self.uraster_instance.visualize_target_mesh(
            sVariable_in='mean',
            sFilename_out='test.png'
        )

        mock_visual.visualize_target_mesh.assert_called_once()
        self.assertTrue(result)

    @patch('uraster.classes.uraster._visual')
    def test_visualize_raster(self, mock_visual):
        """Test raster visualization"""
        mock_visual.visualize_raster.return_value = True

        result = self.uraster_instance.visualize_raster()

        mock_visual.visualize_raster.assert_called_once()
        self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()