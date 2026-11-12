"""
Basic tests for uraster package
"""

import os
import sys
import unittest

# Add the parent directory to the path so we can import uraster
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from uraster import SRaster, URaster
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure uraster package is properly installed or in PYTHONPATH")


class TestUrasterImport(unittest.TestCase):
    """Test basic import functionality"""

    def test_import_uraster(self):
        """Test that URaster can be imported"""
        from uraster import URaster

        self.assertTrue(URaster is not None)

    def test_import_sraster(self):
        """Test that SRaster can be imported"""
        from uraster import SRaster

        self.assertTrue(SRaster is not None)

    def test_uraster_instantiation(self):
        """Test that URaster can be instantiated"""
        try:
            u = URaster()
            self.assertTrue(u is not None)
        except Exception as e:
            self.fail(f"URaster instantiation failed: {e}")

    def test_sraster_instantiation(self):
        """Test that SRaster can be instantiated"""
        try:
            s = SRaster()
            self.assertTrue(s is not None)
        except Exception as e:
            self.fail(f"SRaster instantiation failed: {e}")


if __name__ == "__main__":
    unittest.main()
