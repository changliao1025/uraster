"""
Basic tests for uraster package
"""

import unittest
import sys
import os

# Add the parent directory to the path so we can import uraster
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from uraster import uraster, sraster
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure uraster package is properly installed or in PYTHONPATH")


class TestUrasterImport(unittest.TestCase):
    """Test basic import functionality"""

    def test_import_uraster(self):
        """Test that uraster can be imported"""
        from uraster import uraster
        self.assertTrue(uraster is not None)

    def test_import_sraster(self):
        """Test that sraster can be imported"""
        from uraster import sraster
        self.assertTrue(sraster is not None)

    def test_uraster_instantiation(self):
        """Test that uraster can be instantiated"""
        try:
            u = uraster()
            self.assertTrue(u is not None)
        except Exception as e:
            self.fail(f"uraster instantiation failed: {e}")

    def test_sraster_instantiation(self):
        """Test that sraster can be instantiated"""
        try:
            s = sraster()
            self.assertTrue(s is not None)
        except Exception as e:
            self.fail(f"sraster instantiation failed: {e}")


if __name__ == '__main__':
    unittest.main()