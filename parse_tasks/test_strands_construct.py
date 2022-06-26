import unittest
from .strands_construct import strands_construct


class TestParseTasks(unittest.TestCase):
    def test_strands_construct(self):
        ts = strands_construct(data={
            "arm_number": 4,
            "temp": 20,
            "conc": 0.1,
        })
        print(f"{ts.params}")
