import os.path
import unittest
from result_catch import SL_result_catch as SL, generate_path, save_load, chkdir


class TestResultCatch(unittest.TestCase):
    def test_sl(self):
        print(os.path.abspath("../../"))
        data = {
            "arm_number": 4,
            "temp": 20,
            "conc": 0.1,
        }
        print(f"top_path: {generate_path(data, 'top_path')}, {os.path.isfile(generate_path(data, 'top_path'))}")
        print(f"top_path: {generate_path(data, 'traj_path')}, {os.path.isfile(generate_path(data, 'traj_path'))}")
