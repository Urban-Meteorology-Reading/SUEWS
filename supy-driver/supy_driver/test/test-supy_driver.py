from unittest import TestCase

import supy_driver as sd


# sd.suews_driver.output_name_n(3)[0]

class Test_supy_driver(TestCase):
    def test_is_string(self):
        # print(dir(sd))
        print(dir(sd.suews_driver))
        s = sd.suews_driver.output_name_n(3)[0]
        self.assertTrue(isinstance(s, bytes))
