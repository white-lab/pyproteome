from unittest import TestCase

from pyproteome import pride


class PrideTest(TestCase):
    def test_list_data(self):
        # Raven Data
        ds = pride.list_data_set("PXD003660")

        self.assertEqual(
            len(ds),
            60,
        )
        self.assertIn(
            "20140524_MCF10A_E20VR1_ETP_TMT10.raw",
            [i.get("name") for i in ds],
        )

        # Kristina Data
        ds = pride.list_data_set("PXD006114")

        self.assertEqual(
            len(ds),
            42,
        )
        self.assertIn(
            "20151129_H1975HGF_VehicleControl_TMT6.mgf",
            [i.get("name") for i in ds],
        )
