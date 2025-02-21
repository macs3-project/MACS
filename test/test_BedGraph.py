# Time-stamp: <2025-02-10 08:46:48 Tao Liu>

import pytest

from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.IO.PeakIO import PeakIO

test_regions1 = [(b"chrY", 0, 10593155, 0.0),
                 (b"chrY", 10593155, 10597655, 0.0066254580149)]
overlie_max = [(b"chrY", 0, 75, 20.0),
               (b"chrY", 75, 85, 35.0),
               (b"chrY", 85, 90, 75.0),
               (b"chrY", 90, 150, 10.0),
               (b"chrY", 150, 155, 0.0)]
overlie_mean = [(b"chrY", 0, 70, 6.66667),
                (b"chrY", 70, 75, 9.0),
                (b"chrY", 75, 80, 14.0),
                (b"chrY", 80, 85, 11.66667),
                (b"chrY", 85, 90, 36.66667),
                (b"chrY", 90, 150, 3.33333),
                (b"chrY", 150, 155, 0.0)]
overlie_fisher = [(b"chrY", 0, 70, (0, 0, 20), 92.10340371976183, 1.1074313239555578e-17, 16.9557),
                  (b"chrY", 70, 75, (7, 0, 20), 124.33959502167846, 1.9957116587802055e-24, 23.6999),
                  (b"chrY", 75, 80, (7, 0, 35), 193.41714781149986, 4.773982707347631e-39, 38.3211),
                  (b"chrY", 80, 85, (0, 0, 35), 161.1809565095832, 3.329003070922764e-32, 31.4777),
                  (b"chrY", 85, 90, (0, 75, 35), 506.56872045869005, 3.233076792862357e-106, 105.4904),
                  (b"chrY", 90, 150, (0, 0, 10), 46.051701859880914, 2.8912075645386016e-08, 7.5389),
                  (b"chrY", 150, 155, (0, 0, 0), 0.0, 1.0, 0.0)]
some_peak = [(b"chrY", 50, 80),
             (b"chrY", 95, 152)]


@pytest.fixture
def define_regions():
    test_regions1 = [(b"chrY", 0, 70, 0.0),
                     (b"chrY", 70, 80, 7.0),
                     (b"chrY", 80, 150, 0.0)]
    test_regions2 = [(b"chrY", 0, 85, 0.0),
                     (b"chrY", 85, 90, 75.0),
                     (b"chrY", 90, 155, 0.0)]
    test_regions3 = [(b"chrY", 0, 75, 20.0),
                     (b"chrY", 75, 90, 35.0),
                     (b"chrY", 90, 150, 10.0)]
    bdg1 = bedGraphTrackI()
    bdg2 = bedGraphTrackI()
    bdg3 = bedGraphTrackI()
    for a in test_regions1:
        bdg1.add_loc(a[0], a[1], a[2], a[3])
    for a in test_regions2:
        bdg2.add_loc(a[0], a[1], a[2], a[3])
    for a in test_regions3:
        bdg3.add_loc(a[0], a[1], a[2], a[3])
    return (bdg1, bdg2, bdg3)


def test_add_loc1():
    bdg = bedGraphTrackI()
    for a in test_regions1:
        bdg.add_loc(a[0], a[1], a[2], a[3])


def test_refine_peak():
    bdg = bedGraphTrackI()
    for a in overlie_mean:
        bdg.add_loc(a[0], a[1], a[2], a[3])
    peak = PeakIO()
    for a in some_peak:
        peak.add(a[0], a[1], a[2])
    new_peak = bdg.refine_peaks(peak)
    out = str(new_peak)
    std = "chrom:chrY\tstart:50\tend:80\tname:peak_1\tscore:14\tsummit:77\nchrom:chrY\tstart:95\tend:152\tname:peak_2\tscore:3.33333\tsummit:122\n"
    assert out == std
