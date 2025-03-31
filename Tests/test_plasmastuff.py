from plasmastuff import Sun
from pytest import raises

sun = Sun()

def test_addition():
    assert sun.add(1,3) == 4
    
def test_raises_error():
    with raises(AttributeError):
        sun.spots