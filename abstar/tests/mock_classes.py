#!usr/env/python
# filename: mock_classes.py

from ..utils.mixins import LoggingMixin


class MockAntibody(LoggingMixin):
    def __init__(self):
        super(MockAntibody, self).__init__()
        LoggingMixin.__init__(self)
