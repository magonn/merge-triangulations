#!/usr/bin/env python
# -*- coding: utf-8 -*-

from triangulation import *

if __name__ == '__main__':
    app = ConstructTriangulation(5000, 5000)
    app.ExperimentTime()
    app.root.mainloop()
