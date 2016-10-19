Mutual Information between Continuous Variables in Python
====

An implementation of [a B-Spline-based algorithm](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-118). A nice R version is developed by [Weiyi Cheng](https://github.com/weiyi-bitw/cafr).

## Usage ##

```r
import numpy as np
import mutualinfo as m
m.mi2D(np.linspace(0,1,10), np.linspace(0,1,10), 6, 3)
```

## Prerequisites ##

Python2.7 and NumPy

## License ##

This project is licensed under the MIT license.