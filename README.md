Mutual Information between Continuous Variables in Python
====

This is an implementation of [a B-Spline-based algorithm](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-118). There is a nice R version developed by [Weiyi Cheng](https://github.com/weiyi-bitw/cafr).

## Usage ##

Have `mutualinfo` in the working folder.

```r
import mutualinfo as m
m.mi2D([1,2,7], [0,0,1], 6, 3)
```

`x` and `y` must be numpy arrays in `float32`

## Prerequisites ##

Python 2.7 and 3 above 

NumPy

## License ##

This project is licensed under the MIT license.
