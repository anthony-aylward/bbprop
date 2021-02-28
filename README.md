# bbprop
Test for difference of beta-binomial proportions

## Installation

```sh
pip3 install bbprop
```
or
```
pip3 install --user bbprop
```

## Examples

```python
from bbprop import bbprop_cdf, bbprop_test
help(bbprop_cdf)
help(bbprop_test)
bbprop_cdf(0.1, [30, 12], [4, 5], [6, 5])
bbprop_test(0.1, [30, 12], [4, 5], [6, 5])
```
