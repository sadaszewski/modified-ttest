```diff
- NOTE: Not tested extensively. A few manual checks against outputs
- of the SpatialPack implementation turned out perfectly though.
```

# modified-ttest
Modified t-test reimplementation from R's SpatialPack (paper by Dutilleul)
for Python

This implementation relies almost entirely on the C code from SpatialPack: 
https://github.com/faosorios/SpatialPack/blob/master/src/mod_ttest.c

Kudos to Felipe Osorio.

## Example
```python
import modified_ttest
import numpy as np
from sklearn.datasets import load_breast_cancer
from scipy.stats import f

data = load_breast_cancer()['data']

m1 = np.corrcoef(data[np.random.permutation(len(data))[:5]], rowvar=False)
m2 = np.corrcoef(data[np.random.permutation(len(data))[:5]], rowvar=False)

x = m1.flatten()
y = m2.flatten()
rr, cc = np.meshgrid(range(m1.shape[0]), range(m1.shape[1]), indexing='ij')
coords = np.array([ rr.flatten(), cc.flatten() ]).T

ESS, F, df = modified_ttest.modified_ttest(x, y, coords)
pval = f.cdf(df * F, 1.0, df)

print('pval:', pval)
```

## References

- Clifford, P., Richardson, S., Hemon, D. (1989). Assessing the significance 
   of the correlation between two spatial processes. Biometrics 45, 123--134.
- Dutilleul, P. (1993). Modifying the t-test for assessing the correlation
  between two spatial processes. Biometrics 49, 305--314.
