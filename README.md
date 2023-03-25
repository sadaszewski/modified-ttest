```diff
- NOTE: Not tested extensively. A few manual checks against outputs
- of the SpatialPack implementation checked out perfectly though.
```

# modified-ttest
Modified t-test reimplementation from R\'s SpatialPack (paper by Dutilleul)

This implementation relies almost entirely on the C code from SpatialPack: 
https://github.com/faosorios/SpatialPack/blob/master/src/mod_ttest.c

Kudos to Felipe Osorio.

## References

- Clifford, P., Richardson, S., Hemon, D. (1989). Assessing the significance 
   of the correlation between two spatial processes. Biometrics 45, 123--134.
- Dutilleul, P. (1993). Modifying the t-test for assessing the correlation
  between two spatial processes. Biometrics 49, 305--314.
