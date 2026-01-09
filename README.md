# dnatco-tests
Automatic pytests for libLLKA and DNATCO (https://github.com/cernylab/libLLKA , https://github.com/cernylab/dnatco ).


# Recommended usage:
```bash
ccp4-python -m pytest -s -vv
```

The environmental variable `$DNATCO_ASSETS_PATH` should be set so DNATCO could find path with required CSV files. The variable can be set manually:

```bash
DNATCO_ASSETS_PATH="/xtal/ccp4-10/share/dnatco" ccp4-python -m pytest -s -vv
```

# Tests description

`test_classify_and_write_cif_integration.py` uses the following structures:

 * PDB 3A3A - crystal structure of human selenocystine tRNA, contains insertion codes, low resolution
 * PDB 5JZQ - centrosymmetric crystal structure of Z-DNA, contains alternate conformations, ultra high resolution
 * PDB 9BKD - large structure from SPA cryoEM - human Pdcd4 bound to the 40S small ribosomal subunit


