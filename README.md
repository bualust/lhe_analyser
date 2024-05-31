# LHE file analyser

A simple events analyser based on pylhe to run over LHE files.
In the current configuration it compares different LHE files over several distributions and saves them in the same output directory as png figures.

## To install and run

```bash
source setup.sh
pip install .
lhe-analyser --o output-plots-test --i lhe-files/ --c config.yaml
```
