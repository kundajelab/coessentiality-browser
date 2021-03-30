# coessentiality-browser
Gene browser using coessentiality and related data.
From the paper: 
> A genome-wide almanac of co-essential modules assigns function to uncharacterized genes. (https://www.biorxiv.org/content/10.1101/827071v1)


## Installation/dependencies

```
sudo apt-get install libblas-dev liblapack-dev gfortran
```

In a conda environment, run:
`pip install -r requirements.txt`.

The main requirement is Dash (https://github.com/plotly/dash/), a framework for building the interactive interface.


## Launching interactive browser
From the conda environment, run
`python app.py`

The default IP is 0.0.0.0:5000 .
