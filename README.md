# StretchDyn
A molecular dynamics simulator that simply models stretch energy. Best suited for diatomic or linear molecules.

## Installation

```
conda create -n StretchDyn python=3.8
conda activate StretchDyn
pip install -e .
```

## Development

To run Black for style checking and MyPy for type checking, from the root of the repo run:

``` 
black stretchdyn
mypy stretchdyn
```
