# StretchDyn
A molecular dynamics simulator that simply models stretch energy. Best suited for diatomic or linear molecules.

## Installation

```
conda create -n StretchDyn python=3.8
conda activate StretchDyn
pip install -e .
```

## Development

To run Black for style checking, MyPy for type checking, and pytest for testing run the following commands from the root of the repo:

``` 
black stretchdyn
mypy stretchdyn
pytest stretchdyn/tests
```
