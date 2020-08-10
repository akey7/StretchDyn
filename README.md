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

## Symbolic Differentiation with SymPy

[https://docs.sympy.org/latest/tutorial/calculus.html](https://docs.sympy.org/latest/tutorial/calculus.html)

IPython is great for this:

```
u, k_ab, r_ab, r_e_ab = symbols("u, k_ab, r_ab, r_e_ab") 
diff(0.5 * k_ab * (r_ab - r_e_ab) ** 2, r_ab) 
# 0.5*k_ab*(2*r_ab - 2*r_e_ab)
```
