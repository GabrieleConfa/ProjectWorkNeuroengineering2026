# Notebook generators

The block-study notebooks are written by these scripts, so that the long explanations and the code
stay in one plain Python file each (easy to diff and review). To regenerate a notebook, close it in
Jupyter without saving, then run from the project root:

```
python scripts/notebook_generators/make_volley_nb.py notebooks/block_volley.ipynb
python scripts/notebook_generators/make_design_nb.py notebooks/block_anode_design.ipynb
python scripts/notebook_generators/make_check_nb.py  notebooks/block_model_check.ipynb
```

Regenerating overwrites the notebook and clears its outputs. Edit the generator, not the notebook,
when the change should last.
