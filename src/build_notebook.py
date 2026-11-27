import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

here = Path(__file__).resolve().parent
source = "\n".join(p.read_text(encoding="utf-8") for p in sorted(here.glob("notebook_part*.py")))
cells = [block.strip("\n") for block in source.split("#%%\n") if block.strip()]
notebook = nbformat.v4.new_notebook()
notebook.cells = [nbformat.v4.new_code_cell(cell) for cell in cells]
notebook.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
notebook.metadata["colab"] = {"provenance": []}
target = here.parent / "notebooks" / "Risetku_ParacetamolAnalysis.ipynb"
if "--run" in sys.argv:
    NotebookClient(notebook, timeout=1800, kernel_name="python3", resources={"metadata": {"path": str(here)}}).execute()
nbformat.write(notebook, target)
print(target, len(cells), "cells")
