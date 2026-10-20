from pathlib import Path
import json,time,nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter
O=Path(__file__).resolve().parent;p=O/'Hydrophobic_SAXS_hypothesis.ipynb';n=nbformat.read(p,as_version=4);start=time.time()
NotebookClient(n,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(O)}}).execute();nbformat.write(n,p)
body,_=HTMLExporter().from_notebook_node(n);(O/'Hydrophobic_SAXS_hypothesis.html').write_text(body)
errors=[o for c in n.cells if c.cell_type=='code' for o in c.outputs if o.output_type=='error'];assert not errors
(O/'notebook_execution.json').write_text(json.dumps({'code_cells':sum(c.cell_type=='code' for c in n.cells),'errors':len(errors),'seconds':time.time()-start},indent=2));print('Executed successfully')
