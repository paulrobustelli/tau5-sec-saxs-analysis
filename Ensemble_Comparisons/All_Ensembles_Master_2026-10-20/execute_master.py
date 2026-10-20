from pathlib import Path
import json,time
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter
O=Path(__file__).resolve().parent;p=O/'Tau5_All_Ensembles_Master.ipynb';n=nbformat.read(p,as_version=4)
start=time.time();NotebookClient(n,timeout=180,kernel_name='python3',resources={'metadata':{'path':str(O)}}).execute()
nbformat.write(n,p)
html,_=HTMLExporter().from_notebook_node(n);(O/'Tau5_All_Ensembles_Master.html').write_text(html)
errors=[o for c in n.cells if c.cell_type=='code' for o in c.outputs if o.output_type=='error'];assert not errors
(O/'execution_validation.json').write_text(json.dumps({'code_cells':sum(c.cell_type=='code' for c in n.cells),'errors':len(errors),'elapsed_seconds':time.time()-start,'html_bytes':len(html.encode())},indent=2));print('Executed and exported',flush=True)
