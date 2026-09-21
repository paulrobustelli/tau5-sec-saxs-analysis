from pathlib import Path
import nbformat,json,sys,time,os
from nbclient import NotebookClient
ROOT=Path(__file__).resolve().parents[1]
os.environ['IPYTHONDIR']=str(ROOT/'.ipython')
os.environ['JUPYTER_PATH']=str(ROOT/'.jupyter/share/jupyter')
for path in ([ROOT/name for name in sys.argv[1:]] if len(sys.argv)>1 else sorted(ROOT.glob('0*.ipynb'))):
    print('EXECUTING',path.name,flush=True);t=time.time()
    n=nbformat.read(path,as_version=4)
    client=NotebookClient(n,kernel_name='saxs',timeout=1800,startup_timeout=30,resources={'metadata':{'path':str(ROOT)}})
    client.execute()
    nbformat.write(n,path)
    print('FINISHED',path.name,round(time.time()-t,1),'seconds',flush=True)
