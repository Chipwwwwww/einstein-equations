"""Numerical check of the pinned XCTS Hamiltonian diagnostic normalization."""
import h5py,numpy as np,json,sys
from pathlib import Path
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path('/workspace/scratch/spectre-pilot/id-smoke')
def mat(g,prefix):
 n=len(g['Lapse']);a=np.empty((n,3,3))
 for i,x in enumerate('xyz'):
  for j,y in enumerate('xyz'):a[:,i,j]=g[prefix+'_'+(x+y if i>=j else y+x)][:]
 return a
with h5py.File(root/'BbhVolume0.h5') as f:
 g=max((a for a in f['VolumeData.vol'].values() if isinstance(a,h5py.Group)),key=lambda a:a.attrs['observation_value']);H=g['HamiltonianConstraint'][:];psi=g['ConformalFactor'][:];gi=mat(g,'InverseSpatialMetric');K=mat(g,'ExtrinsicCurvature');R=mat(g,'SpatialRicci');trK=np.einsum('nij,nij->n',gi,K);KK=np.einsum('nik,njl,nij,nkl->n',gi,gi,K,K);C=np.einsum('nij,nij->n',gi,R)+trK**2-KK;converted=-2*H/psi**5
 out={'classification':'NUMERICAL with DOCUMENTED source definition','N':len(H),'raw_saved_RMS':float(np.sqrt(np.mean(H**2))),'raw_saved_max_abs':float(np.max(abs(H))),'physical_ADM_constraint_nodal_RMS':float(np.sqrt(np.mean(C**2))),'physical_ADM_constraint_max_abs':float(np.max(abs(C))),'rescaled_saved_RMS':float(np.sqrt(np.mean(converted**2))),'max_difference_C_ADM_minus_minus2H_over_psi5':float(np.max(abs(C-converted))),'continuum_definition':'H_saved = -psi^5/2 * (R + K^2 - K_ij K^ij - 16 pi rho); vacuum rho=0. L2Norm is unweighted nodal RMS, not a continuum/proper-volume L2 norm.'}
with h5py.File(root/'BbhReductions.h5') as f:out['reported_L2Norm']=float(f['ConstraintDiagnostics.dat'][-1,4])
print(json.dumps(out,indent=2))
