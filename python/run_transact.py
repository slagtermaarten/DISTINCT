import sys
import pandas as pd
from os.path import join
from numpy import abs

from transact.pv_computation import PVComputation
from transact.interpolation import Interpolation
from transact.kernel_computer import KernelComputer
from transact.TRANSACT import TRANSACT

tmp_dir = sys.argv[1]
n_jobs = int(sys.argv[2])
kernel = sys.argv[3]
log_gamma = float(sys.argv[4])
N_source_components = int(sys.argv[5])
N_target_components = int(sys.argv[6])
N_PV = int(sys.argv[7])
center_ref = bool(sys.argv[8])
center_query = bool(sys.argv[9])
id_string = sys.argv[10]

# print(id_string)
# id_string = '{}_{.2f}_{}_{}_{b}_{b}'.format(kernel, log_gamma, N_source_components, N_target_components, center_ref, center_query)

## These files have been written by the wrapping R function
reference_M = pd.read_csv(join(tmp_dir, 'reference.csv'))
query_M = pd.read_csv(join(tmp_dir, 'query.csv'))

try:

  # if kernel in ['linear', 'cosine']:
  if kernel == 'poly':
    kernel_params = { 'degree': abs(log_gamma) }
  elif kernel in ['rbf', 'sigmoid', 'mallow', 'mallows']:
    kernel_params = { 'gamma': 10**log_gamma }
  else:
    kernel_params = {}

  # kernel_params = {'gamma': 10**log_gamma}

  TRANSACT_clf = TRANSACT(
    kernel = kernel,
    kernel_params = kernel_params,
    n_jobs = n_jobs, 
    verbose = True
  )

  TRANSACT_clf.fit(
    reference_M, query_M, 
    n_components = {
      'source' : N_source_components, 
      'target' : N_target_components
    },
    n_pv = N_PV
  )

  reference_CF = pd.DataFrame(TRANSACT_clf.transform(reference_M, center = center_ref))
  reference_CF.to_csv(join(tmp_dir, 'reference_CF_{}.csv'.format(id_string)), index=False)

  query_CF = pd.DataFrame(TRANSACT_clf.transform(query_M, center = center_query))
  query_CF.to_csv(join(tmp_dir, 'query_CF_{}.csv'.format(id_string)), index=False)

  CS_mat = pd.DataFrame(TRANSACT_clf.principal_vectors_.cosine_similarity_)
  CS_mat.to_csv(join(tmp_dir, 'CS_mat_{}.csv'.format(id_string)), index=False)

  Kst_mat = pd.DataFrame(TRANSACT_clf.interpolation_.kernel_values_.k_st)
  Kst_mat.to_csv(join(tmp_dir, 'Kst_mat_{}.csv'.format(id_string)), index=False)

  Ks_mat = pd.DataFrame(TRANSACT_clf.interpolation_.kernel_values_.k_s)
  Ks_mat.to_csv(join(tmp_dir, 'Ks_mat_{}.csv'.format(id_string)), index=False)

  Kt_mat = pd.DataFrame(TRANSACT_clf.interpolation_.kernel_values_.k_t)
  Kt_mat.to_csv(join(tmp_dir, 'Kt_mat_{}.csv'.format(id_string)), index=False)

  EV_mat = pd.DataFrame(TRANSACT_clf.principal_vectors_.dim_reduc_clf_['source'].eigenvectors_)
  EV_mat.to_csv(join(tmp_dir, 'EV_mat_{}.csv'.format(id_string)), index=False)

  alpha_source = pd.DataFrame(TRANSACT_clf.principal_vectors_.alpha_coef['source'])
  alpha_source.to_csv(join(tmp_dir, 'alpha_source_{}.csv'.format(id_string)), index=False)

  from numpy import zeros
  PV_projected_source = pd.DataFrame(TRANSACT_clf.interpolation_.transform(X=reference_M, tau=zeros(N_PV), center=False))
  PV_projected_source.to_csv(join(tmp_dir, 'PV_projected_source_{}.csv'.format(id_string)), index=False)

except Exception as e:
  print('Could not compute results')
  print(e)
