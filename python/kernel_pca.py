from sklearn.decomposition import KernelPCA

def run_kPCA(kernel_file, n_jobs, gamma, N_NLPC):
  import numpy as np
  import pandas as pd
  # from transact.matrix_operations import _center_kernel

  kernel = pd.read_csv(kernel_file)

  # K = _center_kernel(K)
  # K = _right_center_kernel(K)

  dim_reduc_clf = KernelPCA(
    n_components=N_NLPC,
    eigen_solver='dense',
    tol=1e-15,
    max_iter=1e6,
    kernel='precomputed',
    n_jobs=1,
    gamma=gamma,
    fit_inverse_transform=False
  )
  dim_reduc_clf.fit(kernel)

  if False:
    ev_norms = np.linalg.norm(dim_reduc_clf.eigenvectors_, axis=0)
    if np.all(ev_norms <= 1.01):
      print('Well-behaved!')
    else:
      print('Not well-behaved!')
    print(ev_norms)
    print(n_components[t])
    # assert np.all(ev_norms <= 1.01), 'Non-normed eigenvectors detected'

  return dim_reduc_clf
  # return dim_reduc_clf.eigenvectors_
  # alpha_coef = dim_reduc_clf.eigenvectors_ / np.sqrt(dim_reduc_clf.eigenvalues_)
  # return alpha_coef


def run_kPCA_data(M, gamma, N_NLPC, 
    train_idxs=None, 
    kernel = 'mallow', n_jobs=1):
  import numpy as np
  import pandas as pd
  from transact.alternative_kernels import mallow_kernel_wrapper
  from transact.matrix_operations import _center_kernel

  if train_idxs == None:
    train_idxs = np.arange(0, M.shape[0]-1)

  print(len(train_idxs))
  if len(train_idxs) < M.shape[0]:
    apply_idxs = np.arange(len(train_idxs), M.shape[0])
  else:
    apply_idxs = None

  if kernel == 'mallow' or kernel == 'mallows':
    kernel_ = mallow_kernel_wrapper(n_jobs)

    K = kernel_(M[train_idxs,:], gamma = gamma)
    K = _center_kernel(K)

    dim_reduc_clf = KernelPCA(
      n_components=N_NLPC,
      eigen_solver='dense',
      tol=1e-15,
      max_iter=1e6,
      kernel='precomputed',
      n_jobs=n_jobs,
      gamma=gamma,
      fit_inverse_transform=False
    )
    sample_scores = dim_reduc_clf.fit_transform(K)

    if len(apply_idxs) > 0:
      K_a = kernel_(M[apply_idxs,:], M[train_idxs,:], gamma = gamma)
      K_a = _center_kernel(K_a)
      apply_scores = dim_reduc_clf.transform(K_a)
      sample_scores = np.vstack((sample_scores, apply_scores))

  else:
    dim_reduc_clf = KernelPCA(
      n_components=N_NLPC,
      kernel=kernel, 
      fit_inverse_transform=False, 
      gamma=gamma,
      n_jobs=n_jobs,
      tol=1e-15,
      max_iter=1e6
    )
    # print(train_idxs)
    # print(apply_idxs)
    sample_scores = dim_reduc_clf.fit_transform(M[train_idxs,:])

    if len(apply_idxs) > 0:
      apply_scores = dim_reduc_clf.transform(M[apply_idxs,:])
      sample_scores = np.vstack((sample_scores, apply_scores))

  return sample_scores


# if __name__ == '__main__':
#   import sys
#   kernel_file = sys.argv[1]
#   n_jobs = int(sys.argv[2])
#   gamma = float(sys.argv[3])
#   N_NLPC = float(sys.argv[4])
#   run_kPCA(
#     kernel_file=kernel_file, 
#     n_jobs=n_jobs, gamma=gamma, N_NLPC=N_NLPC
#   )
