"""PyWENO WENO reconstructor."""

import numpy as np
import pyweno.cweno

def reconstruct(q, k, points,
                n=None,
                smoothness=None,
                weights=None,
                return_smoothness=False,
                return_weights=False,
                squeeze=True,
                weno_type='JS'):
  """Perform WENO reconstruction of q.

  :param q:                 cell-averaged unknown to reconstruct
  :param k:                 order of reconstruction (odd values from 5 to 11)
  :param points:            reconstruction points
  :param n:                 number of reconstruction points (optional)
  :param smoothness:        use smoothness indicators *smoothness* (optional)
  :param weights:           use non-linear weights *weights* (optional)
  :param return_smoothness: return smoothness indicators? (default: ``False``)
  :param return_weights:    return weights? (default: ``False``)
  :param squeeze:           squeeze the results? (default: ``True``)
  :param weno_type:         Refers to the type of weno used, either 'JS' (Original),
                            'M', or 'Z'
 
  Supported reconstruction points *points* are:

  * ``'left'`` - left edge of each cell

  * ``'right'`` - right edge of each cell

  * ``'middle'`` - the middle of each cell

  * ``'+/-'`` - returns two reconstructions per cell, both at the left
    edge of each cell; the first is considered as a limit from the
    left (-), the second is considered as a limit from the right (+).

  * ``'gauss'`` or ``'gauss_legendre'`` - *n* point Gauss-Legendre
    quadrature points

  * ``'gauss_lobatto'`` - *n* point Gauss-Lobatto quadrature points

  * ``'gauss_radau'`` - *n* point Gauss-Radau quadrature points


  This function wraps several complied WENO kernels which were
  generated by PyWENO.  If you need WENO reconstructions at points not
  supported by this wrapper, please checkout `PyWENO
  <http://pyweno.readthedocs.org/>`_.

  """

  valid_k = (5, 11)
  valid_points = [ 'left', 'right', 'middle', '+/-',
                   'gauss', 'gauss_legendre',
                   'gauss_lobatto',
                   'gauss_radau' ]

  if (k % 2) == 0:
    raise ValueError('even order WENO reconstructions are not supported')

  if k < valid_k[0] or k > valid_k[1]:
    raise ValueError('%d order WENO reconstructions are not supported' % k)

  if not (points in valid_points):
    raise ValueError("'points' must be one of: " + ', '.join(valid_points))

  if points == 'gauss':
    points = 'gauss_legendre'

  N = q.shape[0]
  k = (k+1)/2

  # validate points and n
  if points in [ 'left', 'right', 'middle' ]:
    n = 1
  elif n is None:
    n = k

  if points in [ '+/-' ]:
    n = 2

  assert(n > 0)

  # smoothness
  if smoothness is None:
    smoothness = np.zeros((N,k))
    try:
      func = getattr(pyweno.cweno, 'smoothness%03d' % k)
    except AttributeError:
      raise ValueError('unsupported smoothness indicator: '
                       + 'k=%d not generated' % (2*k-1))
    func(q, smoothness)

  assert(smoothness.shape == (N,k))

  # weights
  if weights is None:
    weights = np.zeros((N,n,k))
    try:
      if weno_type == 'M':
        # This function mapps the weights according to Henick et. al.
        func = getattr(pyweno.cweno,'weights_m_' + points + '%03d%03d' % (k, n))
      elif weno_type == 'Z':
        # This function mapps the weights according to Henick et. al.
        func = getattr(pyweno.cweno,'weights_z_' + points + '%03d%03d' % (k, n))
      else:
        # The regular Shu weights
        func = getattr(pyweno.cweno,'weights_' + points + '%03d%03d' % (k, n))
    except AttributeError:
      raise ValueError('unsupported non-linear weights: '
                       + 'k=%d not generated' % (2*k-1))
    func(smoothness, weights)

  assert(weights.shape == (N,n,k))

  # reconstruct
  qr = np.zeros((N,n))

  if points == '+/-':
    # XXX: finish this off
    raise NotImplementedError, '+/- not implemented yet'

  try:
    func = getattr(pyweno.cweno,
                   'reconstruct_' + points + '%03d%03d' % (k, n))
  except AttributeError:
    raise ValueError('unsupported WENO reconstruction: '
                     + 'points=%s, k=%d, n=%d not generated'
                        % (points, 2*k-1, n))

  func(q, weights, qr)

  # done!
  if squeeze:
    qr = qr.squeeze()
    if weights is not None and return_weights:
      weights = weights.squeeze()

  if return_smoothness and return_weights:
    return (qr, smoothness, weights)

  if return_smoothness:
    return (qr, smoothness)

  if return_weights:
    return (qr, weights)

  return qr
