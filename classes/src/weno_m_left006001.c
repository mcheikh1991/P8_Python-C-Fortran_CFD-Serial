#include <Python.h>
#include <numpy/ndarrayobject.h>

void
weights_m_left006001 (const double *restrict sigma, int n, int ssi, int ssr,
		    double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double acc, sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, omega1, omega0, omega5, omega4,
    omega3, omega2;
  double d0, d1, d2, d3, d4, d5, sum_g, g0, g1, g2, g3, g4, g5;
  for (i = 6; i < n - 6; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      acc = 0.0;
      omega0 = (+0.00216450216450216) / ((sigma0 + 1.0e-6) * (sigma0 + 1.0e-6));
      acc = acc + omega0;
      omega1 = (+0.0649350649350649) / ((sigma1 + 1.0e-6) * (sigma1 + 1.0e-6));
      acc = acc + omega1;
      omega2 = (+0.324675324675325) / ((sigma2 + 1.0e-6) * (sigma2 + 1.0e-6));
      acc = acc + omega2;
      omega3 = (+0.432900432900433) / ((sigma3 + 1.0e-6) * (sigma3 + 1.0e-6));
      acc = acc + omega3;
      omega4 = (+0.162337662337662) / ((sigma4 + 1.0e-6) * (sigma4 + 1.0e-6));
      acc = acc + omega4;
      omega5 = (+0.012987012987013) / ((sigma5 + 1.0e-6) * (sigma5 + 1.0e-6));
      acc = acc + omega5;
      omega0 = (omega0) / (acc);
      omega1 = (omega1) / (acc);
      omega2 = (omega2) / (acc);
      omega3 = (omega3) / (acc);
      omega4 = (omega4) / (acc);
      omega5 = (omega5) / (acc);
      // Mapping the weights using Henrick et.al. method

      // Optimal weights
      d0 = 0.00216450216450216;
      d1 = 0.0649350649350649;
      d2 = 0.324675324675325;
      d3 = 0.432900432900433;
      d4 = 0.162337662337662;
      d5 = 0.012987012987013;
      
      sum_g = 0.0;
      // Mapping Function
      g0 = (omega0*(d0 + d0*d0 - 3.0*d0*omega0 + omega0*omega0))/(d0*d0 + omega0*(1.0-2.0*d0));
      sum_g = sum_g + g0;

      g1 = (omega1*(d1 + d1*d1 - 3.0*d1*omega1 + omega1*omega1))/(d1*d1 + omega1*(1.0-2.0*d1));
      sum_g = sum_g + g1;

      g2 = (omega2*(d2 + d2*d2 - 3.0*d2*omega2 + omega2*omega2))/(d2*d2 + omega2*(1.0-2.0*d2));
      sum_g = sum_g + g2;

      g3 = (omega3*(d3 + d3*d3 - 3.0*d3*omega3 + omega3*omega3))/(d3*d3 + omega3*(1.0-2.0*d3));
      sum_g = sum_g + g3;

      g4 = (omega4*(d4 + d4*d4 - 3.0*d4*omega4 + omega4*omega4))/(d4*d4 + omega4*(1.0-2.0*d4));
      sum_g = sum_g + g4;

      g5 = (omega5*(d5 + d5*d5 - 3.0*d5*omega5 + omega5*omega5))/(d5*d5 + omega5*(1.0-2.0*d5));
      sum_g = sum_g + g5;

      omega0 = g0 / sum_g;
      omega1 = g1 / sum_g;
      omega2 = g2 / sum_g;
      omega3 = g3 / sum_g;
      omega4 = g4 / sum_g;
      omega5 = g5 / sum_g;

      omega[i * wsi + 0 * wsl + 0 * wsr + 0] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr + 0] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr + 0] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr + 0] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr + 0] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr + 0] = omega5;
    }
}

PyObject *
py_weights_m_left006001 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (sigma_py->nd != 2 || sigma_py->descr->type_num != PyArray_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (omega_py->descr->type_num != PyArray_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(omega_py->nd >= 2 && omega_py->nd <= 4))
    {
      PyErr_SetString (PyExc_ValueError, "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = sigma_py->strides[0] / sizeof (double);
  ssr = sigma_py->strides[1] / sizeof (double);

  wsi = omega_py->strides[0] / sizeof (double);
  if (omega_py->nd == 3)
    {
      wsl = omega_py->strides[1] / sizeof (double);
      wsr = omega_py->strides[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = omega_py->strides[1] / sizeof (double);
    }

  weights_m_left006001 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}
