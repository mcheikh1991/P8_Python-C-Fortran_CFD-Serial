#include <Python.h>
#include <numpy/ndarrayobject.h>

void
weights_z_left005001 (const double *restrict sigma, int n, int ssi, int ssr,
		    double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double acc, sigma0, sigma1, sigma2, sigma3, sigma4, omega1, omega3, omega0, omega2, omega4;
  double tau;
  for (i = 5; i < n - 5; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];

      tau = abs(sigma0 - sigma4);

      acc = 0.0;
      omega0 = (+0.00793650793650794)*( 1+ (tau /(sigma0 + 1.0e-6)) * (tau /(sigma0 + 1.0e-6)) );
      acc = acc + omega0;
      omega1 = (+0.158730158730159)*( 1+ (tau /(sigma1 + 1.0e-6)) * (tau /(sigma1 + 1.0e-6)) );
      acc = acc + omega1;
      omega2 = (+0.476190476190476)*( 1+ (tau /(sigma2 + 1.0e-6)) * (tau /(sigma2 + 1.0e-6)) );
      acc = acc + omega2;
      omega3 = (+0.317460317460317)*( 1+ (tau /(sigma3 + 1.0e-6)) * (tau /(sigma3 + 1.0e-6)) );
      acc = acc + omega3;
      omega4 = (+0.0396825396825397)*( 1+ (tau /(sigma4 + 1.0e-6)) * (tau /(sigma4 + 1.0e-6)) );
      acc = acc + omega4;
      
      omega0 = (omega0) / (acc);
      omega1 = (omega1) / (acc);
      omega2 = (omega2) / (acc);
      omega3 = (omega3) / (acc);
      omega4 = (omega4) / (acc);
      omega[i * wsi + 0 * wsl + 0 * wsr + 0] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr + 0] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr + 0] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr + 0] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr + 0] = omega4;
    }
}

PyObject *
py_weights_z_left005001 (PyObject * self, PyObject * args)
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

  weights_z_left005001 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}