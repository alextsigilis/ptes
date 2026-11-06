#
# hosa.py
#
# This file contains functions for Higher Order Spectral Analysis (HOSA)
#
import jax
import numpy as np
import jax.numpy as jnp

# This function estimates the bispectrum on a single segment
@jax.jit
def bisp(x):
    M = x.shape[0]
    y = (x - x.mean()) * jnp.hanning(M)
    Y = jnp.fft.fft(y)
    Y12 = jnp.vstack([jnp.roll(Y,-i) for i in range(M)])
    return jnp.einsum("i,j->ij", Y, Y) * jnp.conj(Y12)

# Map over segments
bisp_seg = jax.vmap(bisp, 0, 0)

# Direct Implementation of the Bispectrum
@jax.jit
def bispd(x):
    B = bisp_seg(x)
    return jnp.mean(B, axis=0)

# Map over channels
bispectrum = jax.jit(jax.vmap(bispd, 0, 0))

@jax.jit
def bic(x):
    B = bisp_seg(x)
    return jnp.abs(jnp.sum(B, axis=0)) / jnp.sqrt(jnp.sum(jnp.abs(B**2)))

# Map over all channels
bicoher = jax.jit(jax.vmap(bic, in_axes=0, out_axes=0))