## No constraints
nop_ct = function (likobj) likobj

## All mu's are zero
nullmu_ct = function (likobj, nstates) {
  ifree = seq_len(nstates * (nstates + 1))[-((nstates+1):(2*nstates))]
  constrain.i(likobj, numeric(nstates*(nstates+1)), ifree)
}
nullmu_from_euclid = function (par, nstates) c(par[1:nstates], numeric(nstates), par[(nstates+1):(length(par))])
nullmu_to_euclid = function (par, nstates) par[-((nstates+1):(2*nstates))]


