printSeed()

NUM_TOTAL_SPECIES     = 3585
NUM_STATES            = 6
NUM_HIDDEN            = 2
NUM_RATES             = NUM_STATES * NUM_HIDDEN
H                     = 0.59

tr                    <- readTrees("ALL_CLEANED_DAT_ST6NOULTRA/tru.newick")[1]
st                    <- readCharacterDataDelimited("ALL_CLEANED_DAT_ST6NOULTRA/tipstate.csv",
                                                    type="NaturalNumbers",
                                                    stateLabels=NUM_STATES,
                                                    delimiter=",",
                                                    headers=TRUE)
st_exp                <- st.expandCharacters( NUM_HIDDEN )
taxa                  <- tr.taxa()
tree_length           <- tr.treeLength()

moves                 = VectorMoves()
monitors              = VectorMonitors()

sig0                  ~ dnExponential( 1.0 / H )
moves.append( mvScale(sig0, lambda=1, tune=true, weight=2.0) )
sig0.setValue(3.678277)  # Init. from previous burn-in run, iter=949
# This burn-in state has posterior = -1197.072
lambda_hid_unorm := fnDiscretizeDistribution(dnLognormal(ln(1.0), sig0), NUM_HIDDEN)
lambda_hid := lambda_hid_unorm / mean(lambda_hid_unorm)

pax_intra <- [0.019, 0.692, 0.13, 0.118, 0.015, 0.026]
c1        <- 6.0
muglob ~ dnLoguniform( 1E-6, 1E2)
muglob.setValue(1.064082e-06)     # Init. from previous burn-in run, iter=949
lambda_sum  ~ dnLoguniform( 1E-6, 1E2)
lambda_sum.setValue(8.022098)     # Init. from previous burn-in run, iter=949
moves.append( mvScale(lambda_sum, lambda=1.0, tune=true, weight=3.0) )
lambda_prop ~ dnDirichlet(c1*pax_intra)
lambda_prop.setValue([0.116017, 0.192851, 0.241205, 0.197265, 0.0763464, 0.176315])  # Init. from previous burn-in run, iter=949
moves.append( mvDirichletSimplex(lambda_prop, tune=true, weight=3.0) )
for (i in 1:NUM_STATES) {
    lambda_obs[i] := lambda_sum * lambda_prop[i]
}

for (j in 1:NUM_HIDDEN) {
    for (i in 1:NUM_STATES) {
	idx = i+(j*NUM_STATES)-NUM_STATES
	lambda[idx] := lambda_obs[i] * lambda_hid[j]
        mu[idx] := muglob
    }
}
moves.append( mvScale(muglob, lambda=1.0, tune=true, weight=3.0) )

qpr := tr.treeLength() / 200
pax_inter <- [0.052,
              0.139, 0.297,
              0.012, 0.086, 0.244,
              0.001, 0.048, 0.014, 0.008,
              0.002, 0.004, 0.044, 0.049, 0.001]
q_obs_sum ~ dnGamma(15.0, qpr)
q_obs_sum.setValue(1.128907)   ## Init. from previous burn-in run, iter=949
moves.append( mvScale(q_obs_sum, lambda=1.0, tune=true, weight=4.0) )
c1q       <- 1.0
q_obs_prop ~ dnDirichlet(c1q*pax_inter)
q_obs_prop.setValue([0.000122084,0.0478,      0.100522,     0.0112694,     0.0941335,      0.259223,
                     0.0957059,     0.0896066,     0.0724381,      0.0909317,      0.0478552,
                     0.00254673,      0.0571963,     0.00935638,      0.0212932]) ## Init. from previous burn-in run, iter=949
moves.append( mvBetaSimplex(q_obs_prop, alpha=10.0, tune=true, weight=15.0) )
moves.append( mvDirichletSimplex(q_obs_prop, alpha=10.0, numCats=1, offset=1, tune=true, weight=15.0) )

## Now I want to loop through pairs of stuff
idx2 = 1
for (i in 2:NUM_STATES) {
    for (j in 1:(i-1)) {
	idx1 = (i-1)*(NUM_STATES-1)+j
        q_obs[idx1]                     := q_obs_sum * q_obs_prop[idx2]
        idx2 = idx2 + 1 
	q_obs[(j-1)*(NUM_STATES-1)+i-1] := q_obs[idx1]
    }
}
  
h ~ dnExp(qpr)
h.setValue(2.331006)  ## Init. from previous burn-in run, iter=949
moves.append( mvScale(h,lambda=0.5,tune=true,weight=5) )
for (i in 1:(NUM_HIDDEN * (NUM_HIDDEN - 1))) q_hid[i] := h

q := fnHiddenStateRateMatrix(q_obs, q_hid, rescaled=false)
pi = rep(0.05/(NUM_RATES-1), NUM_RATES)
pi[2] = 0.95/2 
pi[8] = 0.95/2
timetree ~ dnCDBDP( rootAge           = tr.rootAge(),
                    speciationRates   = lambda,
                    extinctionRates   = mu, 
                    Q                 = q,
                    delta             = 1.0,
                    pi                = pi,
                    rho               = 1.0,
                    condition         = "survival" )
timetree.clamp(tr)
timetree.clampCharData(st_exp)

mymodel = model(q)

### set up the monitors that will output parameter values to file and screen 
monitors.append( mnModel(filename=outdir+"/model"+chnname+".log", printgen=1) )
# monitors.append( mnJointConditionalAncestralState(tree=timetree,
# 						   cdbdp=timetree,
# 						   type="NaturalNumbers",
# 						   printgen=1,
# 						   withTips=true,
# 						   withStartStates=false,
# 						   filename="RES_HISSE_1/ancst.log") )
monitors.append( mnStochasticCharacterMap(cdbdp=timetree,
					  printgen=1,
					  filename=outdir+"/stoch_char_map"+chnname+".log",
					  include_simmap=true) )
monitors.append( mnScreen(printgen=1, lambda_obs, lambda_hid, q_obs, q_hid, muglob) )
mymcmc = mcmc(mymodel, monitors, moves, nruns=1, moveschedule="random", combine="mixed")
files = listFiles(outdir)
for (i in 1:files.size()) {
        if (files[i].find("chkpoint.log") != 0) {
                mymcmc.initializeFromCheckpoint(outdir+"/chkpoint"+chnname+".log")
                break
        }

}
mymcmc.run(generations=4500,
           checkpointFile=outdir+"/chkpoint"+chnname+".log",
           checkpointInterval=100,
           underPrior=!TRUE)
q()


#burnin=25
#n_time_slices = 500
#stoch = readAncestralStateTrace("RES_HISSE_1/stoch_char_map.log")
#summarizeCharacterMaps(stoch, tr, file="RES_HISSE_1/ev.tsv", burnin=0.1)
#char_map_tree = characterMapTree(tree=tr, 
#                 ancestral_state_trace_vector=stoch,
#                 character_file="RES_HISSE_1/margchar.tree", 
#                 posterior_file="RES_HISSE_1/margpstr.tree", 
#                 burnin=burnin, 
#                 num_time_slices=n_time_slices)
#q()






