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

lambda_hid_unorm := fnDiscretizeDistribution(dnLognormal(ln(1.0), sig0), NUM_HIDDEN)
lambda_hid := lambda_hid_unorm / mean(lambda_hid_unorm)

pax_intra <- [0.019, 0.692, 0.13, 0.118, 0.015, 0.026]
c1        <- 6.0
muglob ~ dnLoguniform( 1E-6, 1E2)
lambda_sum  ~ dnLoguniform( 1E-6, 1E2)
lambda_sum.setValue(7.0)
moves.append( mvScale(lambda_sum, lambda=1.0, tune=true, weight=3.0) )
lambda_prop ~ dnDirichlet(c1*pax_intra)
lambda_prop.setValue(pax_intra)
moves.append( mvDirichletSimplex(lambda_prop, tune=true, weight=3.0) )
for (i in 1:NUM_STATES) {
    lambda_obs[i] := lambda_sum * lambda_prop[i]
}

for (j in 1:NUM_HIDDEN) {
    for (i in 1:NUM_STATES) {
        idx         = i+(j*NUM_STATES)-NUM_STATES
        lambda[idx] := lambda_obs[i] * lambda_hid[j]
        mu[idx]     := muglob
    }
}
moves.append( mvScale(muglob, lambda=1.0, tune=true, weight=3.0) )

## I think maybe there are 200 transitions in total...
qpr := tr.treeLength() / 200
pax_inter <- [0.025, 0.071, 0.007, 0.001, 0.001,
              0.026, 0.156, 0.045, 0.023, 0.002,
              0.068, 0.141, 0.121, 0.006, 0.021,
              0.006, 0.041, 0.123, 0.003, 0.024,
              0.001, 0.025, 0.008, 0.004, 0.001,
              0.001, 0.002, 0.023, 0.023, 0.001]
q_obs_sum ~ dnGamma(30.0, qpr)
q_obs_sum.setValue(30.0 / qpr)
moves.append( mvScale(q_obs_sum, lambda=1.0, tune=true, weight=8.0) )
c1q       <- 1.0
q_obs_prop ~ dnDirichlet(c1q*pax_inter)
q_obs_prop.setValue(pax_inter)
moves.append( mvBetaSimplex(q_obs_prop, alpha=10.0, tune=true, weight=30.0) )
moves.append( mvDirichletSimplex(q_obs_prop, alpha=10.0, numCats=1, offset=1, tune=true, weight=30.0) )
for ( i in 1:(NUM_STATES*(NUM_STATES-1)) ) {
    q_obs[i] := q_obs_sum * q_obs_prop[i]
}

h ~ dnExp(qpr)
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
monitors.append( mnScreen(printgen=1, lambda_prop, lambda_sum, lambda_hid, q_obs, q_hid, muglob) )
mymcmc = mcmc(mymodel, monitors, moves, nruns=1, moveschedule="random", combine="mixed")
files = listFiles(outdir)
for (i in 1:files.size()) {
        if (files[i].find("chkpoint.log") != 0) {
                mymcmc.initializeFromCheckpoint(outdir+"/chkpoint"+chnname+".log")
                break
        }

}
mymcmc.run(generations=1470,
           checkpointFile=outdir+"/chkpoint"+chnname+".log",
           checkpointInterval=30,
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




