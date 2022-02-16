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
sig0.setValue(4.731447)   ## Init. from previous burn-in run, iter=1137.
## This state has posterior = -1195.508
moves.append( mvScale(sig0, lambda=1, tune=true, weight=2.0) )

lambda_hid_unorm := fnDiscretizeDistribution(dnLognormal(ln(1.0), sig0), NUM_HIDDEN)
lambda_hid := lambda_hid_unorm / mean(lambda_hid_unorm)

lambda_obscomm ~ dnLoguniform( 1E-6, 1E2)
lambda_obscomm.setValue( 1.758809 ) ## Init. from previous burn-in run, iter=1137.
moves.append( mvScale(lambda_obscomm, lambda=1.0, tune=true, weight=3.0) )

muglob ~ dnLoguniform( 1E-1, 1E2)
muglob.setValue(0.1004459)
for (i in 1:NUM_STATES) lambda_obs[i] := lambda_obscomm
for (j in 1:NUM_HIDDEN) {
    for (i in 1:NUM_STATES) {
	idx = i+(j*NUM_STATES)-NUM_STATES
	lambda[idx] := lambda_obs[i] * lambda_hid[j]
        mu[idx] := muglob
    }
}
moves.append( mvScale(muglob, lambda=1.0, tune=true, weight=3.0) )
qpr := tr.treeLength() / 200
pax_inter <- [0.025, 0.071, 0.007, 0.001, 0.001,
              0.026, 0.156, 0.045, 0.023, 0.002,
              0.068, 0.141, 0.121, 0.006, 0.021,
              0.006, 0.041, 0.123, 0.003, 0.024,
              0.001, 0.025, 0.008, 0.004, 0.001,
              0.001, 0.002, 0.023, 0.023, 0.001]

q_obs_sum ~ dnGamma(30.0, qpr)
q_obs_sum.setValue(4.549051)
moves.append( mvScale(q_obs_sum, lambda=1.0, tune=true, weight=5.0) )
c1q       <- 1.0
q_obs_prop ~ dnDirichlet(c1q*pax_inter)
q_obs_prop.setValue([1.45947e-08,      0.175956,     0.0403727,
                     0.00723752,    0.00898902,   0.000162897,      0.032593,      0.032136,
                     0.0274247,    1.09992e-06,      0.0198763,       0.022911,      0.0540304,
                     0.0246419,      0.0156819,    0.000693605,      0.0122605,      0.0670947,
                     0.0200161,     0.00199223,     0.00221505,       0.017913,      0.0319211,
                     0.121554,     0.00987732,    1.24303e-05,     0.00104503,       0.117556,
                     0.13121,     0.00262573])   ## Init. from previous burn-in run, iter=1137.
moves.append( mvBetaSimplex(q_obs_prop, alpha=10.0, tune=true, weight=30.0) )
moves.append( mvDirichletSimplex(q_obs_prop, alpha=10.0, numCats=1, offset=1, tune=true, weight=30.0) )
for ( i in 1:(NUM_STATES*(NUM_STATES-1)) ) {
    q_obs[i] := q_obs_sum * q_obs_prop[i]
}


h ~ dnExp(qpr)
h.setValue(2.365262)
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
monitors.append( mnModel(filename=outdir+"/model"+chnname+".log", append=TRUE, printgen=1) )
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
monitors.append( mnScreen(printgen=1, lambda_obscomm, lambda_hid, q_obs, q_hid, muglob) )

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






