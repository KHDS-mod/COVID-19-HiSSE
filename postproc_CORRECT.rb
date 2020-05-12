resdir                 = "RES_HISSE_ST6NOULTRA_CORRECT_1"
burnin                 = 5
n_time_slices          = 100

tr                    <- readTrees("ALL_CLEANED_DAT_ST6NOULTRA/tru.newick")[1]
ast                    = readAncestralStateTrace(resdir+"/stoch_char_map.log")
ancestralStateTree(tr, ast, resdir+"/ast.tree", burnin=0.1, reconstruction="marginal")
summarizeCharacterMaps(ast, tr, file=resdir+"/ev.tsv", burnin=0.1)   # What?
char_map_tree = characterMapTree(tree = tr, 
                 ancestral_state_trace_vector = ast, 
                 character_file=resdir+"/marginal_character.tree", 
                 posterior_file=resdir+"/marginal_posterior.tree", 
                 burnin=burnin, 
                 num_time_slices=n_time_slices)
