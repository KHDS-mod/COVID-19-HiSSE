## What columns to select for each models?

contdict = c("Africa", "Asia", "Europe", "N.Amer", "Oceania", "S.Amer")
hidden   = c("A", "B")
qpair_select = c(1:5,7:10,13:15,19:20,25)
names(qpair_select)= c("Africa ↔ Asia","Africa ↔ Europe","Africa ↔ N.Amer","Africa ↔ Oceania","Africa ↔ S.Amer",
                       "Asia ↔ Europe","Asia ↔ N.Amer","Asia ↔ Oceania","Asia ↔ S.Amer",
                       "Europe ↔ N.Amer","Europe ↔ Oceania","Europe ↔ S.Amer",
                       "N.Amer ↔ Oceania","N.Amer ↔ S.Amer",
                       "Oceania ↔ S.Amer")
kdict = list(`full-withmu`         =39,
             `lambsame-withmu`     =34,
             `qpair-withmu`        =24,
             `lambsameqpair-withmu`=19,
             `qeq-withmu`         =10,
             `alleq-withmu`         =5)
coldict = list(
  `full-withmu`=c("Posterior",
                  "lambda_hid[1]","lambda_hid[2]",
                  sprintf("lambda_obs[%d]", 1:6),
                  "muglob",
                  "q_hid[1]",
                  sprintf("q_obs[%d]", 1:30),
                  "sig0"),
  `alleq-withmu`=c("Posterior",
                   "lambda_hid[1]","lambda_hid[2]",
                   "muglob",
                   "lambda_obs[1]",
                   "q_hid[1]",
                   "q_obs[1]",
                   "sig0"),
  `lambsame-withmu`=c("Posterior",
                      "lambda_hid[1]","lambda_hid[2]",
                      "muglob",
                      "lambda_obs[1]",
                      "q_hid[1]",
                      sprintf("q_obs[%d]", 1:30),
                      "sig0"),
  `lambsameqpair-withmu`=c("Posterior",
                           "lambda_hid[1]", "lambda_hid[2]",
                           "lambda_obs[1]",
                           "muglob",
                           "q_hid[1]",
                           sprintf("q_obs[%d]", qpair_select),
                           "sig0"),
  `qeq-withmu`=c("Posterior",
                 "lambda_hid[1]","lambda_hid[2]",
                 "muglob",
                  sprintf("lambda_obs[%d]", 1:6),
                 "q_hid[1]",
                 "q_obs[1]",
                 "sig0"),
  `qpair-withmu`=c("Posterior",
                   "lambda_hid[1]", "lambda_hid[2]",
                   sprintf("lambda_obs[%d]", 1:6),
                   "muglob",
                   "q_hid[1]",
                   sprintf("q_obs[%d]", qpair_select),
                   "sig0")
)
matchind = function (colname) as.integer(gsub("^[^(0-9)]*|\\]$", "", colname))
qcontdirection = function (colname) {
  ind = matchind(colname)
  fromind = ceiling(ind/5)
  toind = (1:6)[-fromind][ind-(fromind-1)*5]
  list(contdict[fromind], contdict[toind])
}
humanise_colnames = function(colselect, colname, qpair=F) {
    if (startsWith(colname, "lambda_hid[")) {
      howmany = sum(startsWith(colselect, "lambda_hid["))
      if (howmany == 1)
        return( "λ: hidden state" )
      else
        return( sprintf("λ: hidden state %s", hidden[matchind(colname)]) )
    } else if (startsWith(colname, "lambda_obs[")) {
      howmany = sum(startsWith(colselect, "lambda_obs["))
      if (howmany == 1)
        return( "λ" )
      else
        return( sprintf("λ: %s", contdict[matchind(colname)]) )
    } else if (startsWith(colname, "muglob")) {
      return( "μ" )
    } else if (startsWith(colname, "q_hid[")) {
      return( "q: hidden state" )
    } else if (startsWith(colname, "q_obs[")) {
      if (!qpair) {
        howmany = sum(startsWith(colselect, "q_obs["))
        if (howmany == 1)
          return( "q" )
        else
          return( do.call(sprintf, c("q: %s → %s", qcontdirection(colname))) )
      } else {
        return( do.call(sprintf, c("q: %s↔%s", qcontdirection(colname))))
      }
    } else if (startsWith(colname, "sig0")) {
      return( "σ" )
    } else {
      return( colname )
    }
}
