f_getPaxPropsFromVolumes<-function(M){
    Minter<-M
        
    vlambda_prop<-diag(Minter)
    vlambda_prop<-vlambda_prop/sum(vlambda_prop)
    diag(Minter)<-0
    


    Minter_sym<-Minter+t(Minter)
    Minter_sym[lower.tri(Minter_sym,diag=TRUE)]<-NA
    Minter_sym<-Minter_sym/sum(Minter_sym,na.rm=TRUE)
    
    Minter<-Minter/sum(Minter,na.rm=TRUE)
    
    list(Minter_sym=Minter_sym,Minter=Minter,vlambda_prop=vlambda_prop)
}




