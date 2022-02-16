source("getPropsFromVolumes.R")


sabreData<-read.csv("sabre.csv",header=TRUE,sep=",")

v_allregions<-unique(c(unique(sabreData$Origin.Region.Name),unique(sabreData$Destination.Region.Name)))


v_regions_remove<-c("ANTARCTICA","USA TERRITORIAL POSSESSION","","UNASSIGNED")
l_regions<-list(Afr=c("SOUTHERN AFRICA","WEST AFRICA","NORTH AFRICA","EAST AFRICA","CENTRAL AFRICA" ),
Asia=c("FAR EAST ASIA","ASIA SUB-CONTINENT","SOUTHEAST ASIA","GULF","MIDDLE EAST","CENTRAL ASIA","Far East Asia" ),
Eur=c("EASTERN EUROPE","WESTERN EUROPE"),
NAm=c("NORTH AMERICA","CARIBBEAN","CENTRAL AMERICA"),
Oc=c("AUSTRALIA","PACIFIC"),
SAm=c("SOUTH AMERICA"))

v_rowremove<-unique(c(which(sabreData$Origin.Region.Name %in% v_regions_remove),which(sabreData$Destination.Region.Name %in% v_regions_remove)))
sabreData<-sabreData[-v_rowremove,,drop=FALSE]

sabreData<-rbind(sabreData,Origin.Continent.Name=NA,Destination.Continent.Name=NA)

for (i in 1:length(l_regions)){
    sabreData[which(sabreData$Origin.Region.Name %in% l_regions[[i]]),"Origin.Continent.Name"]<-names(l_regions)[i]
    sabreData[which(sabreData$Destination.Region.Name %in% l_regions[[i]]),"Destination.Continent.Name"]<-names(l_regions)[i]
}

## Georgia is classified as Asia in NextStrain
sabreData[which(sabreData$Origin.Country.Name == "GEORGIA"),"Origin.Continent.Name"]<- "Asia"
sabreData[which(sabreData$Destination.Country.Name  == "GEORGIA"),"Destination.Continent.Name"]<- "Asia"

l_PaxVol<-vector("list",length(unique(sabreData$year)))
names(l_PaxVol)<-unique(sabreData$year)

for (i in 1:length(l_PaxVol)){
    year<-names(l_PaxVol)[i]
    l_PaxVol[[i]]<-vector("list",length(unique(sabreData$month[which(sabreData$year==year)]))+1)
    names(l_PaxVol[[i]])<-c(paste0("year_",year), intersect(c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"),unique(sabreData$month[which(sabreData$year==year)])))

    l_PaxVol[[i]][[1]]<-matrix(NA,length(l_regions),length(l_regions))
    colnames(l_PaxVol[[i]][[1]])<-names(l_regions)
    rownames(l_PaxVol[[i]][[1]])<-names(l_regions)
    for (k1 in 1:nrow(l_PaxVol[[i]][[1]])){
	for (k2 in 1:ncol(l_PaxVol[[i]][[1]])){
	    l_PaxVol[[i]][[1]][k1,k2]<-sum(sabreData$Passengers[intersect(which(sabreData$year==year),intersect(which(sabreData$Origin.Continent.Name==names(l_regions)[k1]),which(sabreData$Destination.Continent.Name==names(l_regions)[k2])))])
	}
    }   
    for (j in 2:length(l_PaxVol[[i]])){
	month<-names(l_PaxVol[[i]])[j]
	l_PaxVol[[i]][[j]]<-matrix(NA,length(l_regions),length(l_regions))
	colnames(l_PaxVol[[i]][[j]])<-names(l_regions)
	rownames(l_PaxVol[[i]][[j]])<-names(l_regions)
	for (k1 in 1:nrow(l_PaxVol[[i]][[j]])){
	    for (k2 in 1:ncol(l_PaxVol[[i]][[j]])){
		l_PaxVol[[i]][[j]][k1,k2]<-sum(sabreData$Passengers[intersect(intersect(which(sabreData$year==year),which(sabreData$month==month)),intersect(which(sabreData$Origin.Continent.Name==names(l_regions)[k1]),which(sabreData$Destination.Continent.Name==names(l_regions)[k2])))])
	    }
	}
    }
}
l_PropPaxVol<-l_PaxVol

for(i in 1:length(l_PropPaxVol)){
    for(j in 1:length(l_PropPaxVol[[i]])){
	l_PropPaxVol[[i]][[j]]<-f_getPaxPropsFromVolumes(l_PaxVol[[i]][[j]])
    }
}

save(l_PropPaxVol,file="PropPaxVol.RData")
save(l_PaxVol,file="PaxVol.RData")
