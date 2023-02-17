#####################################################
##  This code produces a graph depicting the orthoganol regression
##    between two species. Orthogonal regression is an approach that minimizes
##    the perpendicular distrances between the data points (unlike traditional 
##    linear regressions that only minimize the sum of squared vertical 
##    distances on the y-axis). 

##  Regressions are performed on geometric mean values for each chemical. This
##    code also plots the raw experimental values, allowing the reader to see the
##    full experimental variability for each chemical.

##  This code recreates Figure 3 in Connors et al. 2022. ETC. 41(1)134-147. 

##  Code requires a data.table with columns named "CAS", "Latin.name", and "Effect.value"



library(data.table)
##Sample data produces Figure 3.b in in Connors et al. 2022. ETC. 41(1)134-147. 
#### import file 

d1.c <- fread("rat.csv") #update to appropriate file path



logTicks <- function(side.plot,tck.length=-0.01,tck.lwd=1.5){
  if (side.plot==1) logRange<-par("usr")[1:2]
  if (side.plot==2) logRange<-par("usr")[3:4]
  bottom<-ceiling(logRange[1])
  top<-floor(logRange[2])
  tickPositions<-unlist(sapply(bottom:top,FUN=function(x)x+log10(1:9)))
  tickPositions <-  tickPositions[ tickPositions>logRange[1]  & tickPositions<logRange[2] ]
  axis(side=side.plot,at=10^tickPositions,label=F,tck=tck.length,lwd=tck.lwd)
  tickPositions<-bottom:top
  tickPositions <-  tickPositions[ tickPositions>logRange[1]  & tickPositions<logRange[2] ]
  axis(side=side.plot,at=10^tickPositions,label=F,tck=2*tck.length,lwd=tck.lwd)
}

plotsForPaper <- function(
  chemList,
  xSpecies,
  xData,
  xVar,
  ySpecies,
  yData,
  yVar,
  mainTitle,
  pageText,
  pointsCEX=1.3,
  linesLWD=3,
  grayLWD=2,
  regCOL='blue',
  refCOL='gray',
  meansOnly=FALSE,
  speciesList=SpeciesList,
  noint=FALSE,
  axisLimits=10^(c(-3,8)+c(1,-1)*.3), ##7.28.2021  trimmed from -3,9 to -3,8
  confidenceAlpha=5){
  
  my.pcurve <- function(indata,myRange=c(-5,5)){
    ###indata is assumed to contain columns y and x, with the natural interpretation
    
    ###loadings of the pca give the slope of the orthog. regression line
    ###the line goes through the center of the data
    
    print(summary(pcObject <- princomp(indata)))
    myLoadings <- unclass(loadings(pcObject))[,1]
    OR.slope <- myLoadings["y"]/myLoadings["x"]
    #oreg line goes through the center of the data, so...
    OR.int <- mean(indata$y)-OR.slope*mean(indata$x)
    ###Rsq is empirical in this setting.  The first principal component will always explain at least half
    ###of the total variation, since the first two components will explain 100% of it in this simple setting,
    ###and the first must explain at least as much as the second.
    princomp.rsq <- ((summary(pcObject)$sdev[1]^2)/sum(summary(pcObject)$sdev^2)-.5)/.5
    #princomp.rsq <- (summary(pcObject)$varmat[2,1]-.5)/.5   ### the splus version
    print(c(princomp.rsq=princomp.rsq,corr.rsq=(corrCoef<-cor(indata$x,indata$y))^2))
    charRoots <- apply(pcObject$scores,2,var)
    sphereF <- (nrow(indata)-2)*diff(charRoots)^2/(8*prod(charRoots))
    sphereP <- 1 - pf(sphereF,2,nrow(indata)-2)
    if (sphereP<.05){
      delt <- asin(sqrt(prod(charRoots)/(nrow(indata)-2))*2*qt(1-(confidenceAlpha/100)/2,nrow(indata)-2)/
                     (-diff(charRoots)))/2
      #limits:
      clims<-c(tan(atan(OR.slope)+delt),tan(atan(OR.slope)-delt))
    }
    if (sphereP>=.05) clims <- c(NA,NA)
    print(c(int=OR.int,slope=OR.slope,rsq=princomp.rsq,sphereF=sphereF,sphereP=sphereP,lims=clims))
    c(int=OR.int,slope=OR.slope,rsq=princomp.rsq,sphereF=sphereF,sphereP=sphereP,lims=clims,corr=corrCoef,CI.level=100-confidenceAlpha)
  }
  
  noIntOReg <- function(slopes,y,x,plotTF=F){
    sapply(slopes,FUN=function(slope,plotsTF,x,y){
      #orthogonal slope is -1/slope
      oSlope <- -1/slope
      #orthogonal line:
      oInt <- y - oSlope*x
      #lines meet where y1=a1+b1*x1, and y2=a2+b2*x2 equal each other:
      oProjectionX <- -oInt/(oSlope-slope)
      oProjectionY <- slope*oProjectionX
      if (plotsTF){
        par(pty='s')
        fullRange <- c(0,max(c(x,y)))
        plot(x=c(x,oProjectionX),y=c(y,oProjectionY),xlim=fullRange,ylim=fullRange,type='n')
        points(x=x,y=y)
        segments(x1=x,x2=oProjectionX,y1=y,y2=oProjectionY)
        abline(a=0,b=slope)}
      sum( ((x-oProjectionX)^2 + (y-oProjectionY)^2) )
    },plotsTF=plotTF,x=x,y=y)
  }
  
  if (is.null(axisLimits)) axisRange <- range(na.omit(c(xData[[xVar]],yData[[yVar]])))
  if (!is.null(axisLimits)) axisRange <- axisLimits
  axisLine <- 3.0
  plot(x=1000,y=1000,xlim=axisRange,ylim=axisRange,log='xy',type='n',xlab='',ylab='',axes=F,bty="n")
  abline(a=0,b=1,col=refCOL,lty=1)
  
  mtext(side=1,adj=.5,text=xSpecies,line=axisLine,cex=1.2)
  mtext(side=2,adj=.5,text=ySpecies,line=axisLine-0.5,cex=1.2)
  ###need this at all?
  xunits <- unitAdj(10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]),by=1))
  axis(side=1,
       at=10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]),by=1),
       label=do.call(
         "expression",
         lapply(seq(ceiling(par("usr")[1]),floor(par("usr")[2]),by=1),FUN=function(x)bquote(10^.(x)))),
       #			label=xunits[[1]],
       las=2,adj=1,cex.axis=1)
  
  yunits <- unitAdj(10^seq(ceiling(par("usr")[3]),floor(par("usr")[4]),by=1))
  axis(side=2,
       at=10^seq(ceiling(par("usr")[3]),floor(par("usr")[4]),by=1),
       #			label=yunits[[1]],
       label=do.call(
         "expression",
         lapply(seq(ceiling(par("usr")[1]),floor(par("usr")[2]),by=1),FUN=function(x)bquote(10^.(x)))),
       las=2,adj=1,cex.axis=1)
  #plots are all square, so do units on x- y-axes simultaneously
  redCol <- rgb(213, 62, 79,maxColorValue=256)
  greenCol <- rgb(77, 146, 33,maxColorValue=256)
  yellowCol <- rgb(255,255,150,maxColorValue=256)
  boxHeight <- 0.3
  boxVals <- 10^(par("usr")[3]+c(.0,boxHeight,boxHeight,.0))
  polygon(x=c(10^-3,10^-3,10^0,10^0),y=boxVals,density=-1,col="gray")
  text(x=10^-1.5,y=10^(par("usr")[3]+boxHeight/2),labels=expression(plain(ng)/plain(L)))
  
  polygon(x=c(10^0,10^0,10^3,10^3),y=boxVals,density=-1,col=redCol)
  text(x=10^1.5,y=10^(par("usr")[3]+boxHeight/2),labels=expression(mu*plain(g)/plain(L)))
  polygon(x=c(10^3,10^3,10^6,10^6),y=boxVals,density=-1,col=yellowCol)
  text(x=10^4.5,y=10^(par("usr")[3]+boxHeight/2),labels=expression(plain(mg)/plain(L)))
  polygon(x=c(10^6,10^6,10^9,10^9),y=boxVals,density=-1,col=greenCol)
  text(x=10^7.5,y=10^(par("usr")[3]+boxHeight/2),labels=expression(plain(g)/plain(L)))
  
  polygon(y=c(10^-3,10^-3,10^0,10^0),x=boxVals,density=-1,col="gray")
  text(y=10^-1.5,x=10^(par("usr")[3]+boxHeight/2),labels=expression(plain(ng)/plain(L)),srt=-90)
  
  polygon(y=c(10^0,10^0,10^3,10^3),x=boxVals,density=-1,col=redCol)
  text(y=10^1.5,x=10^(par("usr")[3]+boxHeight/2),labels=expression(mu*plain(g)/plain(L)),srt=-90)
  polygon(y=c(10^3,10^3,10^6,10^6),x=boxVals,density=-1,col=yellowCol)
  text(y=10^4.5,x=10^(par("usr")[3]+boxHeight/2),labels=expression(plain(mg)/plain(L)),srt=-90)
  polygon(y=c(10^6,10^6,10^9,10^9),x=boxVals,density=-1,col=greenCol)
  text(y=10^7.5,x=10^(par("usr")[3]+boxHeight/2),labels=expression(plain(g)/plain(L)),srt=-90)
  if(FALSE){for(iu in unique(xunits[[2]])){
    values <- seq(ceiling(par("usr")[1]),floor(par("usr")[2]),by=1)[xunits[[2]]==iu]
    if (substring(iu,1,1)!='g') values<-c(values,max(values)+log10(8))
    midValue <- mean(range(values))
    labelText <- iu
    if (tolower(labelText)=="ug/l") labelText <- expression(mu*plain(g)/plain(L))
    text(x=10^midValue,y=.027,label=labelText,adj=c(.5,0))
    text(y=10^midValue,x=.035,label=labelText,adj=c(.5,.5),srt=90)
  }}
  
  logTicks(side=1)
  logTicks(side=2)
  #box(lwd=1.5)
  i <- unique(yData[,1])[1]
  xmeans<-ymeans<-xns<-yns<-seq(along=chemList)
  print(rbind(dim(yData),dim(xData)))
  
  for (i in seq(along=chemList)){
    xdata<-xData[[xVar]][as.character(xData[,"CAS"])==as.character(chemList[i])]
    ydata<-yData[[yVar]][as.character(yData[,"CAS"])==as.character(chemList[i])]
    print(str(xData))
    print(str(yData))
    xData<<-xData
    yData<<-yData
    print(chemList[i])
    print(c(x=xdata))
    print(c(y=ydata))
    cat("\n\n")
    xdata[xdata==0] <- NA
    ydata[ydata==0] <- NA
    xdata <- na.omit(xdata)
    ydata <- na.omit(ydata)
    xmean<-10^mean(log10(xdata))
    ymean<-10^mean(log10(ydata))
    xmeans[i]<-xmean
    ymeans[i]<-ymean
    xns[i]<-sum(!is.na(xdata))
    yns[i]<-sum(!is.na(ydata))
    if (!meansOnly){#identify the points
      rgbVal <- 0
      alphaVal <- .2
      points(y=rep(ymean,length(xdata)),x=xdata,cex=grayCEX,col=rgb(rgbVal,rgbVal,rgbVal,alphaVal))
      points(y=ydata,x=rep(xmean,length(ydata)),cex=grayCEX,col=rgb(rgbVal,rgbVal,rgbVal,alphaVal))
      lines(y=rep(ymean,length(xdata)),x=sort(xdata),col=rgb(rgbVal,rgbVal,rgbVal,alphaVal),lwd=grayLWD)
      lines(y=sort(ydata),x=rep(xmean,length(ydata)),col=rgb(rgbVal,rgbVal,rgbVal,alphaVal),lwd=grayLWD)
    }
  }
  print(data.frame(xmeans,ymeans,chemList))
  meanData<-na.omit(data.frame(y=ymeans,x=xmeans,chem=chemList,xsize=xns,ysize=yns,chemNo=1:length(chemList)))
  print(meanData)
  # stop()
  ypos<-rev(seq(par("usr")[1],par("usr")[2]-.1,length=max(45,nrow(meanData))))
  
  if (nrow(meanData)>0){
    if (nrow(meanData)>4){
      logMeanData<-data.frame(x=log10(meanData$x),y=log10(meanData$y),chem=meanData$chem)
      print(na.omit(logMeanData))
      if (!noint){
        lmObject<-lm(y~x,data=na.omit(logMeanData))
        xline<-c(min(c(logMeanData$x,logMeanData$y)),max(c(logMeanData$x,logMeanData$y)))+c(-.25,+.25)
        yline<-predict(lmObject,new=data.frame(x=xline))
        ypos2<-(seq(par("usr")[1],par("usr")[2]-.1,length=25))
        pcurveCoef <-my.pcurve(logMeanData[,c("x","y")])
        print(pcurveCoef)
        newY<-pcurveCoef[1]+pcurveCoef[2]*xline
        lines(x=10^xline,y=10^newY,lwd=linesLWD,col=regCOL)
      }
      if (noint){
        lmObject<-lm(y~-1+x,data=logMeanData)
        xline<-c(min(c(logMeanData$x,logMeanData$y)),max(c(logMeanData$x,logMeanData$y)))+c(-.25,+.25)
        yline<-coef(lmObject)*xline
        lines(x=10^xline,y=10^yline,col='lightblue')
        ypos2<-(seq(par("usr")[1],par("usr")[2]-.1,length=25))
        mtext(side=1,at=10^(par("usr")[2]+3.8),line=2,text=paste("slope=",format(coef(lmObject),digits=3),sep=''),adj=1,col="lightblue")
        testSlopes <- unique(sort(c(seq(1,2,by=.0001),1/seq(1,2,by=.0001))))
        oregSSEs <- noIntOReg(x=logMeanData$x,y=logMeanData$y,slopes=testSlopes)
        oregSlope <- testSlopes[oregSSEs==min(oregSSEs)]
        print(oregSlope)
        yline<-oregSlope*xline
        lines(x=10^xline,y=10^yline,col="blue")
        mtext(side=1,at=10^(par("usr")[2]+7.0),line=2,text=paste("slope=",format(oregSlope,digits=3),sep=''),adj=1,col="blue")
      }
    }
  }
  points(y=ymeans,x=xmeans,col=regCOL,pch=16,cex=pointsCEX)
  points(y=ymeans,x=xmeans,cex=pointsCEX,col=pointsCOL)
  invisible(list(data=logMeanData,results=pcurveCoef))
}

addRegInfo <- function(foo,corrTF,yLC.TF=T,yStr="D. magna",xEC.TF=T,xStr="C. dubia",eqnCOL=regCOL){
  print(foo[[2]])
  mySlope <- format(c(round(foo[[2]]["slope.y"],3),.001))[1]
  myInt <- format(c(round(abs(foo[[2]]["int.y"]),3),.001))[1]
  signInt <- sign(foo[[2]]["int.y"])
  topLine <- -4.10
  xLocation <- 100
  myLowerCI <- format(round(min(foo[[2]][regexpr("lims",names(foo[[2]]))>0]),2),nsmall=2)
  myUpperCI <- format(round(max(foo[[2]][regexpr("lims",names(foo[[2]]))>0]),2),nsmall=2)
  myConfLevel <- format(foo[[2]][regexpr("CI.level",names(foo[[2]]))>0])
  if (signInt==-1)
    mtext(side=1,line=topLine,at=xLocation,adj=0,col=eqnCOL,
          text=bquote(
            paste(
              .(yStr),scriptscriptstyle(" "),
              .(ifelse(yLC.TF,"LC50","LC50"))==.(mySlope),
              scriptscriptstyle("")%*%scriptscriptstyle(""),
              .(xStr),scriptscriptstyle(" "),
              .(ifelse(xEC.TF,"LC50","LC50"))-.(myInt))))
  if (signInt==+1)
    mtext(side=1,line=topLine,at=xLocation,adj=0,col=eqnCOL,
          text=bquote(
            paste(
              .(yStr),scriptscriptstyle(" "),
              .(ifelse(yLC.TF,"LC50","LC50"))==.(mySlope),
              scriptscriptstyle("")%*%scriptscriptstyle(""),
              .(xStr),scriptscriptstyle(" "),
              .(ifelse(xEC.TF,"LC50","LC50"))+.(myInt))))
  if (!corrTF) mtext(side=1,line=topLine+2.2,at=xLocation,adj=0,
                     text=bquote(paste(R^2==.(format(round(foo[[2]]["rsq.Comp.1"],2),digits=2)),", ",.(nrow(foo[[1]]))," chemicals")))
  if (corrTF)  mtext(side=1,line=topLine+2.2,at=xLocation,adj=0,
                     text=bquote(paste(r==.(format(round(foo[[2]]["corr"],2),digits=2)),", ",.(nrow(foo[[1]]))," chemicals")))
  mtext(side=1,line=topLine+1.1,at=xLocation,adj=0,
        text=bquote(paste("Slope ",.(myConfLevel),"% CI: ",
                          group("(",
                                list(
                                  .(myLowerCI),
                                  .(myUpperCI)),")"))))
}

unitAdj <- function(values){
  units <- rep("ug/L",length(values))
  
  resets <- (values>=1000)
  values[resets] <- values[resets]/1000
  units[resets] <- "mg/L"
  
  resets <- (values>=1000)
  values[resets] <- values[resets]/1000
  units[resets] <- "g/L"
  
  list(sapply(values,format,scientific=F),units)
}

##########################################
## Following code executes the graph. Update as needed based on your base file
##########################################


# pdf(file="Paper Figures Test.pdf",height=8,width=8)
doTitles<-FALSE
#pdf(file="Figures/MainFigures99.9.pdf",height=8,width=8,paper="USr")
par2start <- par()
par(pty='s',omi=c(0,0,0,0),mai=c(.8,.8,.1,.1))
SpeciesList<-c("Rodent","Caenorhabditis elegans")
ChemList<-unique(c(d1.c[Latin.name=="Rodent", CAS], d1.c[Latin.name=="Caenorhabditis elegans", CAS]))
isJava<-(substring(names(dev.cur()),1,4)=="java")
pageno <- 0

lineColor1 <- lineColor2 <- rgb(252, 141, 89,maxColorValue=256)
refColor <- "cyan"
pointsCOL <- "red"
grayCEX <- 0.8

dm <- as.data.table(d1.c[Latin.name=="Rodent", list(CAS, Effect.value)])
##  plotting algorithm expects data in units ug/L. My data is in mg/L, adjusting units manually.
dm[, Effect.value :=  Effect.value * 1000]
dm <- as.data.frame(dm)


cd <- as.data.table(d1.c[Latin.name=="Caenorhabditis elegans", list(CAS, Effect.value)])
##  plotting algorithm expects data in units ug/L. My data is in mg/L, adjusting units manually.
cd[, Effect.value :=Effect.value * 1000]
cd <- as.data.frame(cd)

plotOutput <- plotsForPaper(
  ySpecies=expression(paste(italic("Rodent"), " toxicity", phantom(0),LD50,phantom(0), "(mg/kg)")),
  xSpecies=expression(paste(italic("Caenorhabditis elegans"), " toxicity",phantom(0),EC10,phantom(0),"(",mu,"g/L)")),
  chemList = ChemList,
  yData = dm,
  yVar = "Effect.value",
  xData = cd,
  xVar = "Effect.value",
  mainTitle = "",
  pageText = paste("Em vs El, pg.",pageno),
  meansOnly=FALSE,
  refCOL="#80FFFFE6",
  regCOL=lineColor1,
  confidenceAlpha=5
)

# addRegInfo(plotOutput,corr=T,eqnCOL=lineColor1)
dim(plotOutput[[1]])
if(doTitles)mtext(side=3,line=-1,adj=0,text="D.magna and C.dubia Chronic Toxicity")
# dev.off()


