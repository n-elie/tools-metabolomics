# lib.r

# This function retrieve a xset like object
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return (xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, 'xcmsSet'))
        if (is.null(xset@phenoData$sample_group))
            sampclass(xset) = "."
        else
            sampclass(xset) <- xset@phenoData$sample_group
        if (!is.null(xset@phenoData$sample_name))
            rownames(xset@phenoData) = xset@phenoData$sample_name
        return (xset)
    }
}

#@author G. Le Corguille
#The function create a pdf from the different png generated by diffreport
diffreport_png2pdf <- function(filebase) {
    dir.create("pdf")

    pdfEicOutput = paste0("pdf/",filebase,"-eic_pdf.pdf")
    pdfBoxOutput = paste0("pdf/",filebase,"-box_pdf.pdf")

    system(paste0("gm convert ",filebase,"_eic/*.png ",pdfEicOutput))
    system(paste0("gm convert ",filebase,"_box/*.png ",pdfBoxOutput))

}

#@author G. Le Corguille
#This function convert if it is required the Retention Time in minutes
RTSecondToMinute <- function(variableMetadata, convertRTMinute) {
    if (convertRTMinute){
        #converting the retention times (seconds) into minutes
        print("converting the retention times into minutes in the variableMetadata")
        variableMetadata[,"rt"]=variableMetadata[,"rt"]/60
        variableMetadata[,"rtmin"]=variableMetadata[,"rtmin"]/60
        variableMetadata[,"rtmax"]=variableMetadata[,"rtmax"]/60
    }
    return (variableMetadata)
}

#@author G. Le Corguille
#This function format ions identifiers
formatIonIdentifiers <- function(variableMetadata, numDigitsRT=0, numDigitsMZ=0) {
    splitDeco = strsplit(as.character(variableMetadata$name),"_")
    idsDeco = sapply(splitDeco, function(x) { deco=unlist(x)[2]; if (is.na(deco)) return ("") else return(paste0("_",deco)) })
    namecustom = make.unique(paste0("M",round(variableMetadata[,"mz"],numDigitsMZ),"T",round(variableMetadata[,"rt"],numDigitsRT),idsDeco))
    variableMetadata=cbind(name=variableMetadata$name, namecustom=namecustom, variableMetadata[,!(colnames(variableMetadata) %in% c("name"))])
    return(variableMetadata)
}

#The function annotateDiffreport without the corr function which bugs
annotatediff <- function(xset=xset, listArguments=listArguments, variableMetadataOutput="variableMetadata.tsv", dataMatrixOutput="dataMatrix.tsv") {
    # Resolve the bug with x11, with the function png
    options(bitmapType='cairo')

    #Check if the fillpeaks step has been done previously, if it hasn't, there is an error message and the execution is stopped.
    res=try(is.null(xset@filled))

    # ------ annot -------
    listArguments[["calcCiS"]]=as.logical(listArguments[["calcCiS"]])
    listArguments[["calcIso"]]=as.logical(listArguments[["calcIso"]])
    listArguments[["calcCaS"]]=as.logical(listArguments[["calcCaS"]])

    # common parameters
    listArguments4annotate = list(object=xset,
    nSlaves=listArguments[["nSlaves"]],sigma=listArguments[["sigma"]],perfwhm=listArguments[["perfwhm"]],
    maxcharge=listArguments[["maxcharge"]],maxiso=listArguments[["maxiso"]],minfrac=listArguments[["minfrac"]],
    ppm=listArguments[["ppm"]],mzabs=listArguments[["mzabs"]],quick=listArguments[["quick"]],
    polarity=listArguments[["polarity"]],max_peaks=listArguments[["max_peaks"]],intval=listArguments[["intval"]])

    # quick == FALSE
    if(listArguments[["quick"]]==FALSE) {
        listArguments4annotate = append(listArguments4annotate,
            list(graphMethod=listArguments[["graphMethod"]],cor_eic_th=listArguments[["cor_eic_th"]],pval=listArguments[["pval"]],
            calcCiS=listArguments[["calcCiS"]],calcIso=listArguments[["calcIso"]],calcCaS=listArguments[["calcCaS"]]))
        # no ruleset
        if (!is.null(listArguments[["multiplier"]])) {
            listArguments4annotate = append(listArguments4annotate,
                list(multiplier=listArguments[["multiplier"]]))
        }
        # ruleset
        else {
            rulset=read.table(listArguments[["rules"]], h=T, sep=";")
            if (ncol(rulset) < 4) rulset=read.table(listArguments[["rules"]], h=T, sep="\t")
            if (ncol(rulset) < 4) rulset=read.table(listArguments[["rules"]], h=T, sep=",")
            if (ncol(rulset) < 4) {
                error_message="Your ruleset file seems not well formatted. The column separators accepted are ; , and tabulation"
                print(error_message)
                stop(error_message)
            }

            listArguments4annotate = append(listArguments4annotate,
                list(rules=rulset))
        }
    }


    # launch annotate
    xa = do.call("annotate", listArguments4annotate)
    peakList=getPeaklist(xa,intval=listArguments[["intval"]])
    peakList=cbind(groupnames(xa@xcmsSet),peakList); colnames(peakList)[1] = c("name");

    # --- dataMatrix ---
    dataMatrix = peakList[,(make.names(colnames(peakList)) %in% c("name", make.names(sampnames(xa@xcmsSet))))]
    write.table(dataMatrix, sep="\t", quote=FALSE, row.names=FALSE, file=dataMatrixOutput)


    # --- Multi condition : diffreport ---
    diffrepOri=NULL
    if (!is.null(listArguments[["runDiffreport"]]) & nlevels(sampclass(xset))>=2) {
        #Check if the fillpeaks step has been done previously, if it hasn't, there is an error message and the execution is stopped.
        res=try(is.null(xset@filled))
        classes=levels(sampclass(xset))
        x=1:(length(classes)-1)
        for (i in seq(along=x) ) {
            y=1:(length(classes))
            for (n in seq(along=y)){
                if(i+n <= length(classes)){
                    filebase=paste(classes[i],class2=classes[i+n],sep="-vs-")

                    diffrep=diffreport(object=xset,class1=classes[i],class2=classes[i+n],filebase=filebase,eicmax=listArguments[["eicmax"]],eicwidth=listArguments[["eicwidth"]],sortpval=TRUE,value=listArguments[["value"]],h=listArguments[["h"]],w=listArguments[["w"]],mzdec=listArguments[["mzdec"]],missing=0)

                    diffrepOri = diffrep

                    # renamming of the column rtmed to rt to fit with camera peaklist function output
                    colnames(diffrep)[colnames(diffrep)=="rtmed"] <- "rt"
                    colnames(diffrep)[colnames(diffrep)=="mzmed"] <- "mz"

                    # combines results and reorder columns
                    diffrep = merge(peakList, diffrep[,c("name","fold","tstat","pvalue")], by.x="name", by.y="name", sort=F)
                    diffrep = cbind(diffrep[,!(colnames(diffrep) %in% c(sampnames(xa@xcmsSet)))],diffrep[,(colnames(diffrep) %in% c(sampnames(xa@xcmsSet)))])

                    diffrep = RTSecondToMinute(diffrep, listArguments[["convertRTMinute"]])
                    diffrep = formatIonIdentifiers(diffrep, numDigitsRT=listArguments[["numDigitsRT"]], numDigitsMZ=listArguments[["numDigitsMZ"]])

                    if(listArguments[["sortpval"]]){
                        diffrep=diffrep[order(diffrep$pvalue), ]
                    }

                    dir.create("tabular")
                    write.table(diffrep, sep="\t", quote=FALSE, row.names=FALSE, file=paste("tabular/",filebase,"_tsv.tabular",sep=""))

                    if (listArguments[["eicmax"]] != 0) {
                        diffreport_png2pdf(filebase)
                    }
                }
            }
        }
    }

    # --- variableMetadata ---
    variableMetadata=peakList[,!(make.names(colnames(peakList)) %in% c(make.names(sampnames(xa@xcmsSet))))]
    variableMetadata = RTSecondToMinute(variableMetadata, listArguments[["convertRTMinute"]])
    variableMetadata = formatIonIdentifiers(variableMetadata, numDigitsRT=listArguments[["numDigitsRT"]], numDigitsMZ=listArguments[["numDigitsMZ"]])
    # if we have 2 conditions, we keep stat of diffrep
    if (!is.null(listArguments[["runDiffreport"]]) & nlevels(sampclass(xset))==2) {
        variableMetadata = merge(variableMetadata, diffrep[,c("name","fold","tstat","pvalue")],by.x="name", by.y="name", sort=F)
        if(exists("listArguments[[\"sortpval\"]]")){
            variableMetadata=variableMetadata[order(variableMetadata$pvalue), ]
        }
    }

    variableMetadataOri=variableMetadata
    write.table(variableMetadata, sep="\t", quote=FALSE, row.names=FALSE, file=variableMetadataOutput)

    return(list("xa"=xa,"diffrep"=diffrepOri,"variableMetadata"=variableMetadataOri));

}


combinexsAnnos_function <- function(xaP, xaN, listOFlistArgumentsP,listOFlistArgumentsN, diffrepP=NULL,diffrepN=NULL,pos=TRUE,tol=2,ruleset=NULL,keep_meta=TRUE, convertRTMinute=F, numDigitsMZ=0, numDigitsRT=0, variableMetadataOutput="variableMetadata.tsv"){

    #Load the two Rdata to extract the xset objects from positive and negative mode
    cat("\tObject xset from positive mode\n")
    print(xaP)
    cat("\n")

    cat("\tObject xset from negative mode\n")
    print(xaN)
    cat("\n")

    cat("\n")
    cat("\tCombining...\n")
    #Convert the string to numeric for creating matrix
    row=as.numeric(strsplit(ruleset,",")[[1]][1])
    column=as.numeric(strsplit(ruleset,",")[[1]][2])
    ruleset=cbind(row,column)
    #Test if the file comes from an older version tool
    if ((!is.null(xaP)) & (!is.null(xaN))) {
        #Launch the combinexsannos function from CAMERA
        cAnnot=combinexsAnnos(xaP, xaN,pos=pos,tol=tol,ruleset=ruleset)
    } else {
        stop("You must relauch the CAMERA.annotate step with the lastest version.")
    }

    if(pos){
        xa=xaP
        listOFlistArgumentsP=listOFlistArguments
        mode="neg. Mode"
    } else {
        xa=xaN
        listOFlistArgumentsN=listOFlistArguments
        mode="pos. Mode"
    }

    peakList=getPeaklist(xa)
    peakList=cbind(groupnames(xa@xcmsSet),peakList); colnames(peakList)[1] = c("name");
    variableMetadata=cbind(peakList, cAnnot[, c("isotopes", "adduct", "pcgroup",mode)]);
    variableMetadata=variableMetadata[,!(colnames(variableMetadata) %in% c(sampnames(xa@xcmsSet)))]

    #Test if there are more than two classes (conditions)
    if ( nlevels(sampclass(xaP@xcmsSet))==2 & (!is.null(diffrepN)) & (!is.null(diffrepP))) {
        diffrepP = diffrepP[,c("name","fold","tstat","pvalue")]; colnames(diffrepP) = paste("P.",colnames(diffrepP),sep="")
        diffrepN = diffrepN[,c("name","fold","tstat","pvalue")]; colnames(diffrepN) = paste("N.",colnames(diffrepN),sep="")

        variableMetadata = merge(variableMetadata, diffrepP, by.x="name", by.y="P.name")
        variableMetadata = merge(variableMetadata, diffrepN, by.x="name", by.y="N.name")
    }

    rownames(variableMetadata) = NULL
    #TODO: checker
    #colnames(variableMetadata)[1:2] = c("name","mz/rt");

    variableMetadata = RTSecondToMinute(variableMetadata, convertRTMinute)
    variableMetadata = formatIonIdentifiers(variableMetadata, numDigitsRT=numDigitsRT, numDigitsMZ=numDigitsMZ)

    #If the user want to keep only the metabolites which match a difference
    if(keep_meta){
        variableMetadata=variableMetadata[variableMetadata[,c(mode)]!="",]
    }

    #Write the output into a tsv file
    write.table(variableMetadata, sep="\t", quote=FALSE, row.names=FALSE, file=variableMetadataOutput)
    return(variableMetadata);

}

# This function get the raw file path from the arguments
getRawfilePathFromArguments <- function(singlefile, zipfile, listArguments) {
    if (!is.null(listArguments[["zipfile"]]))           zipfile = listArguments[["zipfile"]]
    if (!is.null(listArguments[["zipfilePositive"]]))   zipfile = listArguments[["zipfilePositive"]]
    if (!is.null(listArguments[["zipfileNegative"]]))   zipfile = listArguments[["zipfileNegative"]]

    if (!is.null(listArguments[["singlefile_galaxyPath"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPath"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleName"]]
    }
    if (!is.null(listArguments[["singlefile_galaxyPathPositive"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPathPositive"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleNamePositive"]]
    }
    if (!is.null(listArguments[["singlefile_galaxyPathNegative"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPathNegative"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleNameNegative"]]
    }
    if (exists("singlefile_galaxyPaths")){
        singlefile_galaxyPaths = unlist(strsplit(singlefile_galaxyPaths,","))
        singlefile_sampleNames = unlist(strsplit(singlefile_sampleNames,","))

        singlefile=NULL
        for (singlefile_galaxyPath_i in seq(1:length(singlefile_galaxyPaths))) {
            singlefile_galaxyPath=singlefile_galaxyPaths[singlefile_galaxyPath_i]
            singlefile_sampleName=singlefile_sampleNames[singlefile_galaxyPath_i]
            singlefile[[singlefile_sampleName]] = singlefile_galaxyPath
        }
    }
    for (argument in c("zipfile","zipfilePositive","zipfileNegative","singlefile_galaxyPath","singlefile_sampleName","singlefile_galaxyPathPositive","singlefile_sampleNamePositive","singlefile_galaxyPathNegative","singlefile_sampleNameNegative")) {
        listArguments[[argument]]=NULL
    }
    return(list(zipfile=zipfile, singlefile=singlefile, listArguments=listArguments))
}


# This function retrieve the raw file in the working directory
#   - if zipfile: unzip the file with its directory tree
#   - if singlefiles: set symlink with the good filename
retrieveRawfileInTheWorkingDirectory <- function(singlefile, zipfile) {
    if(!is.null(singlefile) && (length("singlefile")>0)) {
        for (singlefile_sampleName in names(singlefile)) {
            singlefile_galaxyPath = singlefile[[singlefile_sampleName]]
            if(!file.exists(singlefile_galaxyPath)){
                error_message=paste("Cannot access the sample:",singlefile_sampleName,"located:",singlefile_galaxyPath,". Please, contact your administrator ... if you have one!")
                print(error_message); stop(error_message)
            }

            file.symlink(singlefile_galaxyPath,singlefile_sampleName)
        }
        directory = "."

    }
    if(!is.null(zipfile) && (zipfile!="")) {
        if(!file.exists(zipfile)){
            error_message=paste("Cannot access the Zip file:",zipfile,". Please, contact your administrator ... if you have one!")
            print(error_message)
            stop(error_message)
        }

        #list all file in the zip file
        #zip_files=unzip(zipfile,list=T)[,"Name"]

        #unzip
        suppressWarnings(unzip(zipfile, unzip="unzip"))

        #get the directory name
        filesInZip=unzip(zipfile, list=T);
        directories=unique(unlist(lapply(strsplit(filesInZip$Name,"/"), function(x) x[1])));
        directories=directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
        directory = "."
        if (length(directories) == 1) directory = directories

        cat("files_root_directory\t",directory,"\n")

    }
    return (directory)
}

#@TODO: remove this function as soon as we can use xcms 3.x.x from Bioconductor 3.7
# https://github.com/sneumann/CAMERA/issues/33#issuecomment-405168524
# https://github.com/sneumann/xcms/commit/950a3fe794cdb6b0fda88696e31aab3d97a3b7dd
############################################################
## getEIC
getEIC <- function(object, mzrange, rtrange = 200,
                   groupidx, sampleidx = sampnames(object),
                   rt = c("corrected", "raw")) {
    
    files <- filepaths(object)
    grp <- groups(object)
    samp <- sampnames(object)
    prof <- profinfo(object)
    
    rt <- match.arg(rt)
    
    if (is.numeric(sampleidx))
    sampleidx <- sampnames(object)[sampleidx]
    sampidx <- match(sampleidx, sampnames(object))
    
    if (!missing(groupidx)) {
        if (is.numeric(groupidx))
        groupidx <- groupnames(object)[unique(as.integer(groupidx))]
        grpidx <- match(groupidx, groupnames(object, template = groupidx))
    }
    
    if (missing(mzrange)) {
        if (missing(groupidx))
        stop("No m/z range or groups specified")
        if (any(is.na(groupval(object, value = "mz"))))
        warning(
        "`NA` values in xcmsSet. Use fillPeaks() on the object to fill",
        "-in missing peak values. Note however that this will also ",
        "insert intensities of 0 for peaks that can not be filled in.")
        mzmin <- apply(groupval(object, value = "mzmin"), 1, min, na.rm = TRUE)
        mzmax <- apply(groupval(object, value = "mzmax"), 1, max, na.rm = TRUE)
        mzrange <- matrix(c(mzmin[grpidx], mzmax[grpidx]), ncol = 2)
        ## if (any(is.na(groupval(object, value = "mz"))))
        ##     stop('Please use fillPeaks() to fill up NA values !')
        ## mzmin <- -rowMax(-groupval(object, value = "mzmin"))
        ## mzmax <- rowMax(groupval(object, value = "mzmax"))
        ## mzrange <- matrix(c(mzmin[grpidx], mzmax[grpidx]), ncol = 2)
    } else if (all(c("mzmin","mzmax") %in% colnames(mzrange)))
    mzrange <- mzrange[,c("mzmin", "mzmax"),drop=FALSE]
    else if (is.null(dim(mzrange)))
    stop("mzrange must be a matrix")
    colnames(mzrange) <- c("mzmin", "mzmax")
    
    if (length(rtrange) == 1) {
        if (missing(groupidx))
        rtrange <- matrix(rep(range(object@rt[[rt]][sampidx]), nrow(mzrange)),
        ncol = 2, byrow = TRUE)
        else {
            rtrange <- retexp(grp[grpidx,c("rtmin","rtmax"),drop=FALSE], rtrange)
        }
    } else if (is.null(dim(rtrange)))
    stop("rtrange must be a matrix or single number")
    colnames(rtrange) <- c("rtmin", "rtmax")
    
    ## Ensure that we've got corrected retention time if requested.
    if (is.null(object@rt[[rt]]))
    stop(rt, " retention times not present in 'object'!")
    
    ## Ensure that the defined retention time range is within the rtrange of the
    ## object: we're using the max minimal rt of all files and the min maximal rt
    rtrs <- lapply(object@rt[[rt]], range)
    rtr <- c(max(unlist(lapply(rtrs, "[", 1))),
    min(unlist(lapply(rtrs, "[", 2))))
    ## Check if we've got a range which is completely off:
    if (any(rtrange[, "rtmin"] >= rtr[2] | rtrange[, "rtmax"] <= rtr[1])) {
        outs <- which(rtrange[, "rtmin"] >= rtr[2] |
        rtrange[, "rtmax"] <= rtr[1])
        stop(length(outs), " of the specified 'rtrange' are completely outside ",
        "of the retention time range of 'object' which is (", rtr[1], ", ",
        rtr[2], "). The first was: (", rtrange[outs[1], "rtmin"], ", ",
        rtrange[outs[1], "rtmax"], "!")
    }
    lower_rt_outside <- rtrange[, "rtmin"] < rtr[1]
    upper_rt_outside <- rtrange[, "rtmax"] > rtr[2]
    if (any(lower_rt_outside) | any(upper_rt_outside)) {
        ## Silently fix these ranges.
        rtrange[lower_rt_outside, "rtmin"] <- rtr[1]
        rtrange[upper_rt_outside, "rtmax"] <- rtr[2]
    }
    
    if (missing(groupidx))
    gnames <- character(0)
    else
    gnames <- groupidx
    
    eic <- vector("list", length(sampleidx))
    names(eic) <- sampleidx
    
    for (i in seq(along = sampidx)) {
        
        ## cat(sampleidx[i], "")
        flush.console()
        ## getXcmsRaw takes care of rt correction, susetting to scanrage and other
        ## stuff.
        lcraw <- getXcmsRaw(object, sampleidx = sampidx[i], rt=rt)
        currenteic <- xcms::getEIC(lcraw, mzrange, rtrange, step = prof$step)
        eic[[i]] <- currenteic@eic[[1]]
        rm(lcraw)
        gc()
    }
    ## cat("\n")
    
    invisible(new("xcmsEIC", eic = eic, mzrange = mzrange, rtrange = rtrange,
    rt = rt, groupnames = gnames))
}

#@TODO: remove this function as soon as we can use xcms 3.x.x from Bioconductor 3.7
# https://github.com/sneumann/CAMERA/issues/33#issuecomment-405168524
# https://github.com/sneumann/xcms/commit/950a3fe794cdb6b0fda88696e31aab3d97a3b7dd
############################################################
## diffreport
diffreport = function(object,
                      class1 = levels(sampclass(object))[1],
                      class2 = levels(sampclass(object))[2],
                      filebase = character(),
                      eicmax = 0, eicwidth = 200,
                      sortpval = TRUE,
                      classeic = c(class1,class2),
                      value = c("into","maxo","intb"),
                      metlin = FALSE,
                      h = 480, w = 640, mzdec=2,
                      missing = numeric(), ...) {
    
    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }
    
    if (!is.numeric(w) || !is.numeric(h))
        stop("'h' and 'w' have to be numeric")
    ## require(multtest) || stop("Couldn't load multtest")
    
    value <- match.arg(value)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
    stop("No group information found")
    samples <- sampnames(object)
    n <- length(samples)
    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]
    
    values <- groupval(object, "medret", value=value)
    indecies <- groupval(object, "medret", value = "index")
    
    if (!all(c(class1,class2) %in% classlabel))
    stop("Incorrect Class Labels")
    
    ## c1 and c2 are column indices of class1 and class2 resp.
    c1 <- which(classlabel %in% class1)
    c2 <- which(classlabel %in% class2)
    ceic <- which(classlabel %in% classeic)
    if (length(intersect(c1, c2)) > 0)
    stop("Intersecting Classes")
    
    ## Optionally replace NA values with the value provided with missing
    if (length(missing)) {
        if (is.numeric(missing)) {
            ## handles NA, Inf and -Inf
            values[, c(c1, c2)][!is.finite(values[, c(c1, c2)])] <- missing[1]
        } else
        stop("'missing' should be numeric")
    }
    ## Check against missing Values
    if (any(is.na(values[, c(c1, c2)])))
    warning("`NA` values in xcmsSet. Use fillPeaks() on the object to fill",
    "-in missing peak values. Note however that this will also ",
    "insert intensities of 0 for peaks that can not be filled in.")
    
    mean1 <- rowMeans(values[,c1,drop=FALSE], na.rm = TRUE)
    mean2 <- rowMeans(values[,c2,drop=FALSE], na.rm = TRUE)
    
    ## Calculate fold change.
    ## For foldchange <1 set fold to 1/fold
    ## See tstat to check which was higher
    fold <- mean2 / mean1
    fold[!is.na(fold) & fold < 1] <- 1/fold[!is.na(fold) & fold < 1]
    
    testval <- values[,c(c1,c2)]
    ## Replace eventual infinite values with NA (CAMERA issue #33)
    testval[is.infinite(testval)] <- NA
    testclab <- c(rep(0,length(c1)),rep(1,length(c2)))
    
    if (min(length(c1), length(c2)) >= 2) {
        tstat <- mt.teststat(testval, testclab, ...)
        pvalue <- xcms:::pval(testval, testclab, tstat)
    } else {
        message("Too few samples per class, skipping t-test.")
        tstat <- pvalue <- rep(NA,nrow(testval))
    }
    stat <- data.frame(fold = fold, tstat = tstat, pvalue = pvalue)
    if (length(levels(sampclass(object))) >2) {
        pvalAnova<-c()
        for(i in 1:nrow(values)){
            var<-as.numeric(values[i,])
            ano<-summary(aov(var ~ sampclass(object)) )
            pvalAnova<-append(pvalAnova, unlist(ano)["Pr(>F)1"])
        }
        stat<-cbind(stat, anova= pvalAnova)
    }
    if (metlin) {
        neutralmass <- groupmat[,"mzmed"] + ifelse(metlin < 0, 1, -1)
        metlin <- abs(metlin)
        digits <- ceiling(-log10(metlin))+1
        metlinurl <-
        paste("http://metlin.scripps.edu/simple_search_result.php?mass_min=",
        round(neutralmass - metlin, digits), "&mass_max=",
        round(neutralmass + metlin, digits), sep="")
        values <- cbind(metlin = metlinurl, values)
    }
    twosamp <- cbind(name = groupnames(object), stat, groupmat, values)
    if (sortpval) {
        tsidx <- order(twosamp[,"pvalue"])
        twosamp <- twosamp[tsidx,]
        rownames(twosamp) <- 1:nrow(twosamp)
        values<-values[tsidx,]
    } else
    tsidx <- 1:nrow(values)
    
    if (length(filebase))
    write.table(twosamp, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)
    
    if (eicmax > 0) {
        if (length(unique(peaks(object)[,"rt"])) > 1) {
            ## This looks like "normal" LC data
            
            eicmax <- min(eicmax, length(tsidx))
            eics <- getEIC(object, rtrange = eicwidth*1.1, sampleidx = ceic,
            groupidx = tsidx[seq(length = eicmax)])
            
            if (length(filebase)) {
                eicdir <- paste(filebase, "_eic", sep="")
                boxdir <- paste(filebase, "_box", sep="")
                dir.create(eicdir)
                dir.create(boxdir)
                if (capabilities("png")){
                    xcms:::xcmsBoxPlot(values[seq(length = eicmax),],
                    sampclass(object), dirpath=boxdir, pic="png",  width=w, height=h)
                    png(file.path(eicdir, "%003d.png"), width = w, height = h)
                } else {
                    xcms:::xcmsBoxPlot(values[seq(length = eicmax),],
                    sampclass(object), dirpath=boxdir, pic="pdf", width=w, height=h)
                    pdf(file.path(eicdir, "%003d.pdf"), width = w/72,
                    height = h/72, onefile = FALSE)
                }
            }
            plot(eics, object, rtrange = eicwidth, mzdec=mzdec)
            
            if (length(filebase))
            dev.off()
        } else {
            ## This looks like a direct-infusion single spectrum
            if (length(filebase)) {
                eicdir <- paste(filebase, "_eic", sep="")
                boxdir <- paste(filebase, "_box", sep="")
                dir.create(eicdir)
                dir.create(boxdir)
                if (capabilities("png")){
                    xcmsBoxPlot(values[seq(length = eicmax),],
                    sampclass(object), dirpath=boxdir, pic="png",
                    width=w, height=h)
                    png(file.path(eicdir, "%003d.png"), width = w, height = h,
                    units = "px")
                } else {
                    xcmsBoxPlot(values[seq(length = eicmax),],
                    sampclass(object), dirpath=boxdir, pic="pdf",
                    width=w, height=h)
                    pdf(file.path(eicdir, "%003d.pdf"), width = w/72,
                    height = h/72, onefile = FALSE)
                }
            }
            
            plotSpecWindow(object, gidxs = tsidx[seq(length = eicmax)], borderwidth=1)
            
            if (length(filebase))
            dev.off()
        }
    }
    
    invisible(twosamp)
}
