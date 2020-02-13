# lib.r

#@author G. Le Corguille
# solve an issue with batch if arguments are logical TRUE/FALSE
parseCommandArgs <- function(...) {
    args <- batch::parseCommandArgs(...)
    for (key in names(args)) {
        if (args[key] %in% c("TRUE","FALSE"))
            args[key] = as.logical(args[key])
    }
    return(args)
}

#@author G. Le Corguille
# This function will
# - load the packages
# - display the sessionInfo
loadAndDisplayPackages <- function(pkgs) {
    for(pkg in pkgs) suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))

    sessioninfo = sessionInfo()
    cat(sessioninfo$R.version$version.string,"\n")
    cat("Main packages:\n")
    for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
    cat("Other loaded packages:\n")
    for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
}

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
#The function create a zip archive from the different png generated by diffreport
diffreport_png2zip <- function() {
    zip("eic.zip", dir(pattern="_eic"), zip=Sys.which("zip"))
    zip("box.zip", dir(pattern="_box"), zip=Sys.which("zip"))
}

#The function create a zip archive from the different tabular generated by diffreport
diffreport_tabular2zip <- function() {
    zip("tabular.zip", dir(pattern="tabular/*"), zip=Sys.which("zip"))
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
annotatediff <- function(xset=xset, args=args, variableMetadataOutput="variableMetadata.tsv") {
    # Resolve the bug with x11, with the function png
    options(bitmapType='cairo')

    #Check if the fillpeaks step has been done previously, if it hasn't, there is an error message and the execution is stopped.
    res=try(is.null(xset@filled))

    # ------ annot -------
    args$calcCiS=as.logical(args$calcCiS)
    args$calcIso=as.logical(args$calcIso)
    args$calcCaS=as.logical(args$calcCaS)

    # common parameters
    args4annotate = list(object=xset,
        nSlaves=args$nSlaves,sigma=args$sigma,perfwhm=args$perfwhm,
        maxcharge=args$maxcharge,maxiso=args$maxiso,minfrac=args$minfrac,
        ppm=args$ppm,mzabs=args$mzabs,quick=args$quick,
        polarity=args$polarity,max_peaks=args$max_peaks,intval=args$intval)

    # quick == FALSE
    if(args$quick==FALSE) {
        args4annotate = append(args4annotate,
            list(graphMethod=args$graphMethod,cor_eic_th=args$cor_eic_th,pval=args$pval,
            calcCiS=args$calcCiS,calcIso=args$calcIso,calcCaS=args$calcCaS))
        # no ruleset
        if (!is.null(args$multiplier)) {
            args4annotate = append(args4annotate,
                list(multiplier=args$multiplier))
        }
        # ruleset
        else {
            rulset=read.table(args$rules, h=T, sep=";")
            if (ncol(rulset) < 4) rulset=read.table(args$rules, h=T, sep="\t")
            if (ncol(rulset) < 4) rulset=read.table(args$rules, h=T, sep=",")
            if (ncol(rulset) < 4) {
                error_message="Your ruleset file seems not well formatted. The column separators accepted are ; , and tabulation"
                print(error_message)
                stop(error_message)
            }

            args4annotate = append(args4annotate,
                list(rules=rulset))
        }
    }


    # launch annotate
    xa = do.call("annotate", args4annotate)
    peakList=getPeaklist(xa,intval=args$intval)
    peakList=cbind(groupnames(xa@xcmsSet),peakList); colnames(peakList)[1] = c("name");

    # --- Multi condition : diffreport ---
    diffrepOri=NULL
    if (!is.null(args$runDiffreport) & nlevels(sampclass(xset))>=2) {
        #Check if the fillpeaks step has been done previously, if it hasn't, there is an error message and the execution is stopped.
        res=try(is.null(xset@filled))
        classes=levels(sampclass(xset))
        x=1:(length(classes)-1)
        for (i in seq(along=x) ) {
            y=1:(length(classes))
            for (n in seq(along=y)){
                if(i+n <= length(classes)){
                    filebase=paste(classes[i],class2=classes[i+n],sep="-vs-")

                    diffrep=diffreport(
                        object=xset,class1=classes[i],class2=classes[i+n],
                        filebase=filebase,eicmax=args$eicmax,eicwidth=args$eicwidth,
                        sortpval=TRUE,value=args$value,h=args$h,w=args$w,mzdec=args$mzdec,missing=0)

                    diffrepOri = diffrep

                    # renamming of the column rtmed to rt to fit with camera peaklist function output
                    colnames(diffrep)[colnames(diffrep)=="rtmed"] <- "rt"
                    colnames(diffrep)[colnames(diffrep)=="mzmed"] <- "mz"

                    # combines results and reorder columns
                    diffrep = merge(peakList, diffrep[,c("name","fold","tstat","pvalue")], by.x="name", by.y="name", sort=F)
                    diffrep = cbind(diffrep[,!(colnames(diffrep) %in% c(sampnames(xa@xcmsSet)))],diffrep[,(colnames(diffrep) %in% c(sampnames(xa@xcmsSet)))])

                    diffrep = RTSecondToMinute(diffrep, args$convertRTMinute)
                    diffrep = formatIonIdentifiers(diffrep, numDigitsRT=args$numDigitsRT, numDigitsMZ=args$numDigitsMZ)

                    if(args$sortpval){
                        diffrep=diffrep[order(diffrep$pvalue), ]
                    }

                    dir.create("tabular", showWarnings = FALSE)
                    write.table(diffrep, sep="\t", quote=FALSE, row.names=FALSE, file=paste("tabular/",filebase,"_tsv.tabular",sep=""))

                    if (args$eicmax != 0) {
                        if (args$png2 == "pdf")
                            diffreport_png2pdf(filebase)
                    }
                }
            }
        }
        if (args$png2 == "zip")
            diffreport_png2zip()
        if (args$tabular2 == "zip")
            diffreport_tabular2zip()
    }

    # --- variableMetadata ---
    variableMetadata=peakList[,!(make.names(colnames(peakList)) %in% c(make.names(sampnames(xa@xcmsSet))))]
    variableMetadata = RTSecondToMinute(variableMetadata, args$convertRTMinute)
    variableMetadata = formatIonIdentifiers(variableMetadata, numDigitsRT=args$numDigitsRT, numDigitsMZ=args$numDigitsMZ)
    # if we have 2 conditions, we keep stat of diffrep
    if (!is.null(args$runDiffreport) & nlevels(sampclass(xset))==2) {
        variableMetadata = merge(variableMetadata, diffrep[,c("name","fold","tstat","pvalue")],by.x="name", by.y="name", sort=F)
        if(exists("args[[\"sortpval\"]]")){
            variableMetadata=variableMetadata[order(variableMetadata$pvalue), ]
        }
    }

    variableMetadataOri=variableMetadata
    write.table(variableMetadata, sep="\t", quote=FALSE, row.names=FALSE, file=variableMetadataOutput)

    return(list("xa"=xa,"diffrep"=diffrepOri,"variableMetadata"=variableMetadataOri));

}


combinexsAnnos_function <- function(xaP, xaN, diffrepP=NULL,diffrepN=NULL,
    pos=TRUE,tol=2,ruleset=NULL,keep_meta=TRUE, convertRTMinute=F, numDigitsMZ=0,
    numDigitsRT=0, variableMetadataOutput="variableMetadata.tsv"){

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
        mode="neg. Mode"
    } else {
        xa=xaN
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
getRawfilePathFromArguments <- function(singlefile, zipfile, args) {
    if (!is.null(args$zipfile))           zipfile = args$zipfile
    if (!is.null(args$zipfilePositive))   zipfile = args$zipfilePositive
    if (!is.null(args$zipfileNegative))   zipfile = args$zipfileNegative

    if (!is.null(args$singlefile_galaxyPath)) {
        singlefile_galaxyPaths = args$singlefile_galaxyPath;
        singlefile_sampleNames = args$singlefile_sampleName
    }
    if (!is.null(args$singlefile_galaxyPathPositive)) {
        singlefile_galaxyPaths = args$singlefile_galaxyPathPositive;
        singlefile_sampleNames = args$singlefile_sampleNamePositive
    }
    if (!is.null(args$singlefile_galaxyPathNegative)) {
        singlefile_galaxyPaths = args$singlefile_galaxyPathNegative;
        singlefile_sampleNames = args$singlefile_sampleNameNegative
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
    for (argument in c("zipfile", "zipfilePositive", "zipfileNegative",
                        "singlefile_galaxyPath", "singlefile_sampleName",
                        "singlefile_galaxyPathPositive", "singlefile_sampleNamePositive",
                        "singlefile_galaxyPathNegative","singlefile_sampleNameNegative")) {
        args[[argument]]=NULL
    }
    return(list(zipfile=zipfile, singlefile=singlefile, args=args))
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

