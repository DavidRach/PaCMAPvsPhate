#' Runs PaCMAP via Basilisk
#'
#' @param x A flowframe object
#' @param subset Default root
#' @param columns A subset of columns to pass instead
#' @param notcolumns A subset of columns to remove
#' @param outpath Location to store new .fcs files
#' @param returnType Default fcs
#' @param NewName Default PaCMAP
#'
#' @importFrom basilisk basiliskStart
#' @importFrom basilisk basiliskStop
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom Biobase exprs
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom basilisk basiliskRun
#' @importFrom flowCore write.FCS
#' @importFrom reticulate import
#' @importFrom reticulate py_to_r
#' @importFrom Luciernaga Utility_ColAppend
#' @importFrom magrittr %>%
#'
#' @return Two columns for the Dimensionality visualized data
#' @export
#'
#' @examples
#'
#' A <- 2 + 2
Utility_PaCMAP <- function(x, subset="root", columns=NULL, notcolumns=NULL,
    returnType="fcs", outpath=NULL, NewName="PaCMAP"){
  
    proc <- basiliskStart(env1)
    on.exit(basiliskStop(proc))
  
    ff <- gs_pop_get_data(x, subset, inverse.transform=FALSE)
    
    df <- exprs(ff[[1]])
    DF <- as.data.frame(df, check.names = FALSE)

    #Consolidating Columns Going Forward
    CleanedDF <- DF[,-grep("Time|FS|SC|SS|Original|W$|H$|AF|Coereba", names(DF))]

    if (!is.null(columns) && !is.null(notcolumns)) {
      stop("Columns and notcolumns are not currently combinable. Pick one")}

    # If external columns interest specified
    if (!is.null(columns)){CleanedDF1 <- CleanedDF %>% select(all_of(columns))
    } else {CleanedDF1 <- CleanedDF}

    if (!is.null(notcolumns)){CleanedDF1 <- CleanedDF1 %>% select(-all_of(notcolumns))
    } else {CleanedDF1 <- CleanedDF1}

    X <- CleanedDF1

    ThePaCMAP <- basiliskRun(proc, .Internal_PaCMAP, X=X)
  
    newff <- gs_pop_get_data(x, subset, inverse.transform=TRUE)
    df <- exprs(newff[[1]])
    DF <- as.data.frame(df, check.names = FALSE)
  
    TheFileName <- paste0(NewName, ".fcs")
  
    if (is.null(outpath)){outpath <- getwd()}

    fileSpot <- file.path(outpath, TheFileName)
  
    new_fcs <- Utility_ColAppend(ff=newff, DF=DF, columnframe=ThePaCMAP, shift = TRUE)
  
    new_fcs@description$GUID <- TheFileName
    new_fcs@description$`$FIL` <- TheFileName

    if (returnType=="fcs") {write.FCS(new_fcs, filename = fileSpot, delimiter="#")
    } else {return(new_fcs)}
}

.Internal_PaCMAP <-  function(X){
  pacmap <- import("pacmap")
  np <- import("numpy")
  pd <- import("pandas", convert = TRUE)

  # Converting r.x to a pandas DataFrame
  df <- pd$DataFrame(X)
  X <- np$asarray(df)

  # Initializing the pacmap instance
  embedding <- pacmap$PaCMAP(n_components = as.integer(2), n_neighbors=as.integer(15), MN_ratio=0.5, FP_ratio=2.0)

  # Fit the data (The index of transformed data corresponds to the index of the original data)
  X_transformed <- embedding$fit_transform(X, init="pca")

  colnames(X_transformed) <- c("PaCMAP_1", "PaCMAP_2")
  ThePaCMAP <- data.frame(X_transformed)
}
