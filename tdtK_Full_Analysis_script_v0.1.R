###### Complete Rscript to run TdtK Analysis
## Run the entire script from within R or Rstudio
######

### NOTES ########## 
#  If ImageJ/Fiji and/or showinf are NOT found please follows these steps:
#  Specifically for RStudio (Linux) the path to ImageJ-linux64 needs to be set
#  by replacing "/home/geo/Fiji.app/" with the path to YOUR ImageJ-linux64
#  executable. ONLY necessary if the script can't find ImageJ from within
#  Rstudio (but it would work in terminal mode)
#  > Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/home/geo/Fiji.app/", sep=":"))

#  Another example for MacOS (if ImageJ-macosx resides in "/Applications/Fiji.app/Contents/macos/"):
#  > Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/Applications/Fiji.app/Contents/macos/", sep=":"))

#  For Windows, edit System Environment Variables and the location of the ImageJ-win64.exe and showinf.bat
#  to the PATH variable.

# Hard-coded locations if interactive mode is not required (uncomment and replace with relevant paths)
#
# answer_info <- FALSE
# script_dir <- c("~/Rserver_tdtK")
# movie_dir <- c("/home/geo/bacon_mount/Clara/tdtK-data/EpigeneticScreen-7wk/")
# target_dir <- c("~/GIANT/balled_MAYO_only")
# mappings_file <- c("~/GIANT/balled_MAYO_only/excellent traces/mappings.xlsx")
# mappings_file <- c("~/GIANT/balled_others/good traces/mappings.xlsx")

{
  answer_info <- rstudioapi::showQuestion("Enter folder and file?", "\"Yes\" will ask for files and folder locations.", ok = "Yes", cancel = "No")
  
  if (answer_info == TRUE | exists("OS") == FALSE) {
    # Check whether Fiji/ImageJ and bftools are installed
    if (.Platform$OS.type == "windows") { 
      showninf_installed <- system('showinf.bat -version') == 0
      ImageJ_installed <- system('ImageJ-win64.exe --ij2 --headless') == 0
      OS <- "windows"
      
    } else if (Sys.info()["sysname"] == "Darwin") {
      showninf_installed <- system('showinf -version') == 0
      ImageJ_installed <- system('ImageJ-macosx --ij2 --headless') == 0
      OS <- "Darwin"
      
    } else if (.Platform$OS.type == "unix") { 
      showninf_installed <- system('showinf -version') == 0
      ImageJ_installed <- system('ImageJ-linux64 --ij2 --headless') == 0
      OS <- "unix"
      
    } else if (.Platform$OS.type == "Linux") { 
      showninf_installed <- system('showinf -version') == 0
      ImageJ_installed <- system('ImageJ-linux64 --ij2 --headless') == 0
      OS <- "unix"    
      
    } else {
      showninf_installed <- "FALSE"
      ImageJ_installed <- "FALSE"
    }
    
    ifelse(showninf_installed == TRUE & ImageJ_installed == TRUE, paste0("Everything is installed"),
           ifelse(showninf_installed == FALSE & ImageJ_installed == FALSE, stop(paste0("bftools and Fiji/ImageJ are missing !")),
                  ifelse(showninf_installed == FALSE & ImageJ_installed == TRUE, stop(paste0("bftools is missing !")), stop(paste0("ImageJ/Fiji is missing !")))))
    
    
    # Where are the script files?
    script_dir <- rstudioapi::selectDirectory(caption = "Select Directory containing Rscript Files", label = "Rscript Files" )
    IJscript <- normalizePath(file.path(script_dir,"macro.ijm"))
    
    if(file.exists(IJscript) == FALSE)
    {
      stop("This folder does not contain the required file(s)")
    }
    
    # Where are the movie files?
    movie_dir <-rstudioapi::selectDirectory(caption = "Select Directory containing tdtK Movie Files", label = "tdtK Movie Files" )
    
    # Target directory - this is where all processed files will be analyzed
    target_dir <- rstudioapi::selectDirectory(caption = "Select Directory for Final Processing", label = "Processing Folder" )
    
    # Point to the Mappings file
    mappings_file <-  rstudioapi::selectFile(caption = "Select Genotype Mappings File", label = "mappings.xlsx -or- .csv")
    
    if (basename(mappings_file) == "mappings.xlsx")
    {
      mappings_file_type <- "XLSX"
    } 
    
    if (basename(mappings_file) == "mappings.csv")
    {
      mappings_file_type <- "CSV"
    } 
    
    if (basename(mappings_file) != "mappings.xlsx" && basename(mappings_file) != "mappings.csv")
    {
      
      stop("The mappings file is not named correctly")
    }
  }
  
  
  # Run the first part of the script - this creates kymographs from each CXD file
  # MAYO Screen Script No.1
  # Imports CXD movies of tdtk hearts and writes TIFF files of Kymographs and overview images
  
  # Required libraries and functions
  options(java.parameters = "-Xmx11g" )
  
  library(EBImage)
  library(matrixStats)
  library(RBioFormats)
  library(pracma)
  library(zoo)
  library(stringr)
  library(reshape2)
  library(plyr)
  
  
  # Define contrast function
  RMSE <- function(m, o)  {
    sqrt(mean((m - o)^2))
  }
  
  # Define Slope function
  rollingSlope.lm.fit <- function(vector) {
    
    a <- coef(.lm.fit(cbind(1, seq(vector)), vector))[2]
    return(a)
  }
  
  # New find peaks from https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
  find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
      z <- i - m + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + m + 1
      w <- ifelse(w < length(x), w, length(x))
      if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks <- unlist(pks)
    pks
  }
  
  # And define a different find.max function
  which.maxN <- function(x, N=2){ #N=2 means find the 2nd max
    len <- length(x)
    if(N>len){
      warning('N greater than length(x).  Setting N=length(x)')
      N <- length(x)
    }
    which(x == sort(x,partial=len-N+1)[len-N+1])
  }
  
  
  # Switch directory containing CXD files
  setwd(movie_dir)
  
  # Get file list, for all necessary file types (.cxd)
  
  filelist <- NULL
  filelist$cxd <- list.files(".", pattern = "\\.cxd$", recursive = TRUE)
  
  total <- length(filelist$cxd)
  pb <- txtProgressBar(max = total, style = 3)
  
  for (f in 1:length(filelist$cxd))
  {
    setTxtProgressBar(pb, f)
    
    # Check if tiff files have already been generated - skip to next file
    tiffs <- list.files(path = dirname(filelist$cxd[f]), pattern = paste0(basename(filelist$cxd[f]),'_peak_'))
    if(length(tiffs) >= 1)
    {
      next
    }
    
    
    if(file.size(filelist$cxd[f]) < 150000000) # OS-independent version
    {
      next
    }   
    
    # Check if metadata file is present
    
    if(file.exists(paste0(substr(filelist$cxd[f],1,gregexpr('.cxd',filelist$cxd[f])[[1]][1]),"_new_meta_data.csv")))
    {
      meta.data <- read.csv(paste0(substr(filelist$cxd[f],1,gregexpr('.cxd',filelist$cxd[f])[[1]][1]),"_new_meta_data.csv"))
    } else
    {
      # Read metadata of cxd file and store separately for later use. Also,
      # check if file can be imported or is defective and skip if so
      
      if(Sys.info()['sysname'] == "Windows"){
        meta.data <- as.data.frame(system(paste0('showinf.bat -nopix ', shQuote(filelist$cxd[f])), intern = TRUE), stringsAsFactors = FALSE)
      } else {
        meta.data <- as.data.frame(system(paste0("showinf -nopix ", shQuote(filelist$cxd[f])), intern = TRUE), stringsAsFactors = FALSE)
      }
      
      meta.data[,1] <- gsub("\t", "", meta.data[,1])
      meta.data <- as.data.frame(str_split_fixed(meta.data[,1], "[:=]", 2), stringsAsFactors = FALSE)
      
      # Correct for missing date/time info in file
      test <- tryCatch(meta.data[c(which(meta.data$V1 == "Capture Region") + 1):c(which(meta.data$V1 == "File Info Last Field Date  Time")-1),], error=function(e) e)
      if(inherits(test, "error")) 
      {
        
        meta.data[length(meta.data$V1) + 1,1] <- c("File Info Last Field Date  Time")
        meta.data[length(meta.data$V1),2] <- 11644444801
        
      }
      
      
      
      # Get timecodes for all frames
      g <- meta.data
      names(g) <- c("L1", "value")
      test <- tryCatch(
        g$value <- as.numeric(g$value), error=function(e){})
      
      
      g$L1 <- paste0(g$L1,"    ")
      g$Field <- unlist(lapply(g$L1, function(x) substring(x,gregexpr(' ',x)[[1]][1]+1,gregexpr(' ',x)[[1]][2]-1)))
      g$Type <- unlist(lapply(g$L1, function(x) substring(x,gregexpr(' ',x)[[1]][[2]],gregexpr(' ',x)[[1]][[3]])))
      g$Type2 <- unlist(lapply(g$Type, function(x)gsub("[[:space:]]", "", x)))
      
      cxd_fields_timecodes <- g[which(g$Type2 == "Time_From_Start"),]
      cxd_fields_timecodes <- cxd_fields_timecodes[order(cxd_fields_timecodes$Field),]
      cxd_fields_intervals <- g[which(g$Type2 == "Time_From_Last"),]
      cxd_fields_intervals <- cxd_fields_intervals[order(cxd_fields_intervals$Field),]
      
      # Get time interval
      total_time <- cxd_fields_timecodes[which.max(cxd_fields_timecodes$Field),2]
      
      cxd_meta.data <- NULL
      cxd_meta.data$sizeX <- as.numeric(meta.data$V2[which(meta.data$V1 == "Width ")])
      cxd_meta.data$sizeY <- as.numeric(meta.data$V2[which(meta.data$V1 == "Height ")])
      cxd_meta.data$sizeZ <- as.numeric(meta.data$V2[which(meta.data$V1 == "SizeZ ")])
      cxd_meta.data$sizeC <- 1
      cxd_meta.data$sizeT <- as.numeric(meta.data$V2[which(meta.data$V1 == "SizeT ")])
      
      
      # Check if files are very short (less than 200 frames) - skip to next file
      
      if(cxd_meta.data$sizeT < 200)
      {
        next
      }  
      
      cxd_meta.data$pixelType <- gsub(" ", "", meta.data$V2[which(meta.data$V1 == "Pixel type ")])
      cxd_meta.data$bitsPerPixel <- as.numeric(meta.data$V2[which(meta.data$V1 == "Valid bits per pixel ")])
      cxd_meta.data$imageCount <- as.numeric(cxd_fields_timecodes[which.max(cxd_fields_timecodes$Field),3])
      cxd_meta.data$dimensionOrder <- gsub(" ", "", meta.data$V2[which(meta.data$V1 == "Dimension order ")])
      cxd_meta.data$orderCertain <- gsub(" ", "", meta.data$V2[which(meta.data$V1 == "Dimension order ")])
      cxd_meta.data$rgb <- gsub(" ", "", meta.data$V2[which(meta.data$V1 == "RGB ")])
      cxd_meta.data$littleEndian <- gsub(" ", "", meta.data$V2[which(meta.data$V1 == "Endianness ")])
      cxd_meta.data$interleaved <- gsub(" ", "", meta.data$V2[which(meta.data$V1 == "Interleaved ")])
      cxd_meta.data$falseColor <- 0
      cxd_meta.data$metadataComplete <- 1
      cxd_meta.data$thumbnail <- 0
      cxd_meta.data$series <- 1
      cxd_meta.data$resolutionLevel <- 1
      cxd_meta.data$time_interval <- mean(as.numeric(cxd_fields_intervals$value), na.rm = TRUE)
      
      # Read Scale factor and binning from file - adjust if missing
      scale_factor <- as.numeric(str_split_fixed(meta.data$V2[which(meta.data$V1 == "factor")], ";", 2)[1])
      magnification <- as.numeric(str_split_fixed(meta.data$V2[which(meta.data$V1 == "magnification")], ";", 2)[1])
      
      if(is.na(magnification))
      {
        magnification <- 1  # if no magnification is given
      }
      
      
      if(scale_factor == 1 | is.null(scale_factor))
      {
        scale_factor <- 0.65 # if not calibrated and default to 1pxl - change to 10x setting for 6.5um CMOS chip.
      }  
      
      if(is.na(scale_factor))
      {
        stop(paste0("No calibration data found for file ", filelist$cxd[f],". Please move this file into a different folder and re-run script again without this file."))
      } 
      
      if(magnification == 2 & scale_factor == 0.65)   # To catch potentially wrongly set calibration in HCImage
      {  
        meta.data$coreMetadata$resolution <- 0.52
      }
      
      
      cxd_meta.data$resolution  <- scale_factor * magnification
      
      # CXD files from HCImage have a funky date of birth (Jan 1, 1601, 8:00am). We have to substract these seconds to get the Unix epoch
      cxd_meta.data$created_unix_from_file <- c(as.numeric(meta.data$V2[which(meta.data$V1 == "File Info Last Field Date  Time")]) - 11644444800) 
      
      
      write.csv(melt(cxd_meta.data), file = paste0(substr(filelist$cxd[f],1,gregexpr('.cxd',filelist$cxd[f])[[1]][1]),"_new_meta_data.csv"), row.names=FALSE)
      
      meta.data <- melt(cxd_meta.data)
      
    }
    
    
    
    X_ <- meta.data$value[1]
    Y_ <- meta.data$value[2]
    T_ <- meta.data$value[5]
    timepoints <- 500
    
    if (T_ < timepoints)
    {
      timepoints <- T_
    }
    
    # Adjustment for high-FPS movies - consider more timepoints
    if(meta.data[which(meta.data$L1 == "time_interval"),1] > 0.010)
    {
      next
    } else
    {
      timepoints <- T_
    }
    
    
    # Use maximum number of pixels x/y and default 
    # img2 <- read.image(filelist$cxd[f], normalize = TRUE, subset = list(x = 1:X_, y = 1:Y_, t=1:timepoints))
    img2 <- read.image(filelist$cxd[f], read.metadata = FALSE, normalize = FALSE)
    Xpos <- filelist$Xpos[f]
    image_ <- imageData(img2)
    rm(img2)
    # Normalize image 0 to 1
    gc()
    image_ <- image_ - min(image_)
    image_ <- image_ / max(image_)
    
    # standard deviation of each pixel over time to indicate which change the most over time
    y <- matrix(image_, X_ * Y_, T_)
    y1 <- rowSds(y)
    x <- matrix(y1, X_, Y_)
    
    # x <- apply(image_, c(1,2), sd)
    x <- x / max(x)     # Normalize
    
    
    # Reset all dark values to a minimum
    x[which(x < quantile(x)[4])] <- 0 
    
    # Intensity profile along X-axis, then find the ones above threshold
    v <- apply(x, 1, sum)
    
    above_thresh <- which(v >= mean(v))
    below_thresh <- which(v < mean(v))
    
    # What are the largest gaps in the above_tresh indices
    brights_ <- rle(diff(c(1,above_thresh)))
    bright_borders <- data.frame(brights_$values)
    bright_borders$lengths <- brights_$lengths
    bright_borders$steps <- bright_borders$brights_.values * bright_borders$lengths
    bright_borders$edge <- cumsum(bright_borders$steps)
    
    # below_thresh identifies the dark stripes - find the edges of them:
    stripes_ <- rle(diff(c(1,below_thresh)))
    stripe_borders <- data.frame(stripes_$values)
    stripe_borders$lengths <- stripes_$lengths
    stripe_borders$steps <- stripe_borders$stripes_.values * stripe_borders$lengths
    stripe_borders$edge <- cumsum(stripe_borders$steps)
    
    # Dark /  bright matrix
    border_mask <- data.frame(Xpos = c(1:dim(image_)[1]))  
    border_mask$type[border_mask$Xpos %in% above_thresh] <- "bright"
    border_mask$type[border_mask$Xpos %in% below_thresh] <- "dark"
    border_mask$edge[border_mask$Xpos %in% above_thresh] <- 100
    border_mask$edge[border_mask$Xpos %in% below_thresh] <- 1
    
    # Duplicate the last n (30) rows to fill
    
    border_mask[dim(image_)[1]:c(dim(image_)[1]+29),] <- border_mask[dim(image_)[1],]
    borders_ <- rollapply(border_mask$edge, 30, median)
    edges <- which(borders_ == 50.5)
    
    if (length(unique(border_mask$type[1:edges[1]])) != 1)
    {
      border_mask$type[1:edges[1]] <- border_mask$type[edges[1]]
      border_mask$edge[1:edges[1]] <- border_mask$edge[edges[1]]
    }
    
    if (border_mask$type[1] == "dark")
    {
      borders_ <- edges
    } else
    {
      borders_ <- c(1,edges)
    }
    
    if(length(borders_) < 2)
    {
      next
    }
    
    # Exclude a narrow (60pxl) bright anterior stripe
    
    if (is.na(diff(borders_[1:2])))
    {
      next
    }
    
    while (diff(borders_[1:2]) < 60)
    {
      borders_ <- borders_[3:length(borders_)]
    }
    
    # Compute the intensities along X for each bordered area
    bright_area <- apply(x, 1, trapz)
    bright_area1 <- bright_area[borders_[1]:borders_[2]]
    
    if(length(borders_) <= 3)
    {
      
      bright_area2 <- 0
      
    } else {
      
      bright_area2 <- bright_area[borders_[3]:borders_[4]]
      
      
    }
    
    # Reset minimum and find top 4 peaks in the 2 anterior stripes
    bright_area1[which(bright_area1 < quantile(bright_area1)[4])] <- quantile(bright_area1)[4]
    bright_area2[which(bright_area2 < quantile(bright_area2)[4])] <- quantile(bright_area2)[4]
    
    
    peaks1 <- unique(find_peaks(bright_area1))
    peaks2 <- unique(find_peaks(bright_area2))
    
    # Check if any peaks in peaks1 - if not, peaks2 becomes peaks1; reset borders as well
    if (length(peaks1) == 0)
    {
      peaks1 <- peaks2
      peaks2 <- peaks2[-c(1:length(peaks2))]
      bright_area1 <- bright_area2
      
      
      borders_[1] <- borders_[3]
      borders_[2] <- borders_[4]
      
      
    }
    
    # If no peaks are found - skip to the next file
    if (length(peaks1) == 0)
    {
      
      next  
      
    }
    
    peaklist <- list()
    peaklist[1] <- peaks1[which.max(bright_area1[peaks1])] + borders_[1]
    
    if(length(peaks1) > 1)
    { 
      peaklist[2] <- peaks1[which.maxN(bright_area1[peaks1], N=2)] + borders_[1]
    }
    
    if (length(peaks1) > 2)
    {
      peaklist[3] <- peaks1[which.maxN(bright_area1[peaks1], N=3)] + borders_[1]
    }
    
    
    
    if(length(peaks2) != 0)
    {
      peaklist[4] <- peaks2[which.max(bright_area2[peaks2])] + borders_[3]
      
      if(length(peaks2) > 1)
      {
        peaklist[5] <- peaks2[which.maxN(bright_area2[peaks2], N=2)] + borders_[3]
      }
      if (length(peaks2) > 2)
      {
        
        peaklist[6] <- peaks2[which.maxN(bright_area2[peaks2], N=3)] + borders_[3]
        
      }
      
    }
    
    peaklist <- unlist(peaklist)
    
    
    # Check if any peak is set outside the X-range - if so, set to max.X
    if(any(peaklist > dim(image_)[1], na.rm = TRUE))
    {
      peaklist[which(peaklist > dim(image_)[1])] <- dim(image_)[1]
    }
    
    # Highlight the Kymograph position in the image with a white line
    
    q <- x
    q[peaklist,] <- 1
    writeImage(q, file=paste0(substr(filelist$cxd[f],1,gregexpr('.cxd',filelist$cxd[f])[[1]][1]),'_SD and peaklines.tiff'))  
    
    peaklist <- peaklist[which(!is.na(peaklist))]
    
    if (length(peaklist) == 0)
    {
      next
    }
    
    # Create kymographs for peaks
    # Make sure to have unique peaks only !
    peaklist <- unique(peaklist)
    
    kymograph <- matrix(nrow = Y_, ncol = timepoints)
    graphs <- array(kymograph,dim = c(dim(image_)[2],dim(image_)[3],length(peaklist)))
    
    for (i in 1:length(peaklist)) # i-th peak-position {
    {
      for (j in 1:dim(image_)[3]) # all timepoints
      {
        kymograph[,j] <- image_[peaklist[i],,j]  
      }	
      graphs[,,i] <- kymograph
    }
    
    
    for (k in 1:length(peaklist))
    {
      
      writeImage(t(graphs[,,k]), file=paste0(filelist$cxd[f],'_peak_',k,'_at Xpos_', peaklist[k],'.tiff'))  
    }
    
    ## Add directionality 
    
    # Peaks
    k <- peaklist 
    k <- as.numeric(levels(k))[k]
    k <- c(k, seq(from= floor(dim(image_)[1] / 24), to=floor(dim(image_)[1] / 2), by=floor(dim(image_)[1] / 24))) # add 6 regularly spaced positions in addition
    k <- sort(k)
    k[which(k < 21)] <- 21
    k[which(k > (dim(image_)[1]) - 21)] <- dim(image_)[1] - 21
    k <- unique(k)
    
    
    # n pixel-wide Stripe for each peak  
    
    
    stripe_pic <- function(x) {
      n=20
      
      stripe_coordinates <- c(x-n):c(x+n)
      extent <- floor(dim(image_)[2] * 0.35)
      apply(image_[stripe_coordinates, c(1:extent, c(dim(image_)[2] - extent) : dim(image_)[2]) ,], c(3), mean)
    }
    
    Xpos_averages <- llply(k, stripe_pic)
    
    # Apply low-pass filter
    filter_peaks <- function(x){
      a <- x
      b <- 1:length(x)
      lowpass.spline <- smooth.spline(b,a, spar = 0.6)
      lowpass.loess <- loess(a ~ b, data = data.frame(x = a, y = b))
      highpass <- a - predict(lowpass.loess, b)
      
    }
    
    x <- llply(Xpos_averages, function (x) rollapply(x, 10, mean))
    
    # Normalize
    norma <- function(data){
      data <- data - min(data)
      data <- data / max(data)
    }
    
    x <- llply(x, norma)
    x <- llply(x, function (x) (x - 1) * -1)
    
    # Identify peaks
    findPeaks <- function (x, m = 10){
      shape <- diff(sign(diff(x, na.pad = FALSE)))
      pks <- sapply(which(shape < 0), FUN = function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
      })
      pks <- unlist(pks)
      pks
    }
    
    
    peaks_position_list <- llply(x, findPeaks)
    names(peaks_position_list) <- k
    
    # Aggregate all peaks into a single vector
    peaks_positions <- sort(unlist(peaks_position_list))
    
    # Define a window that collects all peaks within. Identify all windows
    # that show missing peaks.
    
    window_size <- floor(sd(diff(peaks_positions)))
    
    peak_groups <- list()
    for (i in 1:c(length(x[[1]]) - window_size))
    {
      peak_groups[[i]] <- sum(peaks_positions %in% i:c(i + window_size), na.rm = TRUE)
    }
    
    
    counted <- unlist(peak_groups)
    
    
    slide_and_count <- function(x) {
      
      peak_max <- NULL
      n = 1
      # Expected input is a vector with counts that has multiple series of 0s separating
      # the counts.
      
      for (i in 2:c(length(x) - 1))
      {
        if(x[i] == 0 && x[i + 1] != 0){
          n = i + 2
        }
        if(i == c(length(x) - 1) && x[i] != 0)
        {
          peak_max <- c(peak_max, max(x[n:i]))
        }
        
        if(x[i] == 0 && x[i - 1] != 0){
          peak_max <- c(peak_max, max(x[n:i-1]))
          next
        } 
        
      }    
      return(peak_max)   # returns the highest number for each interval - this is the expected peak number
    }     
    
    capture_peaks <- length(slide_and_count(counted))
    peaks_per_kymograph <- llply(peaks_position_list, length)
    
    duplicates_registry <- data.frame(sort(unlist(peaks_per_kymograph)))
    duplicates_registry$Xpos <- row.names(duplicates_registry)
    duplicates_registry <- split(duplicates_registry, duplicates_registry$sort.unlist.peaks_per_kymograph..)
    
    length_duplicates <- unlist(llply(duplicates_registry, function(x) length(unlist(x))[1]))/2
    
    if(any(length_duplicates < 2))
    {    
      duplicates_registry <- duplicates_registry[-which(length_duplicates < 2)]
    }  
    
    # Find all sets of peaks with identical number, and remove those that have no pairs
    
    to_use <- llply(duplicates_registry, function(x) names(peaks_position_list) %in% as.character(x[,2]))
    
    grouped_peaks <- llply(to_use, function (x) peaks_position_list[x])
    delta_list <- list()
    xpos_first <- list()
    xpos_last <- list()
    first_names <- list()
    last_names <- list()
    
    
    for (z in 1: length(grouped_peaks))
    {
      delta_list[[z]] <-  grouped_peaks[[z]][[1]] - grouped_peaks[[z]][[length(grouped_peaks[[z]])]]
      xpos_first[[z]] <-  grouped_peaks[[z]][[1]]
      first_names[[z]] <- duplicates_registry[[z]]$Xpos[1]
      
      xpos_last[[z]] <- grouped_peaks[[z]][[length(grouped_peaks[[z]])]]
      last_names[[z]] <- duplicates_registry[[z]]$Xpos[length(duplicates_registry[[z]]$Xpos)]
      
    }
    names(xpos_first) <- first_names
    names(xpos_last) <- last_names
    
    # Try to get as many kymographs as possible to maximize chances of good deltas
    
    all_coordinates_1 <- data.frame(cbind(xpos_first[[1]], xpos_last[[1]]))
    names(all_coordinates_1) <- c(first_names[[1]], last_names[[1]])
    
    X_diff <- diff(as.numeric(names(all_coordinates_1)))
    all_coordinates <- all_coordinates_1
    
    if(length(xpos_first) > 1)
    {
      all_coordinates_2 <- data.frame(cbind(xpos_first[[2]], xpos_last[[2]]))
      names(all_coordinates_2) <- c(first_names[[2]], last_names[[2]])
      X_diff_2 <- diff(as.numeric(names(all_coordinates_2)))
      if (X_diff_2 > X_diff){
        X_diff <- X_diff_2
        all_coordinates <- all_coordinates_2
      }
    }
    
    
    all_coordinates$delta <- all_coordinates[,1] - all_coordinates[,2]
    
    midpoint <- mean(unique(all_coordinates$delta))
    all_coordinates$direction <- "n.d."
    all_coordinates$direction[which(all_coordinates$delta >= midpoint)] <- "anterograde"
    all_coordinates$direction[which(all_coordinates$delta < midpoint)] <- "retrograde"
    
    write.csv(all_coordinates, file=paste0(filelist$cxd[f],'_directionmarks','.csv'), row.names = FALSE)  
    
  }
  
### Libraries necessary for the next batch of scripts
  
  
  library(broom)
  library(doParallel)
  library(plyr)
  library(dplyr)
  library(e1071)
  library(fifer)
  library(foreach)
  library(ggplot2)
  library(parallel)
  library(purrr)
  library(quantmod)
  library(signal)
  library(tidyr)
  library(tidyverse)
  library(tools)
  library(TTR)
  library(xlsx)
  library(baseline)
  
  
  
  # Move to target directory and create key folders
  setwd(normalizePath(target_dir))
  
  dir.create('TIFFs', showWarnings = FALSE, recursive = FALSE, mode = "0777")
  dir.create('balled', showWarnings = FALSE, recursive = FALSE, mode = "0777")
  
  # Copy all kymograph and metadata files from the movie directory to the target directories
  flist <- list.files(movie_dir, "csv$", full.names = TRUE, recursive = TRUE)
  if (length(flist) == 0)
  {
    stop("No CSV files found. Please select a different folder")
  }
  file.copy(flist, "balled")
  
  flist <- list.files(movie_dir, "*Xpos*", full.names = TRUE, recursive = TRUE)
  
  # Remove any Xpos files that are CSV or JPG
  flist <- flist[-grep(".csv", flist, fixed = TRUE)]
  flist <- flist[-grep(".jpg", flist, fixed = TRUE)]
  
  
  if (length(flist) == 0)
  {
    stop("No TIFF files found. Please select a different folder")
  }
  file.copy(flist, "TIFFs")
  
  
  # Apply rolling ball algorithm via Fiji/ImageJ macro (OS-dependent) on kymographs
  if(OS == "Darwin")
  {
    system(paste0("ImageJ-macosx --ij2 --headless --console --run \"", IJscript, "\" \'input=\"", normalizePath(target_dir), "/TIFFs/\",output=\"", normalizePath(target_dir), "/balled/\"\' "))  
  } else if(OS == "windows")
  {
    system(paste0("ImageJ-win64.exe --ij2 --headless --run \"", IJscript, "\"  \"input='", normalizePath(target_dir),"/TIFFs/', output='", normalizePath(target_dir),"/balled/'\""))
  } else if(OS == "unix")
  {
    system(paste0("ImageJ-linux64 --ij2 --headless --console --run \"", IJscript, "\" \'input=\"", normalizePath(target_dir), "/TIFFs/\",output=\"", normalizePath(target_dir), "/balled/\"\' "))  
  }

  
  
  # Tracing script - all kymographs become traced    
  setwd(normalizePath(file.path(target_dir, "balled")))
  # MAYO Screen Script No.2
  # Imports TIFF Kymographs (post ImageJ/FIJI batch 'Background substraction')
  # find . -iname *.tiff -exec /bin/cp "{}" /Volumes/HEART_DATA/ \; (for MacOSX)
  # find . -iname '*.tiff' -exec cp -n {} /mnt/MAYO\ data/TIFFs/ \; (for Linux)
  # find . -iname '*peaklines*' -exec mv {} peakline_tiffs/ \; 
  # and traces these. Exports JPG files with tracing overlay, and CSV table of
  # transients
  
  
  
  
  
  filelist <- NULL
  filelist$file <- list.files(".", pattern = "\\.tiff$")
  
  # Populate the list with Genotypes, derived from the filename
  filelist$Genotype <- factor(substr(filelist$file,1,regexpr("_",filelist$file[])-1))
  
  # Set averaging value for smoothing step
  avgx <- 4
  
  total <- length(filelist$file)
  pb <- txtProgressBar(max = total, style = 3)
  
  for (f in 1:length(filelist$file))
    
  {
    setTxtProgressBar(pb, f)
    
    if(file.exists(paste0(filelist$file[f],"_traced.jpg")))
    {
      next
    } 
    
    img2 <- read.image(filelist$file[f])
    
    if (length(dim(img2)) > 2) # Make sure we only use images that are x,y (not RGB)
      next(f)
    
    image_ <- imageData(img2)
    image_ <- t(image_)
    
    image_ <- image_ / max(image_) 
    image_ <- gblur(image_, sigma = 1.5) # improves edge trnaemsacing
    
    
    # remove bright edges with sine half-wave
    p <- seq(0,pi,pi/c(length(image_[,1])-1))
    sin_p <- sin(p)
    
    modifier_2 <- replicate(length(image_[1,]), sin_p)
    
    # Adjust image with this modifier
    image_ <- image_ * modifier_2
    
    # Rolling mean for each time slice
    image_mean <- NULL
    dd <- function(x) {unlist(rollmean(x, avgx))}
    image_mean <- apply(image_,1,dd)
    
    # Normalize everything from to 0..1
    image_mean <- sweep(image_mean,1,apply(image_mean , 1, max),`/`)  
    
    # Calculate the value range - substract the bottom quantile
    background_ <- quantile(image_mean)[1]
    image_mean_bg <- image_mean - background_
    
    
    # Transpose image
    image_mean_bg <- t(image_mean_bg)
    
    
    # Edge tracing
    
    up <- list()
    down <- list()
    
    for (i in 1:length(image_mean_bg[1,])) # first run to get the outlines
    {
      x <- as.vector(mean(image_mean_bg[,i]))
      test1 <- which(image_mean_bg[,i] > x)
      up_ <- test1[1] # upper (left) boundary
      
      test2 <- which(rev(image_mean_bg[,i]) > x)
      down_left <- test2[1] # lower (right) boundary
      down_ <- c(length(image_mean_bg[,1])) - down_left
      up[[i]] <- unlist(up_)	
      down[[i]] <- unlist(down_)
      
    }
    
    up <- unlist(up)
    down <- unlist(down)
    positions <- 1:length(image_mean_bg[1,])
    
    
    ## Sanity check - make sure that the values are not too different between sequential x-positions
    limits <- mean(up) - sd(up)
    
    to_be_corrected <- which(up <= limits)
    g = 1
    
    while(length(to_be_corrected) > 1)
    {
      to_be_corrected <- which(up <= limits)
      g = g + 1
      for (i in to_be_corrected)
      {
        
        x <- as.vector(mean(image_mean_bg[,i], na.rm=TRUE))
        
        test1 <- which(image_mean_bg[,i] > x)
        
        up_ <- test1[g] # second boundary
        
        up[i] <- unlist(up_)
        
      }
    }
    
    limits_2 <- mean(down) + sd(down)
    
    to_be_corrected <- which(down >= limits_2)
    g = 1
    
    while(length(to_be_corrected) > 1)
    {
      to_be_corrected <- which(down >= limits_2)
      g = g + 1
      for (i in to_be_corrected)
      {
        
        x <- as.vector(mean(image_mean_bg[,i], na.rm=TRUE))
        test2 <- which(rev(image_mean_bg[,i]) > x)
        down_left <- test2[g] # second boundary
        down_ <- c(length(image_mean_bg[,1])) - down_left
        down[i] <- unlist(down_)
        
      }
    }
    
    
    diameters <- data.frame(cbind(up, down))
    diameters$distance <- diameters[,2] - diameters[,1]
    
    # New data table that uses temporal and spatial resolution data to calculate real distances
    meta.data <- read.csv(paste0(substr(filelist$file[f],1,gregexpr('.cxd',filelist$file[f])[[1]][1]),"_new_meta_data.csv"))
    
    time_interval <- meta.data[which(meta.data$L1 == "time_interval"),1]
    resolution_meta <- meta.data[which(meta.data$L1 == "resolution"),1]
    
    diameters$time <- c(0, time_interval * 1:c(length(diameters$distance)-1))
    
    for_export <- diameters[,c(4,3,1,2)]
    for_export[,2] <- for_export[,2] * resolution_meta 
    
    
    # Define upper half of the movie
    up_midline <- as.integer(max(up) - sd(up))
    Half_up_image <- image_mean_bg[1:up_midline,,drop = F]
    
    down_midline <- as.integer(max(down) - sd(down))
    Half_down_image <- image_mean_bg[length(image_[,1]):down_midline,,drop = F]
    
    
    # Brightness changes as edge moves into the frame of the halfmovie
    auc_function <- function(x) {unlist(trapz(x))}
    up_brightness <- apply(Half_up_image,2,auc_function)
    down_brightness <- apply(Half_down_image,2,auc_function)
    
    for_export$up_bright <- up_brightness
    for_export$down_bright <- down_brightness
    
    for_export$down_norm <- c(for_export$down - min(for_export$down)) / c(max(for_export$down) - min(for_export$down))[-length(for_export$down)]
    for_export$up_norm <- c(for_export$up - min(for_export$up)) / c(max(for_export$up) - min(for_export$up))[-length(for_export$up)]
    for_export$distance_norm <- for_export$down_norm - for_export$up_norm
    
    write.table(for_export, file = paste0(filelist$file[f],".csv"), row.names = FALSE, col.names = TRUE, sep = ",")
    
    
    image_diameters <- matrix(nrow=length(image_mean_bg[,1]) , ncol = length(image_mean_bg[1,]))
    image_diameters[] <- 0 
    image_diameters[cbind(up, positions)] <- 255
    image_diameters[cbind(down, positions)] <- 255
    
    # x = image_diameters
    # y = image_mean
    # Since the image would turn upside-down we need to re-arrange the rows to be correct
    x = t(image_diameters) #[nrow(image_diameters):1,, drop = FALSE]
    y = image_mean
    
    x <- x[,ncol(x):1, drop = FALSE]
    y <- y[,ncol(y):1, drop = FALSE]
    
    
    red_green = rgbImage(green=x, red=y, blue = y)
    #red_green <- transpose(red_green)
    writeImage(red_green, file=paste0(filelist$file[f],"_traced.jpg"), type = "JPEG")
  }
  close(pb)
  
  
  # Tracing quality control step - copies files into a 'good' and 'bad traces' folder
  # MAYO Screen Script No.3
  # Imports transients stored in CSV files of traced Kymographs 
  # and performs quality control. Data will be split into '
  # Good traces' and 'Bad traces'. Exports summary statistics as well.
  
  
  
  # Get file list, for all necessary file types
  
  filelist <- NULL
  filelist$csv <- list.files(".", pattern = "\\.tiff.csv$")
  filelist$smoothie_up <- 1:length(filelist$csv)
  filelist$smoothie_up_norm <- filelist$smoothie_up
  filelist$smoothie_down <- filelist$smoothie_up
  filelist$smoothie_down_norm <- filelist$smoothie_up
  filelist$sd <- filelist$smoothie_up
  filelist$dist_cor <- filelist$smoothie_up
  filelist$norm_dist_cor <- filelist$smoothie_up
  
  total <- length(filelist$csv)
  pb <- txtProgressBar(max = total, style = 3)
  
  errorlist <- NULL
  # https://stats.stackexchange.com/questions/24607/how-to-measure-smoothness-of-a-time-series-in-r
  # Standard deviation of the time-to-time differences, should be minimal for continuous data
  
  for (f in 1:length(filelist$csv))
  {
    setTxtProgressBar(pb, f)
    x <- read.csv(filelist$csv[f])
    
    if (length(colnames(x)) < 4)
    {
      next
    }
    
    
    #filelist$smoothie_up_norm[f] <-  sd(diff(x$up))/abs(mean(diff(x$up), na.rm=TRUE))
    
    if(sd(diff(x$down)) != 0 & sd(diff(x$up)) != 0)
    {
      filelist$smoothie_up[f] <-  sd(diff(x$up))
      filelist$smoothie_down[f] <-  sd(diff(x$down))
      filelist$autocorr_up[f] <- cor(diff(x$up)[-length(diff(x$up))],diff(x$up)[-1])
      filelist$autocorr_down[f] <- cor(diff(x$down)[-length(diff(x$down))],diff(x$down)[-1])
      filelist$autocorr_sd[f] <- sd(c(filelist$autocorr_up[f], filelist$autocorr_down[f]))
      filelist$sd[f] <- sd(c(filelist$smoothie_down[f], filelist$smoothie_up[f]))
      filelist$dist_cor[f] <- cor(diff(x$distance)[-length(diff(x$distance))],diff(x$distance)[-1])
      filelist$norm_dist_cor[f] <- cor(diff(x$distance_norm)[-length(diff(x$distance_norm))],diff(x$distance_norm)[-1])
      
      
    } else {
      filelist$smoothie_up[f] <-  10
      filelist$smoothie_down[f] <- 10
      #filelist$smoothie_down_norm[f] <-  sd(diff(x$down))/abs(mean(diff(x$down), na.rm=TRUE))
      filelist$autocorr_up[f] <- 0
      filelist$autocorr_down[f] <- 0
      filelist$autocorr_sd[f] <- 10
      filelist$sd[f] <- 100
      filelist$dist_cor[f] <- 0
      filelist$norm_dist_cor[f] <- 0
      
    }
    
    
    
    
    
    
  }
  
  close(pb)
  
  
  # Add genotype info to allow for sorting by genotype
  x <- data.frame(filelist)
  x$CODE <- as.factor(unlist(lapply(filelist$csv, function(x) substring(x,1,gregexpr('_',x)[[1]]-1)[1])))
  x$ID <- unlist(lapply(filelist$csv, function(x) substring(x,gregexpr('_',x)[[1]][1]+1,gregexpr('_',x)[[1]][2]-1)))
  x$Xpos <- as.numeric(unlist(lapply(filelist$csv, function(x) substring(x,gregexpr('Xpos_',x)[[1]]+5,gregexpr('.tiff',x)[[1]]-1)[1])))
  x$quality <- 'n.t.'
  x$quality[which(x$autocorr_up > median(x$autocorr_up, na.rm = TRUE) * 0.9 & x$autocorr_down > median(x$autocorr_down, na.rm = TRUE) * 0.9 & x$sd < median(x$sd, na.rm = TRUE) * 1.5)] <- 'excellent'
  x$quality[which(x$quality == 'n.t.' & x$dist_cor > 0.6 & x$sd < 4 & x$norm_dist_cor > 0.6)] <- 'good'
  x$tiff <- unlist(lapply(filelist$csv, function(x) substring(x,1,gregexpr('.tiff',x)[[1]]+4)[1]))
  x$jpeg <- unlist(lapply(filelist$csv, function(x) paste0(substring(x,1,gregexpr('.tiff',x)[[1]]+4)[1], "_traced.jpg")))
  
# Separate data into 3 folders - excellent, 'good and 'bad' traces. Copy csv files of excellent traces as well.
  
  mainDir <- getwd()
  subDir <- paste0("bad traces", collapse =' ')
  subDirPath <- file.path(mainDir, subDir)
  dir.create(subDirPath, showWarnings = FALSE)
  file.copy(as.character(x$csv[which(x$quality == "n.t.")]),subDirPath, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date=FALSE)
  file.copy(x$jpeg[which(x$quality == "n.t.")],subDirPath, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date=FALSE)
  
  subDir2 <- paste0("good traces", collapse =' ')
  subDirPath2 <- file.path(mainDir, subDir2)
  dir.create(subDirPath2, showWarnings = FALSE)
  file.copy(as.character(x$csv[which(x$quality == "good")]),subDirPath2,overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date=FALSE)
  file.copy(x$jpeg[which(x$quality == "good")],subDirPath2, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date=FALSE)
  
  subDir3 <- paste0("excellent traces", collapse =' ')
  subDirPath3 <- file.path(mainDir, subDir3)
  dir.create(subDirPath3, showWarnings = FALSE)
  file.copy(as.character(x$csv[which(x$quality == "excellent")]),subDirPath3, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date=FALSE)
  file.copy(x$jpeg[which(x$quality == "excellent")],subDirPath3, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date=FALSE)
  
  write.table(x, file = "Quality_control.csv", row.names = FALSE, col.names = TRUE, sep = ",")
  
  
# Load QC data
  # QC <- read.csv(file = "Quality_control.csv")
  QC <- x
  # Identify for each file if at least 1 file is of excellent quality
  QC$cxd <- unlist(lapply(QC$jpeg, function(x) substring(x, 1, regexpr("cxd", x) + 2)))
  
  QC_files <- QC[,c(15,18)]
  
  QC_files <- melt(table(QC_files))
  
  QC_files <- split(QC_files, QC_files$quality)
  
  # No excellent traces:
  need_more_traces <- unique(QC_files$excellent$cxd[which(QC_files$excellent$value == 0)])
  
  # Select these movies from good traces:
  rescue_traces <- QC_files$good[QC_files$good$cxd[need_more_traces],]
  rescue_traces <- rescue_traces[which(rescue_traces$value > 0), ]
  
  QC_good <- QC[which(QC$quality == "good"),]
  rescue_traces_QC <- merge(rescue_traces, QC_good, by.x = "cxd", by.y = "cxd", all.x = TRUE)
  
  file.copy(as.character(rescue_traces_QC$csv),subDirPath3, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date=FALSE)
  file.copy(rescue_traces_QC$jpeg,subDirPath3, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date=FALSE)
  
  
  
  # Final run - several scripts that count beats and aggregate the data  
  setwd(normalizePath(file.path(target_dir, "balled", "excellent traces")))
  # setwd("/home/geo/GIANT/balled_MAYO_only/excellent traces")
  file.copy(normalizePath(mappings_file), ".")
  
  
  # Transient analysis
  # MAYO Screen Script No.4
  # Imports transients stored in CSV files of good Kymographs 
  # and identifies intervals.
  # A file that maps the Codes with genotypes should be provided
  # names 'mappings.csv' with the following structure:
  #
  # CODE	    cross	                                type	      fly	    human
  # MAYO0001	Hand4.2 tdtK x BL-41698 - CG4747 RNAi	experiment	CG4747	GLYR1
  # ...
  
  # This file will export .Rdata files that contain all raw data and will
  # be used subsequently for downstream analysis, i.e. graphs, statistics...
  
  # Adjust findpeaks parameters for Peaks and Valleys to better match acutal peaks before rollmean
  
  findPeaksP <- function (x, thresh = 0) 
  {
    pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 
                   0) + 1
    if (!missing(thresh)) {
      if (sign(thresh) < 0) 
        thresh <- -thresh
      pks[x[pks - 1] - coredata(x[pks]) > thresh]
    }
    else pks
  }
  
  findPeaksV <- function (x, thresh = 0) 
  {
    pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 
                   0) + 1
    if (!missing(thresh)) {
      if (sign(thresh) < 0) 
        thresh <- -thresh
      pks[x[pks - 1] - coredata(x[pks]) > thresh]
    }
    else pks
  }
  options(warn=2)
  
  filelist <- NULL
  filelist$jpg <- list.files(".", pattern = "\\.jpg$")
  filelist$file <- gsub("_traced.jpg", ".csv", filelist$jpg)
  filelist$Genotype <- as.factor(unlist(lapply(filelist$jpg, function(x) substring(x,1,gregexpr('_',x)[[1]]-1)[1])))
  
  # Extract more information from the filenames
  # Genotype is:
  filelist$Genotype_filename <- as.factor(unlist(lapply(filelist$file, function(x) substring(x,1,gregexpr('_',x)[[1]]-1)[1])))
  
  # Age is:
  filelist$Age_filename <- as.factor(unlist(lapply(filelist$file, function(x) substring(x,gregexpr('_',x)[[1]][2]+1,gregexpr('_',x)[[1]]+1)[2])))
  
  # Sex is:
  filelist$Sex_filename <- as.factor(unlist(lapply(filelist$file, function(x) substring(x,gregexpr('_',x)[[1]][3]-1,gregexpr('_',x)[[1]]-1)[3])))
  
  # Comment is:
  filelist$comment <- unlist(lapply(filelist$file, function(x) substring(x,gregexpr('_',x)[[1]][5]+1, gregexpr('\\.',x)[[1]]-1)[1]), use.names = TRUE)
  filelist$comment[which(is.na(filelist$comment))] <- "none"

  # Populate the list with Genotypes, derived from the filename
  filelist$Genotype <- factor(substr(filelist$file,1,regexpr("_",filelist$file[])-1))
  
  # Export the Genotype list (can be used to create Instructions.csv file)
    if(mappings_file_type == "CSV")
  {
    library('openxlsx')
    write.xlsx(read.csv("mappings.csv"), "mappings.xlsx")
  }
  
  crosses <- read.xlsx2(file="mappings.xlsx", sheetIndex = 1)
  crosses$CODE <- gsub("_", "", crosses$CODE)
  
  z <- crosses[,c(1,2,3,6,7)]
  names(z) <- c("CODE",  "cross", "type",  "fly",   "human")
  
  # Is sorting and/or grouping present? If not, group by symbols and type
  if (length(names(z)) == 5)
  {
    symbols_ <- list(as.factor(levels(z$human)))
    symbols_[[2]] <- 1:length(levels(z$human))
    z$sorting <- symbols_[[2]][match(z$human, symbols_[[1]])]
    
    groups_ <- list(levels(z$type))
    groups_[[2]] <- 1:length(levels(z$type))
    z$group <- groups_[[2]][match(z$type, groups_[[1]])]
  }
  z$sorting <- as.factor(z$sorting)
  z$group	<- as.factor(z$group)
  z_original <- z
  
  
  filelist$Genotype <- factor(substr(filelist$file,1,regexpr("_",filelist$file[])-1))
  filelist$index <- c(1:length(filelist$file))

  
  # Extract more information from the filenames
  # Genotype is:
  filelist$Genotype_filename <- as.factor(unlist(lapply(filelist$file, function(x) substring(x,1,gregexpr('_',x)[[1]]-1)[1])))
  #write.table(filelist$Genotype_filename, file = "Genotype_filename.csv", sep = ",", row.names = FALSE)
  
  # Age is:
  filelist$Age_filename <- as.factor(unlist(lapply(filelist$file, function(x) substring(x,gregexpr('_',x)[[1]][2]+1,gregexpr('_',x)[[1]]+1)[2])))
  
  # Sex is:
  filelist$Sex_filename <- as.factor(unlist(lapply(filelist$file, function(x) substring(x,gregexpr('_',x)[[1]][3]-1,gregexpr('_',x)[[1]]-1)[3])))
  
  
  # Comment is:
  filelist$comment <- unlist(lapply(filelist$file, function(x) substring(x,gregexpr('_',x)[[1]][5]+1, gregexpr('\\.',x)[[1]]-1)[1]), use.names = TRUE)
  filelist$comment[which(is.na(filelist$comment))] <- "none"
  
  # Reconstruct CODE/Genotype to include age and sex (adjust z)
  filelist$CODED <- factor(paste0(filelist$Genotype,filelist$Age_filename, filelist$Sex_filename))
  filelist$CODE <- factor(filelist$Genotype)
  filelist$cross <- z$cross[match(filelist$CODE,z$CODE)]
  filelist$type <- z$type[match(filelist$CODE,z$CODE)]
  filelist$symbol_fly <- z$fly[match(filelist$CODE,z$CODE)]
  filelist$symbol_human <- z$human[match(filelist$CODE,z$CODE)]
  filelist$sorting <- z$sorting[match(filelist$CODE,z$CODE)]
  filelist$group <- z$group[match(filelist$CODE,z$CODE)]
  filelist$Genotype <- factor(filelist$CODED) # New Genotype value
  
  z_new <- NULL
  z_new$CODE <- filelist$Genotype
  z_new$cross <- filelist$cross[match(z_new$CODE, filelist$Genotype)]
  z_new$type <- filelist$type[match(z_new$CODE, filelist$Genotype)]
  z_new$fly <- filelist$symbol_fly[match(z_new$CODE, filelist$Genotype)]
  z_new$human <- filelist$symbol_human[match(z_new$CODE, filelist$Genotype)]
  z_new$origCode <- filelist$CODE[match(z_new$CODE, filelist$Genotype)]
  z_new$sorting <- filelist$sorting[match(z_new$CODE, filelist$Genotype)]
  z_new$group <- filelist$group[match(z_new$CODE, filelist$Genotype)]
  z_new <- data.frame(z_new)
  z <- z_new
  
  save(z, file = "z.Rdata")
  save(filelist, file="filelist.Rdata")
  
print("Codes assigned")
  
  
  # Create master table with overview of all files and indices
  genotype_indices <- data.frame(aggregate(filelist$index, by= list(filelist$Genotype), function (x) x))
  
  names(genotype_indices) <- c("Code", "Index")
  genotype_indices$n <- sapply(genotype_indices$Index, function (x) length(x))
  save(genotype_indices, file="genotype_indices.Rdata")
  
  filelist_with_index <- data.frame(filelist)
  filelist_with_index$IDnumber <- as.numeric(unlist(lapply(filelist_with_index$file, function(x) as.numeric(substring(x,gregexpr('_',x)[[1]][1]+1,gregexpr('_',x)[[1]]-1)[2]))))
  filelist_with_index$flyID <- as.factor(sapply(1:length(filelist_with_index$file), function(i) substring(filelist_with_index$file[i], gregexpr('_',filelist_with_index$file[i])[[1]] +1 , gregexpr('w[fm]_',filelist_with_index$file[i])[[1]]+1)[[1]]))
  filelist_with_index$flyID <- as.factor(sapply(1:length(filelist_with_index$file), function(i) paste0(filelist_with_index$Genotype_filename[i], "_", filelist_with_index$flyID[i])))
  save(filelist_with_index, file="filelist_with_index.Rdata")
  
  
  direction_files <- data.frame("direction" = list.files(path = "../", pattern = "_directionmarks.csv"))
  direction_files$CXD <- unlist(llply(direction_files[,1], function(x) substr(x, 1, gregexpr('.cxd', x)[[1]][1]+3)))
  direction_jpg_matched_files <- data.frame("traced_csv" = filelist_with_index$file)
  direction_jpg_matched_files$CXD <- unlist(llply(direction_jpg_matched_files[,1], function(x) substr(x, 1, gregexpr('.cxd', x)[[1]][1]+3)))
  directions <- merge(direction_files, direction_jpg_matched_files, by.x = "CXD", by.y = "CXD", all.y = TRUE )
  
print("Filelist done")

# Work through the transients now
# By genotype (Code, i), then by fly (j)

avgx = 5 	# Steps to average in rolling mean for transient
avgt = 400   # Steps for rollmean in threshold transient


# Process each file and gather spatiotemporal data

number_of_codes <- length(genotype_indices$Code)

# To save memory we will split the codes into chunks
# and store these in a temporary sub-folder
split_codes <- split(1:number_of_codes, ceiling(seq_along(1:number_of_codes)/20))
mainDir <- getwd()
subDir <- paste0("temp", collapse =' ')
subDirPath <- file.path(mainDir, subDir)
dir.create(subDirPath, showWarnings = FALSE)

for(p in 1:length(split_codes))
{
  message("Calculating Genotypes in chunks of 20, part ", p)
  intervals_final <- NULL
  
  # Internal chunk counter:
  cc <- 0
  
  # Get current sets
  current_code_range <- split_codes[[p]]
  
  # Initialize List to store data
  transients_experiment <- vector("list", length(current_code_range))
  timescale <- vector("list", length(current_code_range))
  total_time <- vector("list", length(current_code_range))
  intervals_final <- vector("list", length(current_code_range))
  transients_final <- vector("list", length(current_code_range))
  identicals_final <- vector("list", length(current_code_range))
  time_used_final <- vector("list", length(current_code_range))
  transient_metrics_final <- vector("list", length(current_code_range))
  spline_all <- vector("list", length(current_code_range))
  EDD_this_genotype <- vector("list", length(current_code_range))
  ESD_this_genotype <- vector("list", length(current_code_range))
  FS_this_genotype <- vector("list", length(current_code_range))
  meanEDD <- vector("list", length(current_code_range))
  meanESD <- vector("list", length(current_code_range))
  directionality_this_genotype <- vector("list", length(current_code_range))
  
  
  
  for (i in current_code_range[1]:current_code_range[length(current_code_range)])
  {
    
    cc <- cc + 1
    
    # Set and reset variables to store loop data mylist <- vector("list", length(genotype_indices$Index[[i]]))
    interval_this_genotype <- vector("list", length(genotype_indices$Index[[i]]))
    transients_this_genotype <- vector("list", length(genotype_indices$Index[[i]]))
    identicals_this_genotype <- vector("list", length(genotype_indices$Index[[i]]))
    time_used <- vector("list", length(genotype_indices$Index[[i]]))
    transient_metrics_this_genotype <- vector("list", length(genotype_indices$Index[[i]]))
    EDD_this_fly <- vector("list", length(genotype_indices$Index[[i]]))
    ESD_this_fly <- vector("list", length(genotype_indices$Index[[i]]))
    FS_this_fly <- vector("list", length(genotype_indices$Index[[i]]))
    mean_EDD_this_fly <- vector("list", length(genotype_indices$Index[[i]]))
    mean_ESD_this_fly <- vector("list", length(genotype_indices$Index[[i]]))
    directionality_this_fly <- vector("list", length(genotype_indices$Index[[i]]))
    spline_this_genotype <- vector("list", length(genotype_indices$Index[[i]]))
    
    print(paste0(i, " of ",number_of_codes, " done (", round(c(i/number_of_codes * 100), digits = 2), "%)"))
    for (j in 1:(genotype_indices$n[i]))
      
    {
      v <- genotype_indices$Index[[i]][j]
      transients_experiment <- read.csv(file = as.character(filelist_with_index$file[match(v, filelist_with_index$index)]), header=TRUE, sep = ",")[,1:2]
      names(transients_experiment) <- c("time", "distance")
      timescale[[v]] <- median(diff(transients_experiment[,1]))
      total_time[[v]] <- length(transients_experiment[,1]) * timescale[[v]]
      
      directions_this_file <- read.csv(file = paste0("../" , as.character(directions$direction[which(directions$traced_csv %in% as.character(filelist_with_index$file[match(v, filelist_with_index$index)]))])), header=TRUE, sep = ",", stringsAsFactors = FALSE)
      beat_directions <- melt(table(directions_this_file$direction))
      total_directions <- sum(beat_directions$value)
      anterograde_beats <- beat_directions$value[which(beat_directions$Var1 == "anterograde")]
      directionality_this_fly[[j]] <- anterograde_beats / total_directions
      
      
      x2 <- transients_experiment
      x2[2] <- x2[2] * -1 # Current transient, inverted
      x2[2] <- round(x2[2], digits = 2)
      
      
      # Baseline correction
      corr_x2 <- data.frame(distance = c(baseline(rbind(x2[,2]), wm=20, ws=20, method = 'rollingBall')@corrected))
      corr_x2$time <- x2[,1]
      corr_x2 <- corr_x2[,2:1]
      
      # Generate smooth spline for original data (for best peaks)    
      
      dd <- data.frame(corr_x2)														
      
      k <- merge(x2, dd[,1:2], by.x='time', by.y = 'time', all.y = TRUE)
      
      # Identify peaks and valleys
      
      valleys <- findPeaksV(-dd$distance)
      hills <- findPeaksP(dd$distance)
      
      # Check if Valley/Peak indices are assigned beyond length of transient
      if(last(hills) > length(dd$time))
      {
        hills[length(hills)] <-  length(dd$time)
      }
      
      if(last(valleys) > length(dd$time))
      {
        valleys[length(valleys)] <-  length(dd$time)
      }
      
      dd$peak <- NA
      dd$peak[hills] <- "Peak"
      dd$peak[valleys] <- "Valley"
      
      k$peak <- NA
      
      # Correct for supershort 
      if(length(dd$time) < avgt)
      {
        avgt <- as.integer(length(dd$time) / 3)
      }
      
      # define a threshold distance that separates valleys from peaks
      threshold_dd <- data.frame(unlist(rollmeanr(dd[2], avgt)))
      
      # Reset avgt if it was changed for shor movies
      avgt <- 400
      
      # Add median distance to the trailing values of the threshold
      extension_ <- c(length(dd$time) - length(threshold_dd[,1])) # positions to fill
      threshold_dd[length(threshold_dd[,1]):c(length(threshold_dd[,1])+extension_),] <- median(tail(threshold_dd[,1])) # fill with median
      
      threshold_dd$time <- dd[1:c(length(dd[,1])),1]
      threshold_dd <- threshold_dd[,c(1,2)]
      threshold_dd <- threshold_dd[,2:1] # Order columns	
      
      # Peak and Valley IDs
      V_id <- which(dd$peak == "Valley")
      P_id <- which(dd$peak == "Peak")
      
      # Find above-threshold valley and remove it
      
      h <- merge(dd[V_id,],  threshold_dd, by.x = "time", by.y="time")
      t <- h[which(h[,2] >= h[,4]),1]
      dd$peak[dd$time %in% t] <- NA
      
      # Find below-threshold hill and remove it
      
      h <- merge(dd[P_id,],  threshold_dd, by.x = "time", by.y="time")
      t <- h[which(h[,2] <= h[,4]),1]
      dd$peak[dd$time %in% t] <- NA
      
      # Determine peak/valley pairs (find a valley, followed by a peak, followed by a valley)
      # Keep only those as real valley/peaks
      
      all_pos <- which(!is.na(dd$peak))
      all_pos_kinds <- dd$peak[all_pos]
      
      peaky <- function(x)
      {
        return(all(x == c("Peak",  "Peak")))
      }
      
      
      if(length(which(rollapply(all_pos_kinds, width = 2, by = 1, FUN = peaky, align = "left"))) == 0)
      {
        all_pos_cleaned <- all_pos
      } else
      {
        all_pos_cleaned <- all_pos[-c(which(rollapply(all_pos_kinds, width = 2, by = 1, FUN = peaky, align = "left")) + 1)]
      }
      
      
      new_peaklist <- dd$peak[all_pos_cleaned]
      
      
      correct_peak_pattern <- function(x)
      {
        return(all(x == c("Valley",  "Peak",  "Valley")))
      }
      
      start_positions <- all_pos_cleaned[c(which(rollapply(new_peaklist, width = 3, by = 1, FUN = correct_peak_pattern, align = "left")))]
      peak_postions <-  all_pos_cleaned[c(which(rollapply(new_peaklist, width = 3, by = 1, FUN = correct_peak_pattern, align = "left")) + 1)]
      tail_positions <-  all_pos_cleaned[c(which(rollapply(new_peaklist, width = 3, by = 1, FUN = correct_peak_pattern, align = "left")) + 2)]
      
      position_index <- data.frame(starts = start_positions)
      position_index$peaks <- peak_postions
      position_index$ends <- tail_positions
      
      time_index <- data.frame(matrix(x2$time[as.matrix(position_index)],nrow = length(position_index$ends),ncol = 3))
      names(time_index) <- names(position_index)
      
      # quartz()
      # plot(x2, cex = 0.5, type = "l")
      # abline(v=x2$time[position_index$starts], col = 'red')
      # abline(v=x2$time[position_index$peaks], col = 'blue')
      # abline(v=x2$time[position_index$ends], col = 'green')
      
      # Create working copy
      transients_all <- k[,c(1,2,4)]
      names(transients_all) <- names(dd)
      
      
      # Diameter-pairs
      EndDiameters <- data.frame(DD = transients_all[start_positions,2])
      EndDiameters$SD <- transients_all[peak_postions,2]
      EndDiameters$FS <- c(EndDiameters$DD - EndDiameters$SD) / EndDiameters$DD
      
      EDD_this_fly[[j]] <- median(EndDiameters$DD)
      ESD_this_fly[[j]] <- median(EndDiameters$SD)
      FS_this_fly[[j]] <-  median(EndDiameters$FS)
      
      mean_EDD_this_fly <- mean(EDD_this_fly[[j]] * -1)
      mean_ESD_this_fly <- mean(ESD_this_fly[[j]] * -1)
      
      
      # Set beginning with first valley and end with valley
      transients_all <- transients_all[start_positions[1]:length(transients_all[,1]),]
      transients_all <- transients_all[1:last(tail_positions),]
      row.names(transients_all) <- NULL 
      
      
      # Get intervals lengths for rhythmicity:
      intervals <- time_index
      
      intervals$SI <- intervals$ends - intervals$starts
      intervals$deltaSI <- c(diff(intervals$SI), NA)
      intervals$deltaSIplus <- NA
      intervals$deltaSIplus[1:c(length(intervals$SI)-1)] <- intervals$deltaSI[2:length(intervals$SI)] 
      intervals$HP <- c(diff(intervals$starts), NA)
      intervals$DI <- intervals$HP - intervals$SI
      intervals$i <- i
      intervals$j <- j
      intervals$Index <- v
      interval_this_genotype[[j]] <- intervals
      
      intervals_final[[cc]] <- dplyr::bind_rows(interval_this_genotype)
      
      # Extract all transients
      transients <- list()
      interval_mark <- list()
      
      for (n in 1:length(position_index$starts))
      {
        
        transients[[n]] <- x2[position_index$starts[n]:position_index$ends[n],]
        transients[[n]]$diameter <- transients[[n]]$distance * -1
        transients[[n]]$timestamp <- transients[[n]]$time
        transients[[n]]$time <- transients[[n]]$time - transients[[n]]$time[1] # set t to t0
        transients[[n]]$distance <- transients[[n]]$distance - transients[[n]]$distance[1] # set distance start to 0
        transients[[n]]$distance <- transients[[n]]$distance / max(transients[[n]]$distance) # normalize distance
        transients[[n]]$delta <- c(transients[[n]]$diameter - transients[[n]]$diameter[1]) * -1
        transients[[n]]$from_EDD <- mean_EDD_this_fly - transients[[n]]$delta
        transients[[n]]$velocity <- c(0,diff(transients[[n]]$delta) / diff(transients[[n]]$time))
        
        # Modify start index to skip leading periods of non-contraction
        transients[[n]] <- transients[[n]][c(which(transients[[n]]$velocity != 0)[1] - 1):length(transients[[n]]$velocity),]
        transients[[n]]$time <- transients[[n]]$time - transients[[n]]$time[1]
        transients[[n]]$i <- i
        transients[[n]]$j <- j
        transients[[n]]$index <- v
        transients[[n]]$beat <- n
      }
      
      transients_this_genotype[[j]] <- bind_rows(transients)
      transients_final[[cc]] <- bind_rows(transients_this_genotype)
      
      
      AUC <- list() # Area under the curve - to identify super small transients
      for (m in 1:(length(transients)))
      {
        AUC[[m]] <- trapz(transients[[m]][,2])
      }
      AUC <- unlist(AUC)
      
      # Select all transients (place to select transients of certain size)
      upper <- quantile(AUC, na.rm = TRUE)[5] 
      lower <- quantile(AUC, na.rm = TRUE)[1]
      
      identicals <- transients[which(AUC >= lower & AUC <= upper)]
      
      # Remove very short transients
      transient_lengths <- sapply(identicals, nrow)
      
      identicals <- identicals[which(transient_lengths > 10)]
      
      identicals_this_genotype[[j]] <- identicals
      names(identicals_this_genotype) <- rep(genotype_indices$Code[i], length(identicals_this_genotype))
      
      identicals_final[[cc]] <- identicals_this_genotype
      
      
      # Length of each transient
      used_time <- list()
      
      for (n in 1:(length(identicals)))
      {
        used_time[n] <- max(identicals[[n]][,1])
      }
      
      time_used[[j]] <- sum(unlist(used_time))
      time_used_final[[cc]] <- unlist(time_used)
      
      
      
      # Analysis
      transient_metrics <- list()
      velocities <- list()
      for (q in 1:length(identicals))
      {
        # Create smooth spline
        transient_metrics$spl[[q]] <- smooth.spline(identicals[[q]][,1:2])
        
        # Determine the timepoints with max positive and negative slopes
        velocities$max_velocity[[q]] <- identicals[[q]]$velocity[which.max(identicals[[q]]$velocity)]
        velocities$max_neg_velocity[[q]]  <- identicals[[q]]$velocity[which.min(identicals[[q]]$velocity)]
        
      }
      
      transient_metrics_this_genotype[[j]] <- velocities
      transient_metrics_this_genotype[[j]]$genotype <- genotype_indices$Code[i]
      transient_metrics_final[[cc]] <- transient_metrics_this_genotype
      
      spline_this_genotype[[j]] <- transient_metrics$spl
      spline_all[[cc]] <- spline_this_genotype
      
      
      EDD_this_genotype[[cc]] <- EDD_this_fly
      ESD_this_genotype[[cc]] <- ESD_this_fly
      FS_this_genotype[[cc]] <-  unlist(FS_this_fly)
      mean_EDD_this_fly <- mean(EDD_this_fly[[j]] * -1)
      mean_ESD_this_fly <- mean(ESD_this_fly[[j]] * -1)
      meanEDD[[cc]] <- mean_EDD_this_fly
      meanESD[[cc]] <- mean_ESD_this_fly
      directionality_this_genotype[[cc]] <- directionality_this_fly
      
      
    }
    
  }
  
  intervals_final <- bind_rows(intervals_final)
  
  # Store each set individually
  fileName_intervals_final <- paste0("temp/Genotype_chunks_intervals_final_part_", p,".Rds")
  saveRDS(intervals_final, file=fileName_intervals_final)
  
  fileName_timescale <- paste0("temp/Genotype_chunks_timescale_part_", p,".Rds")
  saveRDS(timescale, file=fileName_timescale)
  
  fileName_total_time <- paste0("temp/Genotype_chunks_total_time_part_", p,".Rds")
  saveRDS(total_time, file=fileName_total_time)
  
  fileName_transients_final <- paste0("temp/Genotype_chunks_transients_final_part_", p,".Rds")
  saveRDS(total_time, file=fileName_transients_final)
  
  fileName_identicals_final <- paste0("temp/Genotype_chunks_identicals_final_part_", p,".Rds")
  saveRDS(identicals_final, file=fileName_identicals_final)
  
  fileName_time_used_final <- paste0("temp/Genotype_chunks_time_used_final_part_", p,".Rds")
  saveRDS(time_used_final, file=fileName_time_used_final)
  
  fileName_transient_metrics_final <- paste0("temp/Genotype_chunks_transient_metrics_final_part_", p,".Rds")
  saveRDS(transient_metrics_final, file=fileName_transient_metrics_final)
  
  fileName_spline_all <- paste0("temp/Genotype_chunks_spline_all_part_", p,".Rds")
  saveRDS(spline_all, file=fileName_spline_all)
  
  fileName_EDD_this_genotype <- paste0("temp/Genotype_chunks_EDD_this_genotype_part_", p,".Rds")
  saveRDS(EDD_this_genotype, file=fileName_EDD_this_genotype)
  
  fileName_ESD_this_genotype <- paste0("temp/Genotype_chunks_ESD_this_genotype_part_", p,".Rds")
  saveRDS(ESD_this_genotype, file=fileName_ESD_this_genotype)
  
  fileName_FS_this_genotype <- paste0("temp/Genotype_chunks_FS_this_genotype_part_", p,".Rds")
  saveRDS(FS_this_genotype, file=fileName_FS_this_genotype)
  
  fileName_meanEDD <- paste0("temp/Genotype_chunks_meanEDD_part_", p,".Rds")
  saveRDS(meanEDD, file=fileName_meanEDD)
  
  fileName_meanESD <- paste0("temp/Genotype_chunks_meanESD_part_", p,".Rds")
  saveRDS(meanESD, file=fileName_meanESD)
  
  fileName_directionality_this_genotype <- paste0("temp/Genotype_chunks_directionality_this_genotype_part_", p,".Rds")
  saveRDS(directionality_this_genotype, file=fileName_directionality_this_genotype)
  
}  

# Combine Spline_all chunks to compute more temporal data:
{
  
  spline_all <- list()
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_spline_all_part_", p,".Rds")
    data_read <- readRDS(fileName)
    spline_all <- c(spline_all, data_read)
  } 
  
}

# Store complete data of splines
{
  save(spline_all, file="spline_all.Rdata")
}

# Transient analysis pt.1
metrics <- function(y){
  
  max.value <- which.max(y$y)
  
  # Contraction range
  Contraction_range <- y$y[1:max.value]
  
  # Relaxation range
  Relaxation_range <- y$y[c(max.value + 1): length(y$y)]
  
  # Time to 10%/25%... peak
  tt10p <- y$x[which(Contraction_range > 0.1)[1]]
  tt25p <- y$x[which(Contraction_range > 0.25)[1]]
  tt50p <- y$x[which(Contraction_range > 0.5)[1]]
  tt75p <- y$x[which(Contraction_range > 0.75)[1]]
  tt90p <- y$x[which(Contraction_range > 0.9)[1]]
  
  # Time to 15%/25%... return
  tt90r <- y$x[which(Relaxation_range < 0.9)[1] + max.value]
  tt75r <- y$x[which(Relaxation_range < 0.75)[1] + max.value]
  tt50r <- y$x[which(Relaxation_range < 0.5)[1] + max.value]
  tt25r <- y$x[which(Relaxation_range < 0.25)[1] + max.value]
  tt10r <- y$x[which(Relaxation_range < 0.10)[1] + max.value]
  
  list(y, flatten(list(peaked = y$x[max.value], tt10p = tt10p, tt25p = tt25p, tt50p = tt50p, tt75p = tt75p, 
                       tt90p = tt90p, tt90r = tt90r, tt75r = tt75r, tt50r = tt50r, tt25r = tt25r, tt10r = tt10r)))
  
}

# Split spline numerical data evenly into chunks of 20 and save individually
split_spline_indices <- split(1:length(spline_all), ceiling(seq_along(1:length(spline_all))/20))
for(i in 1:length(split_spline_indices))
{
  metrics_all <- list()
  message("Calculating Splines part ", i)
  
  y <- spline_all[split_spline_indices[[i]]]
  
  for (k in 1:length(y))
  {
    u <- list()  
    values <- list()
    for (j in 1:length(y[[k]]))
    {
      u[[j]] <- lapply(y[[k]][[j]], function(x) predict(x, x=seq(0,max(x$x), by=0.0002)))
      
      values[[j]] <- lapply(u[[j]], function(x) metrics(x))
      values[[j]]$name <- genotype_indices$Code[split_spline_indices[[i]][k]]
      values[[j]]$id <- j
      
    }
    metrics_all[[k]] <- values  
    
  }
  
  fileName <- paste0("temp/Splines_part_", i,".Rds")
  saveRDS(metrics_all, file=fileName)  
  
}  

# Combine chunks back into single object
metrics_read <- list()
complete <- list()
for(i in 1:length(split_spline_indices))
{
  
  fileName <- paste0("temp/Splines_part_", i,".Rds")
  metrics_read <- readRDS(fileName)
  complete <- c(complete, metrics_read)
} 

save(metrics_all, file = "metrics_all_splined.Rdata")


# Load remaining data chunks:
# Combine chunks back into single object
{
  timescale <- list()
  total_time <- list()
  intervals_final <- list()
  transients_final <- list()
  identicals_final <- list()
  time_used_final <- list()
  transient_metrics_final <- list()
  EDD_this_genotype <- list()
  ESD_this_genotype <- list()
  FS_this_genotype <- list()
  meanEDD <- list()
  meanESD <- list()
  directionality_this_genotype <- list()
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_intervals_final_part_", p,".Rds")
    data_read <- readRDS(fileName)
    intervals_final <- c(intervals_final, data_read)
  } 
  intervals_final <- as.data.frame(intervals_final)
  
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_FS_this_genotype_part_", p,".Rds")
    data_read <- readRDS(fileName)
    FS_this_genotype <- c(FS_this_genotype, data_read)
  } 
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_identicals_final_part_", p,".Rds")
    data_read <- readRDS(fileName)
    identicals_final <- c(identicals_final, data_read)
  } 
  
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_time_used_final_part_", p,".Rds")
    data_read <- readRDS(fileName)
    time_used_final <- c(time_used_final, data_read)
  } 
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_EDD_this_genotype_part_", p,".Rds")
    data_read <- readRDS(fileName)
    EDD_this_genotype <- c(EDD_this_genotype, data_read)
  } 
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_transient_metrics_final_part_", p,".Rds")
    data_read <- readRDS(fileName)
    transient_metrics_final <- c(transient_metrics_final, data_read)
  } 
  
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_timescale_part_", p,".Rds")
    data_read <- readRDS(fileName)
    timescale <- c(timescale, data_read)
  } 
  
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_meanEDD_part_", p,".Rds")
    data_read <- readRDS(fileName)
    meanEDD <- c(meanEDD, data_read)
  } 
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_ESD_this_genotype_part_", p,".Rds")
    data_read <- readRDS(fileName)
    ESD_this_genotype <- c(ESD_this_genotype, data_read)
  } 
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_total_time_part_", p,".Rds")
    data_read <- readRDS(fileName)
    total_time <- c(total_time, data_read)
  } 
  
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_transients_final_part_", p,".Rds")
    data_read <- readRDS(fileName)
    transients_final <- c(transients_final, data_read)
  } 
    
  for(p in 1:length(split_codes))
  {
    fileName <-  paste0("temp/Genotype_chunks_directionality_this_genotype_part_", p,".Rds")
    data_read <- readRDS(fileName)
    directionality_this_genotype <- c(directionality_this_genotype, data_read)
  }   
  
  
  
}
# Store all data
{
  save(timescale, file="timescale.Rdata")
  save(total_time, file="total_time.Rdata")
  save(intervals_final, file="intervals_final.Rdata")
  save(transients_final, file="transients_final.Rdata")
  save(identicals_final, file="identicals_final.Rdata")
  save(time_used_final, file="time_used_final.Rdata")
  save(transient_metrics_final, file="transient_metrics_final.Rdata")
  save(EDD_this_genotype, file="EDD_this_genotype.Rdata")
  save(ESD_this_genotype, file="ESD_this_genotype.Rdata")
  save(FS_this_genotype, file="FS_this_genotype.Rdata")
  save(meanEDD, file="meanEDD.Rdata")
  save(meanESD, file="meanESD.Rdata")
  save(directionality_this_genotype, file="directionality_this_genotype.Rdata")
  
}


# Transient analysis pt.2


# load(file="genotype_indices.Rdata")
# load(file = "metrics_all_splined.Rdata")

final_transients <- matrix(nrow = tail(cumsum(unlist(lapply(metrics_all, function(x) lengths(x) -2   ))), n=1), ncol = 14) # 11 columns for 11 measurements + 3 indices
final_transients <- data.frame(final_transients)

j=1
i=1


all_js <- c(0, cumsum(lengths(metrics_all)))

starts <- c(0, cumsum(unlist(lapply(metrics_all, function(x) lengths(x) -2)))) + 1
ends <- c(0, cumsum(unlist(lapply(metrics_all, function(x) lengths(x) -2))))


for (i in 1:length(metrics_all))
{
  a <- metrics_all[[i]]
  
  
  for (j in 1:length(a))
  {
    b <- a[[j]]
    
    index_ <- genotype_indices$Index[[i]][j]
    
    
    data <- unlist(lapply(lapply(b[1], '[', 2), '[[', 1))
    
    for (k in 2:c(length(b)-2))
    {
      data <- rbind(data,unlist(lapply(lapply(b[k], '[', 2), '[[', 1)))
    } 
    
    row.names(data) <- NULL
    colnames(final_transients) <- c("i", "j", "index", colnames(data))
    
    first <- starts[all_js[i] + j]
    last <-  ends[all_js[i] + j+1]
    
    final_transients[first:last,] <- data.frame(cbind(i = i, j = j, index = index_, data))
  }
  
}

save(final_transients, file = "final_transients.Rdata")


# Metadata and final aggregate
# Meta data extraction and aggregation

filelist <- NULL
filelist$csv <- list.files("../", pattern = "\\_new_meta_data.csv$", recursive = TRUE, full.names = TRUE)

total <- length(filelist$csv)


# Create empty matrix to be filled with metadata
baseline <- read.csv(filelist$csv[1], sep = ",")
meta_data <- matrix(nrow = length(filelist$csv), ncol = 21) 
meta_data <- data.frame(meta_data)
names(meta_data) <- baseline[,2]
meta_data$filename <- NA


pb <- txtProgressBar(max = total, style = 3)

for (f in 1:length(filelist$csv))
{
  setTxtProgressBar(pb, f)
  
  baseline <- read.csv(filelist$csv[f], sep = ",", stringsAsFactors = F)
  
  
  meta_data[f,] <- c(t(as.numeric(baseline[1:21,1])) , filelist$csv[f])
} 

meta_data$filename <- gsub("..//", "", meta_data$filename)
meta_data$CODE <- NA
meta_data$CODE <- factor(substr(meta_data$filename,1,regexpr("_",meta_data$filename[])-1))

meta_data$flyID <- NA
meta_data$flyID <- as.factor(unlist(lapply(meta_data$filename, function(x) substring(x,1,gregexpr('_',x)[[1]][3]-1))))

meta_data$recordingday <- as.Date(as.POSIXct(as.numeric(meta_data$created_unix_from_file), origin = "1970-01-01"))
meta_data$week_year <- format(meta_data$recordingday, format="%Y-%U")

# Fix dates that are '1' - use the date from the creation date 
# MAYO data/Screen Data$      find . -iname '*.cxd' -exec stat -f "%m%t%Sm %N" "{}" > unix_data.csv \;

# In R, load this file into a data frame:
# - if present!
if(file.exists("test_unix_data.csv")){
  
  
  date_created_raw <- read.table("test_unix_data.csv", fill =  TRUE, as.is=TRUE, sep="\t")
  # First row is the unix timecode, so the time the file was modified the last time (i.e. created)
  # Problem is to get the filename from the second column
  # 
  # Solution: define a function that splits a string according to a character, then unlist this string
  # and determine its length (i.e. how many substrings = depth of the path).
  # Then take the last string from the split string list (which is the cxd.file), and replace the .cxd
  # ending with .xml to allow for comparison with xml filelist
  
  cxdfile <- function(x) {
    y <- length(unlist(strsplit(x, "/")))
    x <- unlist(strsplit(x, "/"))[y]
    #gsub("cxd", "xml", x)
  } 
  
  # Apply this function to all lines of the file list, then add the filename as extra column
  date_created_raw$file <- lapply(date_created_raw[,1], cxdfile)
  date_created_raw <- date_created_raw[!duplicated(date_created_raw$V1),]
  date_created_raw <- date_created_raw[!duplicated(date_created_raw$V1),]
  date_created_raw$flyID <- as.factor(unlist(lapply(date_created_raw$file, function(x) substring(x,1,gregexpr('_',x)[[1]][3]-1))))
  date_created_raw$V1 <- unlist(lapply(date_created_raw$V1, function(x) substring(x,1,gregexpr(' ',x)[[1]][1]-1)))
  date_created_raw$file <- unlist(date_created_raw$file )
  
  meta_data$cxd_date <- date_created_raw$V1[match(meta_data$flyID, date_created_raw$flyID)]
  meta_data$cxd_file <- date_created_raw$file[match(meta_data$flyID, date_created_raw$flyID)]
  
} else
{
  meta_data$cxd_date <- NA
  meta_data$cxd_file <- meta_data$filename
}


for (i in 1:length(meta_data$cxd_date))
{
  if(!is.na(meta_data$cxd_date[i]))
  {
    next
  } else
  {
    meta_data$cxd_date[i] <- meta_data$created_unix_from_file[i]
  }
}


meta_data$cxd_recordingday <- as.Date(as.POSIXct(as.integer(meta_data$cxd_date), origin = "1970-01-01"))
meta_data$cxd_week_year <- format(meta_data$cxd_recordingday, format="%U %Y")

# Fix filenames
meta_data$flyID <- gsub("..//", "", meta_data$flyID)
meta_data$filename <- gsub("..//", "", meta_data$filename)
meta_data$CODE <- gsub("..//", "", meta_data$CODE)

meta_data$cxd_file <- gsub("_new_meta_data.csv", "cxd", meta_data$cxd_file)
write.csv(meta_data, file="meta_data_all.csv", row.names = F)


# MAYO Screen Script No.5
# Imports .Rdata files and performs downstream analysis 
# to determine phenotypes and statistical differences.


# Define functions
row.compare <- function (x)	{
  x[order(x$ID, -x$Length),]
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Import Master cross list
crosses <- z

crosses <- crosses[,c(6,2,3,4,5)]
colnames(crosses) <- c("CODE", "cross", "type", "gene", "human.ortholog")


genotype_indices <- data.frame(aggregate(filelist$index, by= list(filelist$Genotype), function (x) x))
names(genotype_indices) <- c("Code", "Index")
genotype_indices$n <- sapply(genotype_indices$Index, function (x) length(x))
genotype_indices$cross <- z$cross[match(genotype_indices$Code, z$CODE)]
genotype_indices$human <- z$human[match(genotype_indices$Code, z$CODE)]


genotypes <- z$CODE

number_of_codes <- length(genotype_indices$Code)

# Data quality
coverage <- list()

for (i in 1:number_of_codes)
{
  for (j in 1:(genotype_indices$n[i]))
  {
    # Determine how many seconds of the whole transient have been used
    coverage[unlist(genotype_indices$Index[[i]][j])] <- unlist(time_used_final[[i]][[j]]) / unlist(total_time[unlist(genotype_indices$Index[[i]][j])])
  }}
coverage <- unlist(coverage)

# Compute coverage and beat statistics per genotype
master_table <- genotype_indices[,c(1,3)]
flies <- function (x) {genotype_indices$Index[[x]]	}


master_table$Cross <- z$cross[match(master_table$Code, z$CODE)]
master_table$type <- z$type[match(master_table$Code, z$CODE)]
master_table$human <- z$human[match(master_table$Code, z$CODE)]
master_table$fly <- z$fly[match(master_table$Code, z$CODE)]

master_table$AgeSex <- as.factor(unlist(lapply(as.character(master_table$Code), function (x) substrRight(x,2))))
master_table$Sex <- as.factor(unlist(lapply(as.character(master_table$Code), function (x) substrRight(x,1))))
master_table$Age <- as.factor(unlist(lapply(as.character(master_table$AgeSex), function (x) substr(x,1,1))))
master_table$coverage <- sapply(1:length(master_table$Code), function (x) mean(coverage[flies(x)]))
master_table$SDcoverage <- sapply(1:length(master_table$Code), function (x) sd(coverage[flies(x)]))



# Data to extract
# Heart Period = average transient length
# Heart Rate = beats per minute
# Analysis

HR_all <- list()
EDD <- list()
ESD <- list()
FS <- list()
HR_this_fly <- list()
SIs <- list()
DIs <- list()
FS_SD <- list()
files_ <- list()
min_velocity <- list()
max_velocity <- list()


for (i in 1:number_of_codes) # loop through genotypes
{
  
  HR_this_fly <- vector("list", genotype_indices$n[i])
  EDD_this_fly <- list()
  ESD_this_fly <- list()
  min_velocity_ <- list()
  max_velocity_ <- list()
  files_this_genotype <- list()
  
  for (j in 1:(genotype_indices$n[i])) # loop through flies
  {
    
    HR_this_fly[[j]] <- length(identicals_final[[i]][[j]]) / round(unlist(total_time[unlist(genotype_indices$Index[[i]][j])]), digits = 2)
    
    EDD_this_fly[[j]] <- mean(EDD_this_genotype[[i]][[j]], na.rm = TRUE) * -1
    ESD_this_fly[[j]] <- mean(ESD_this_genotype[[i]][[j]], na.rm = TRUE) * -1
    
    files_this_genotype[[j]] <- filelist_with_index[genotype_indices$Index[[i]][j],1]
    min_velocity_[[j]] <- mean(transient_metrics_final[[i]][[j]]$max_neg_velocity, na.rm = TRUE)
    max_velocity_[[j]] <- mean(transient_metrics_final[[i]][[j]]$max_velocity, na.rm = TRUE)
    
    
  }
  
  files_[[i]] <- unlist(files_this_genotype)
  
  EDD[[i]]	<- unlist(EDD_this_fly)
  names(EDD)[[i]] <- levels(genotype_indices$Code)[i]
  
  ESD[[i]]	<- unlist(ESD_this_fly)
  names(ESD)[[i]] <- levels(genotype_indices$Code)[i]
  
  FS[[i]] <- FS_this_genotype[[i]]
  names(FS)[[i]] <- levels(genotype_indices$Code)[i]
  
  FS_SD[[i]] <- sd(FS_this_genotype[[i]], na.rm = TRUE)
  names(FS_SD)[[i]] <- levels(genotype_indices$Code)[i]
  
  names(FS_this_genotype)[[i]] <- levels(genotype_indices$Code)[i]
  
  min_velocity[[i]] <- unlist(min_velocity_)
  names(min_velocity)[[i]] <- levels(genotype_indices$Code)[i]
  
  max_velocity[[i]] <- unlist(max_velocity_)
  names(max_velocity)[[i]] <- levels(genotype_indices$Code)[i]
  
  
  HR_all[[i]] <- unlist(HR_this_fly)
  names(HR_all)[[i]] <- levels(genotype_indices$Code)[i]
  
}

master_table$HR <- sapply(1:length(master_table$Code), function (x) mean(HR_all[[x]]))
master_table$SDHR <- sapply(1:length(master_table$Code), function (x) sd(HR_all[[x]]))

master_table$EDD <- sapply(1:length(master_table$Code), function (x) mean(EDD[[x]]))
master_table$SDEDD <- sapply(1:length(master_table$Code), function (x) sd(EDD[[x]]))

master_table$ESD <- sapply(1:length(master_table$Code), function (x) mean(ESD[[x]]))
master_table$SDESD <- sapply(1:length(master_table$Code), function (x) sd(ESD[[x]]))

master_table$FS <- sapply(1:length(master_table$Code), function (x) mean(FS[[x]]))
master_table$SDFS <- sapply(1:length(master_table$Code), function (x) sd(FS[[x]]))

master_table$max_velocity <- sapply(1:length(master_table$Code), function (x) mean(max_velocity[[x]]))
master_table$SDMV <- sapply(1:length(master_table$Code), function (x) sd(max_velocity[[x]]))

master_table$min_velocity <- sapply(1:length(master_table$Code), function (x) mean(min_velocity[[x]]))
master_table$SDMinV <- sapply(1:length(master_table$Code), function (x) sd(min_velocity[[x]]))

write.table(master_table, file = "master_table.csv", sep = ",", row.names = FALSE)



# Files
files_used <- reshape2::melt(files_)


# EDD
EDD_data <- reshape2::melt(EDD)
EDD_data$CODE <- as.factor(unlist(lapply(EDD_data$L1, function(x) substring(x,1,gregexpr('[mf]',x)[[1]]-2)[1])))
EDD_data$human <- crosses$human.ortholog[match(EDD_data$CODE, crosses$CODE)]
EDD_data$fly <- crosses$gene[match(EDD_data$CODE, crosses$CODE)]
EDD_data$type <- crosses$type[match(EDD_data$CODE, crosses$CODE)]
names(EDD_data)[1] <- "EDD"


# ESD
ESD_data <- reshape2::melt(ESD)
ESD_data$CODE <- as.factor(unlist(lapply(ESD_data$L1, function(x) substring(x,1,gregexpr('[mf]',x)[[1]]-2)[1])))
ESD_data$human <- crosses$human.ortholog[match(ESD_data$CODE, crosses$CODE)]
ESD_data$fly <- crosses$gene[match(ESD_data$CODE, crosses$CODE)]
ESD_data$type <- crosses$type[match(ESD_data$CODE, crosses$CODE)]
names(ESD_data)[1] <- "ESD"


# FS
FS_data <- reshape2::melt(FS)
FS_data$CODE <- as.factor(unlist(lapply(FS_data$L1, function(x) substring(x,1,gregexpr('[mf]',x)[[1]]-2)[1])))
FS_data$human <- crosses$human.ortholog[match(FS_data$CODE, crosses$CODE)]
FS_data$fly <- crosses$gene[match(FS_data$CODE, crosses$CODE)]
FS_data$type <- crosses$type[match(FS_data$CODE, crosses$CODE)]
names(FS_data)[1] <- "FS"


# HR
HR_data <- reshape2::melt(HR_all)
HR_data$CODE <- as.factor(unlist(lapply(HR_data$L1, function(x) substring(x,1,gregexpr('[fm]',x)[[1]]-2)[1])))
HR_data$human <- crosses$human.ortholog[match(HR_data$CODE, crosses$CODE)]
HR_data$fly <- crosses$gene[match(HR_data$CODE, crosses$CODE)]
HR_data$type <- crosses$type[match(HR_data$CODE, crosses$CODE)]
names(HR_data)[1] <- "HR"



# neg velocity
min_velocity_data <- reshape2::melt(min_velocity)
min_velocity_data$CODE <- as.factor(unlist(lapply(min_velocity_data$L1, function(x) substring(x,1,gregexpr('[mf]',x)[[1]]-2)[1])))
min_velocity_data$human <- crosses$human.ortholog[match(min_velocity_data$CODE, crosses$CODE)]
min_velocity_data$fly <- crosses$gene[match(min_velocity_data$CODE, crosses$CODE)]
min_velocity_data$type <- crosses$type[match(min_velocity_data$CODE, crosses$CODE)]
names(min_velocity_data)[1] <- "min_velocity"

# max_velocity
max_velocity_data <- reshape2::melt(max_velocity)
max_velocity_data$CODE <- as.factor(unlist(lapply(max_velocity_data$L1, function(x) substring(x,1,gregexpr('[fm]',x)[[1]]-2)[1])))
max_velocity_data$human <- crosses$human.ortholog[match(max_velocity_data$CODE, crosses$CODE)]
max_velocity_data$fly <- crosses$gene[match(max_velocity_data$CODE, crosses$CODE)]
max_velocity_data$type <- crosses$type[match(max_velocity_data$CODE, crosses$CODE)]
names(max_velocity_data)[1] <- "max.velocity"


write.table(EDD_data, file = "EDD_table.csv", sep = ",", row.names = FALSE)
write.table(ESD_data, file = "ESD_table.csv", sep = ",", row.names = FALSE)
write.table(FS_data, file = "FS_table.csv", sep = ",", row.names = FALSE)
write.table(HR_data, file = "HR_table.csv", sep = ",", row.names = FALSE)
write.table(min_velocity_data, file = "min_velocity_data.csv", sep = ",", row.names = FALSE)
write.table(max_velocity_data, file = "max_velocity_table.csv", sep = ",", row.names = FALSE)


all_data <- cbind(files_used[,1], EDD_data[,1], ESD_data[,1], FS_data[,1], HR_data[,1], min_velocity_data[,1], max_velocity_data)

names(all_data) <- c("file", "EDD", "ESD", "FS", "HR","min.velocity", "max.velocity", "CODE_long", "CODE", "human_gene", "fly_gene", "type")
all_data$cross <- crosses$cross[match(all_data$CODE, crosses$CODE)]

all_data$flyID <- as.factor(sapply(1:length(all_data$file), function(i) substring(all_data$file[i], gregexpr('_',all_data$file[i])[[1]] +1 , gregexpr('w[fm]_',all_data$file[i])[[1]]+1)[[1]]))
all_data$flyID <- as.factor(sapply(1:length(all_data$file), function(i) paste0(all_data$CODE[i], "_", all_data$flyID[i])))
all_data$ID_number <- as.factor(unlist(lapply(all_data$file, function(x) as.numeric(substring(x,gregexpr('_',x)[[1]][1]+1,gregexpr('_',x)[[1]]-1)[2]))))
all_data$Xpos <-  as.numeric(sapply(1:length(all_data$file), function(i) substring(all_data$file[i], gregexpr('Xpos_',all_data$file[i])[[1]] + 5 , gregexpr('.tiff',all_data$file[i])[[1]]-1)[[1]]))


meta_data <- read.csv(file="meta_data_all.csv", sep = ",")

all_data <- merge(all_data, meta_data, by.x = "flyID", by.y = "flyID", all.x = TRUE)


# MAD - based on final_transients calculate MAD_tt10r

final_transients_grouped <- group_by(final_transients, i, j, index)
MAD_final_transients <- dplyr::select(final_transients_grouped, 1:3, 14)
MAD_tt10r <- MAD_final_transients %>% summarise_all(c("mad"), na.rm = TRUE)

# Arrhythmia Index (AI)
intervals_grouped_AI <- group_by(intervals_final, i, j, Index)
intervals_grouped_AI <- dplyr::select(intervals_grouped_AI, 7, 9:11)
AI_index <- intervals_grouped_AI %>% summarise_all(c("sd", "median"), na.rm = TRUE)
AI_index$AI <- AI_index$sd / AI_index$median

# Combine MAD_tt10r and AI
merged_AI_MADtt10r <- left_join(AI_index, MAD_tt10r, by = c("i" = "i", "j" = "j", "Index" = "index"))
merged_AI_MADtt10r_files <- left_join(filelist_with_index, merged_AI_MADtt10r, by = c("index" = "Index"))
names(merged_AI_MADtt10r_files)[24] <- "MAD_tt10r"

# Stroke volume - based on all_data table result calculate SV
Volume_SV <- function(x,y){  
  pi / 4 * (x^2 - y^2)
}

all_data_grouped_SV <- group_by(all_data, flyID, file)
all_data_grouped_SV <- dplyr::select(all_data_grouped_SV, 1:4, 6, 16)
all_data_grouped_SV$SV <- Volume_SV(all_data_grouped_SV$EDD, all_data_grouped_SV$ESD)
all_data_grouped_SV$CO <- all_data_grouped_SV$SV * all_data_grouped_SV$HR

# Identify segment based on project

# Typical value is about 330px for a 1000pxl-wide movie.
threshold_Xpos <- floor(mean(unique(all_data_grouped_SV$Xpos)))

abdominal_segment <- function(x) {
  if(x > threshold_Xpos) {
    return("posterior")
  } else {
    return("anterior")
  }}

all_data_grouped_SV$segment <- sapply(all_data_grouped_SV$Xpos, abdominal_segment)
all_data_grouped_SV$direction_anterior <- unlist(directionality_this_genotype)
new_parameters <- left_join(merged_AI_MADtt10r_files, all_data_grouped_SV, by = c("jpg" = "file"))
new_parameters <- new_parameters[,c(1:4,19:20,23:24,26:32)]
new_parameters$jpg_noX <- new_parameters$jpg
new_parameters$anterograde_percent <- unlist(directionality_this_genotype)
new_parameters <- new_parameters[,-c(2,5,6,9:12)]

# Remove Peak-Number from filename
new_parameters$jpg_noX <- as.factor(unlist(lapply(new_parameters$jpg_noX, function(x) paste0(substring(x, 1, gregexpr('peak_',x)[[1]][1]+4), "X",substring(x,gregexpr('_at',x)[[1]][1])))))

# Merge all data with final tables

all_data_added <- left_join(all_data, new_parameters, by = c("file" = "jpg"))
write.table(all_data_added, file = "all_data_table.csv", sep = ",", row.names = FALSE)





# Transient analysis
meta_data <- as_tibble(read.csv(file="meta_data_all.csv", sep = ","))
crosses <- read.xlsx2(file="mappings.xlsx", sheetIndex = 1)
crosses$CODE <- gsub("_", "", crosses$CODE)
crosses <- crosses[,c(1,4,5,9)]


# Organize by genotype and flyID and index
grouped <- group_by(final_transients, i, j, index)
grouped_per_fly <- grouped %>% summarise_all(c("median","sd"), na.rm = TRUE)

# Add the actual genotype name via index
grouped_per_fly$CODE <- filelist_with_index$Genotype_filename[grouped_per_fly$index %in% filelist_with_index$index]
grouped_per_fly$jpg <- filelist_with_index$jpg[match(grouped_per_fly$index, filelist_with_index$index)]
grouped_per_fly <- grouped_per_fly %>% dplyr::select("CODE", "jpg", everything()) # Re-order columns

intervals_grouped <- group_by(intervals_final, i, j, Index)
#intervals_grouped_per_fly <- summarise_all(intervals_grouped, funs(median,sd), na.rm = TRUE)
intervals_grouped_per_fly <- intervals_grouped %>% summarise_all(c("median","sd"), na.rm = TRUE)




MAD_intervals <- dplyr::select(intervals_grouped, 4,7:11)
# MADs <- summarise_all(MAD_intervals, funs(mad), na.rm = TRUE)
MADs <- MAD_intervals %>% summarise_all(c("mad"), na.rm = TRUE)


names(MADs) <- c("i", "j", "Index", "MAD_SI", "MAD_HP", "MAD_DI")




merged_transients <- left_join(grouped_per_fly, intervals_grouped_per_fly, by = c("i" = "i", "j" = "j", "index" = "Index"))
merged_transients$relaxtime <- merged_transients$HP_median - merged_transients$peaked_median

merged_transients <- left_join(merged_transients, MADs, by = c("i" = "i", "j" = "j", "index" = "Index"))

# Summarize per kymograph
# grouped <- group_by(grouped, CODE, jpg, i) # Add Code as grouping variable


# Collapse all transients into a per fly average
#grouped_per_fly <- summarise_all(grouped, funs(mean,sd), na.rm = TRUE)
#grouped_per_fly <- grouped_per_fly[,c(1:4,6:28)]



# Read all data from previous analysis
all_data <- as_tibble(read.csv(file = "all_data_table.csv", sep = ","))
all_data <- group_by(all_data, file)
# all_data <- all_data[-which(is.na(all_data$Xpos)),] # Bugfix


# Combine all data and timestamps
final_all_data <- left_join(all_data, merged_transients, by = c("file" = "jpg"))
final_all_data <- unique(final_all_data)

# add Age:
final_all_data$Age <- as.factor(unlist(lapply(final_all_data$file, function(x) substring(x,gregexpr('_',x)[[1]][2]+1,gregexpr('_',x)[[1]]+1)[2])))
# add Sex:
final_all_data$Sex <- as.factor(unlist(lapply(final_all_data$file, function(x) substring(x,gregexpr('_',x)[[1]][3]-1,gregexpr('_',x)[[1]]-1)[3])))

# Find all cxd files with duplicated peaks (bug from previous runs) and keep only the first
final_all_data_distinct <- final_all_data %>% dplyr::select(cxd_file, Xpos, file) %>% distinct()
final_all_data <- final_all_data[match(final_all_data_distinct$file, final_all_data$file),]

# Remove Peak-Number from filename
final_all_data$file <- as.factor(unlist(lapply(final_all_data$file, function(x) paste0(substring(x, 1, gregexpr('peak_',x)[[1]][1]+4), "X",substring(x,gregexpr('_at',x)[[1]][1])))))

# Column clean-up
final_all_data <- dplyr::select(final_all_data, -c("CODE.x", "sizeX", "sizeY", "sizeZ", "sizeC", "sizeT", "pixelType", "bitsPerPixel" , "imageCount" , "dimensionOrder", 
                                                   "orderCertain", "rgb", "littleEndian", "interleaved", "falseColor", "metadataComplete", "thumbnail", "series" , "resolutionLevel", 
                                                   "created_unix_from_file", "CODE.y", "recordingday", "week_year" , "i", "j", "index", "peaked_sd", "tt10p_sd", "tt25p_sd",
                                                   "tt50p_sd", "tt75p_sd", "tt90p_sd", "tt90r_sd", "tt75r_sd", "tt50r_sd", "tt25r_sd", "tt10r_sd", "starts_median" , "peaks_median", 
                                                   "ends_median", "starts_sd", "peaks_sd","ends_sd" , "SI_sd", "deltaSI_sd", "deltaSIplus_sd", "HP_sd", "DI_sd" )) # remove superfluous columns

final_all_data <- dplyr::select(final_all_data, 23, 9, 1, 12, 13, 20:22, 10, 11, 15, 44, 45, 3:6, 38:40, 35, 24, 7:8, 41:43, 25:34, 36:37, 2, 14, 16:19, everything())
final_all_data <- final_all_data[,1:42]
final_all_data$stocks <- crosses$Stock.collection[match(final_all_data$CODE, crosses$CODE)]
final_all_data <- unique(final_all_data)
final_all_data <- group_by(final_all_data, flyID)

  write.csv(final_all_data, file = "final_all_data.csv", row.names = F)


# Collapse all transients into a per genotype average
data_per_fly <- final_all_data[,c(3,12:21, 23:31,34,36:39)]
data_per_fly <- group_by(data_per_fly, flyID)

final_all_data_per_fly <- data_per_fly %>% summarise_all(c("median"), na.rm = TRUE)
final_all_data_per_fly <-left_join(final_all_data_per_fly, dplyr::select(final_all_data, c(1:10, 22, 35, 41:43)), by = c("flyID" = "flyID"))
final_all_data_per_fly <- unique(final_all_data_per_fly)

  write.csv(final_all_data_per_fly, file = "final_all_data_per_fly.csv", row.names = F)

}





