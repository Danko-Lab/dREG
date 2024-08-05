install.packages("data.table", repos = "http://cran.us.r-project.org")
install.packages("e1071", repos = "http://cran.us.r-project.org")
install.packages("mvtnorm", repos = "http://cran.us.r-project.org")
install.packages("parallel", repos = "http://cran.us.r-project.org")
install.packages("randomForest", repos = "http://cran.us.r-project.org")
install.packages("rmutil", repos = "http://cran.us.r-project.org")
#install.packages("rphast", repos = "http://cran.us.r-project.org")
install.packages("snowfall", repos = "http://cran.us.r-project.org")
install.packages("devtools", repos = "http://cran.us.r-project.org")
library(devtools)
# Install RPHAST via github since the CRAN package has been pulled.
devtools::install_github("CshlSiepelLab/RPHAST")
# Here a bug in current devtools 1.12.0 if the installation starts from a sub-folder.
r <- try( install_github("andrelmartins/bigWig", subdir="bigWig", force=TRUE) )
if ( class( r )== "try-error" )
	stop("!!!!!Failed to install bigWig from github using install_github method.!!!!\n 
	You can install manually using the following commands:\n
	---------------------------------------------------\n
	git clone https://github.com/andrelmartins/bigWig.git\n
	cd bigWig \n
	R CMD INSTALL bigWig\n
	---------------------------------------------------");



