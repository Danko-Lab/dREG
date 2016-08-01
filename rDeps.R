install.packages("boot", repos = "http://cran.us.r-project.org")
install.packages("grid", repos = "http://cran.us.r-project.org")
install.packages("e1071", repos = "http://cran.us.r-project.org")
install.packages("rphast", repos = "http://cran.us.r-project.org")
install.packages("devtools", repos = "http://cran.us.r-project.org")
library(devtools)
# Here a bug in current devtools 1.12.0 if the installation starts from a sub-folder.
if (class(try( install_github("andrelmartins/bigWig", subdir="bigWig", force=TRUE) ) )== "try-error")
	stop("!!!!!Failed to install bigWig from github using install_github method.!!!!\n 
	You can install manually as following commands\n
	---------------------------------------------------\n
	git clone https://github.com/andrelmartins/bigWig.git\n
	cd bigWig
	R CMD INSTALL bigWig\n
	---------------------------------------------------\n
	");



